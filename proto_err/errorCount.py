#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time
from utils import *
import difflib
import itertools
# import re

class error():
    """ Information about the errors in a read """
    def __init__(self,true,emission,read,readPos):
        self.true = true
        self.emission = emission
        self.read = read
        self.readPos = readPos # position on read where error starts 
        ## Was the read aligned correctly? 
        try:
            ## messy way of extracting read length
            s = list(self.read.qname.split('st=')[1])
            curInt = '0'
            tempList = []
            while curInt.isdigit():
                curInt = s.pop(0)
                tempList.append(curInt)
            sampPos = int("".join(tempList[:-1]))
            self.alignedDist =  sampPos - self.read.positions[0]
        except:
            self.alignedDist = None
        if self.alignedDist is None:
            self.alignedCorrectly = None
        elif not self.alignedDist is None and (self.alignedDist < len(self.read.seq)):
            self.alignedCorrectly = True
        else:
            self.alignedCorrectly = False

    def __str__(self):
        return "%s error(%s to %s)" %(self.errorType,self.true,self.emission)

    def before(self,j):
        """Return the preceding j bases,return N when bases missing"""
        i = self.readPos
        b = self.read.seq[i-j:i]
        while len(b) < j:
            b = 'N' +b
        return b
    def after(self,j):
        """Return the following j bases,return N when bases missing"""
        i = self.readPos
        a = self.read.seq[i+1:j+i+1]
        while len(a) < j:
            a = a + 'N'
        return a

    def qscore(self,i):
        """Return the quality score at a base +i i from the error start position"""
        return asciiToInt(self.read.qqual[self.readPos+1])

    @property 
    def trueSeq(self):
        """Sequence of the truth"""
        return  str(self.read.seq)
    @property 
    def emissionSeq(self):
        """Sequence emmited"""
        trSeq = list(self.read.seq)
        trSeq[self.readPos]  = self.emission
        return  "".join(trSeq)

    @property 
    def errorType(self):
        if self.isSnp:
            return 'SNP'
        elif self.isInsertion:
            return 'Insertion'
        elif self.isDeletion:
            return 'Deletion'
        else:
            return 'Unknown Error Type'
    @property 
    def isSnp(self):
        """Is the error a SNP"""
        return len(self.true) == len(self.emission)
    @property 
    def isIndel(self):
        """Is the error an INDEL"""
        return len(self.true) != len(self.emission)
    @property 
    def isInsertion(self):
        """Is the error an insertion"""
        return len(self.true) < len(self.emission)
    @property 
    def isDeletion(self):
        """Is the error a deletion"""
        return len(self.true) > len(self.emission)
    @property
    def qual(self):
        """Quality score of error base"""
        return asciiToInt(self.read.qqual[self.readPos])


class errorReader():
    """Iterable over errors in aligned reads"""
    def __init__(self, samfile,ref):
        self.samfile = pysam.Samfile( samfile ).fetch()
        self.ref = ref
        self.alphabet = getAlphabet()
        self.readCounter = {}
        self.readCounter['Total'] = 0
        self.readCounter['UnMapped'] = 0
        self.readCounter['Mapped'] = 0
        self.readCounter['perfectAlignments'] = 0
        self.readCounter['mismatchedAlignments'] = 0
        self.errorList = []

    def __iter__(self):
        return self
    @property
    def refRead(self):
        """
        Function to return the reference sequence a read is aligned to
        """
        refRead = [str(self.ref[pos]) for pos in self.read.positions]
        return "".join(refRead)

    def readNext(self):
        ## If there are errors left in the read 
        self.read = self.samfile.next()
        self.readCounter['Total'] += 1
        if self.read.is_unmapped:
            self.readCounter['UnMapped'] += 1
        else:
            self.readCounter['Mapped'] += 1
            self.checkRead()

    def next(self):
        ## if the error list is empty get the next read
        while not self.errorList:
            self.readNext()
        else:
            return self.errorList.pop(0)

        # self.readCounter['totalErrorBases'] = len(self.errorList)

    def checkRead(self):
        """
        A function which take a read and generates some error objects
        """
        if self.read.opt('NM') == 0:
            self.readCounter['perfectAlignments'] += 1
        else:
            self.readCounter['mismatchedAlignments'] += 1
        self.currentRefReadList = list(self.refRead)
        self.currentReadList = list(str(self.read.seq))
        self.readPos = 0
        for tup in self.read.cigar:
            cigarInt = tup[0]
            numBases = tup[1]
            if cigarInt == 0:
                ## Match or mismatch
                self.checkSNPs(N=numBases)
                self.readPos += numBases
            elif cigarInt == 1:
                ## Insertion to the reference
                self.checkInsertion(N=numBases)
                self.readPos += numBases
            elif cigarInt == 2:
                ## Deletion from the reference
                self.checkDeletion(N=numBases)
            elif cigarInt == 3:
                ## skipped region from the reference
                self.checkSkipped(N=numBases)
            elif cigarInt == 4:
                ##  soft clipping (clipped sequences present in SEQ)
                self.checkSoftClipped(N=numBases)
            elif cigarInt == 5:
                ##  hard clipping (clipped sequences NOT present in SEQ)
                self.checkHardClipped(N=numBases)
            elif cigarInt == 6:
                ## padding (silent deletion from padded reference)
                self.checkPadding(N=numBases)
            elif cigarInt == 7:
                ## sequence match
                self.checkSeqMatch(N=numBases)
            elif cigarInt == 8:
                ## sequence mismatch
                self.checkSeqMismatch(N=numBases)
        assert self.currentRefReadList == []
        assert self.currentReadList == []
    def checkSNPs(self,N):
        readSeg = popLong(self.currentReadList,0,N)
        refSeg = popLong(self.currentRefReadList,0,N)
        assert len(readSeg) == len(refSeg)
        for i,tupl in enumerate(itertools.izip(refSeg,readSeg)):
            true,emission = tupl
            ## if it's an error, check the preceding  opt.maxOrder bases
            if not true == emission:
                ## Check preceding bases
                self.errorList.append(error(true,emission,self.read,readPos=i+self.readPos))   
    def checkInsertion(self,N):
        insSeg = popLong(self.currentReadList,0,N)
        self.errorList.append(error(true='',emission="".join(insSeg),read=self.read,readPos=self.readPos))
    def checkDeletion(self,N):
        i = self.read.positions[self.readPos-1] + 1
        j =  self.read.positions[self.readPos]        
        delSeg = self.ref[i:j]
        self.errorList.append(error(true=delSeg,emission="",read=self.read,readPos=self.readPos))
    def checkSkipped(self,N):
        logging.error("Haven't written a handler for this case yet")
        0/0
    def checkSoftClipped(self,N):
        del self.currentReadList[0:N]
    def checkHardClipped(self,N):
        pass
    def checkPadding(self,N):
        logging.error("Haven't written a handler for this case yet")
        0/0
    def checkSeqMismatch(self,N):
        logging.error("Haven't written a handler for this case yet")
        0/0
    def checkSkipped(self,N):
        logging.error("Haven't written a handler for this case yet")
        0/0





class counter(): 
    """Takes a list of errors and does some kmer counting"""
    def __init__(self,ref,errorList=None,samfile=None,opt=None):
        self.logger     = logging.getLogger()
        if (not errorList and not samfile) or (errorList and samfile):
            self.logger.error("counter takes errorList or samfile, at least one and not both")
        if samfile and not errorList:
            self.errorList = []
            for error in errorReader(samfile,ref):
                self.errorList.append(error)
        else:
            self.errorList = errorList
        # the results dictionary
        self.res = {}
        self.res['kmerCounts'] = AutoVivification()
        self.res['RefCounts'] =  AutoVivification()
        self.res['errorMode'] = AutoVivification()
        self.res['qualCounter'] = AutoVivification()
        # set toplevel name 
        self.ref = ref
        self.numErrors = len(self.errorList)
        self.countErrorKmerRun = False
        if opt:
            self.opt = opt


    def setup(self,opt):
        """
        Function to set up dictionary structure for outputed stats
        """
        self.opt = opt
        ## Create the dictonary of error modes
        ## res[A][T] = counts of A -> T SNPs
        # self.res['errorMode'] = dict(zip(getAlphabet(),[dict(zip(getAlphabet(),[0]*4)) for i in getAlphabet()]))

    def kmerFreq(self,seq,kmer):
        """
        Function to calculate the number of times a kmer appears in a sequence
        """
        count = seq.count(kmer)
        return count,float(count)/len(seq)

    def countRefKmer(self,maxKmerLength=None):
        """
        Function to count all kmers in long sequence (reference) of length kmerLen or below
        """
        if not maxKmerLength:
            maxKmerLength = self.opt.maxKmerLength
        alpabet = getAlphabet()
        for klen in [i+1 for i in range(self.opt.maxKmerLength)]:
            # generate all kmers of length klen
            kmerList = kmerCombo(alpabet,klen)
            self.res['RefCounts'] = dict(zip(kmerList,[0]*len(kmerList)))
            for kmer in kmerList:
                self.res['RefCounts'][str(kmer)] = self.ref.count(kmer)
        return self.res['RefCounts']

    def countErrorKmer(self,maxKmerLength=None):
        """
        Function which takes a list of errors and counts the kmers before  
        and after and error
        """
        self.countErrorKmerRun = True
        if not maxKmerLength:
            maxKmerLength = self.opt.maxKmerLength
        # 
        for error in self.errorList:
            try:
                self.res['errorMode'][error.true][error.emission] += 1
                self.res['qualCounter'][error.qual] +=1
            except:
                self.res['errorMode'][error.true][error.emission] = 1
                self.res['qualCounter'][error.qual] =1
            for j in [k +1 for k in range(self.opt.maxKmerLength)]:
                try:
                    self.res['kmerCounts']['before'][error.true][error.emission][error.before(j)] += 1
                except:
                    self.res['kmerCounts']['before'][error.true][error.emission][error.before(j)] = 1
                try:
                    self.res['kmerCounts']['after'][error.true][error.emission][error.after(j)] += 1
                except:
                    self.res['kmerCounts']['after'][error.true][error.emission][error.after(j)]= 1
        return self.res['errorMode'],self.res['kmerCounts']

    def readDiff(self,read,ref):
        return difflib.ndiff(read,ref)

    def countAlignedBases(self,read):
        count = 0
        for t in read.cigar:
            count += t[1]
        return count

    def getCount(self,truth=None,emission=None,kmer='',after=False):
        """Gets the count for a given {truth,emmision,kmer}"""
        ## This function is a bit hacky, come back and rewrite later.
        ## If flexible querys are priorty, maybe create mongo db of error objects?
        if not self.countErrorKmerRun:
            self.countErrorKmer()
        if after:
            dic = self.res['kmerCounts']['after']
        else:
            dic = self.res['kmerCounts']['before']

        if truth and emission and kmer:
            ##For particular transistion and kmer
            return dic[truth][emission][kmer]
        elif truth is None and emission is None:
            ## We have to iterate through everything 
            countOut = 0 
            for emmitedDic in dic.values():
                for kmerDic in emmitedDic.values():
                    if kmer:
                        ## If only one kmer
                        countOut += kmerDic.get(kmer,0)
                    else:
                        ## Otherwise count them all
                        for count in kmerDic.values():
                            countOut += count

            return countOut
        elif not truth is None and emission is None:
            tmpDic = dic[truth]
            countOut = 0 
            for kmerDic in tmpDic.values():
                    if kmer:
                        ## If only one kmer
                        countOut += kmerDic.get(kmer,0)
                    else:
                        ## Otherwise count them all
                        for count in kmerDic.values():
                            countOut += count
            if countOut == {}:
                return 0
            else:
                return countOut
        elif truth is None and not emission is None:
            countOut = 0 
            for emmitedDic in dic.values():
                kmerDic =  emmitedDic[emission]
                if kmer:
                    ## If only one kmer
                    countOut += kmerDic.get(kmer,0)
                else:
                    ## Otherwise count them all
                    for count in kmerDic.values():
                        countOut += count

            return countOut






           


class aggregator(): 
    """Takes a counter object and does some aggreation"""
    def __init__(self,counter):
        self.counts = counter.res
        self.precedingKmers = AutoVivification()

    def countPrecedingKmers(self):
        """Count all the kmers preceding any error"""
        for truth, emmitedDic in self.counts['kmerCounts']['before'].iteritems():
            for emmited,kmerDic in emmitedDic.iteritems():
                for kmer,count in kmerDic.iteritems():
                    try:
                        self.precedingKmers[kmer] += count
                    except:
                        self.precedingKmers[kmer] = count
    def precedingKmersCount(self,kmer):
        if not len(self.precedingKmers):
            self.countPrecedingKmers()
        return self.precedingKmers[kmer]









