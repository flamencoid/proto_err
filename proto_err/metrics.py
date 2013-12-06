#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time
from utils import *
import difflib
import itertools

class error():
    """ Information about the errors in a read """
    def __init__(self,true,emission,read,readPos):
        self.true = true
        self.emission = emission
        self.read = read
        self.readPos = readPos # position on read where error starts 
        # self.leftFlank = leftFlank 
        # self.rightFlank = rightFlank
    def __str__(self):
        return "%s error(%s to %s)" %(self.errorType,self.true,self.emission)

# read.seq[i-self.opt.maxOrder:i],read.seq[i+1:self.opt.maxOrder+i+1]

    def before(self,j):
        """Return the preceding j bases,return N when bases missing"""
        # b = self.leftFlank[-j:
        i = self.readPos
        b = self.read.seq[i-j:i]
        while len(b) < j:
            b = 'N' +b
        return b
    def after(self,j):
        """Return the following j bases,return N when bases missing"""
        # a = self.rightFlank[0:j]
        i = self.readPos
        a = self.read.seq[i+1:j+i+1]
        while len(a) < j:
            a = a + 'N'
        return a

    @property 
    def trueSeq(self):
        """Sequence of the truth"""
        return  self.leftFlank + self.true + self.rightFlank
    @property 
    def emissionSeq(self):
        """Sequence emmited"""
        return  self.leftFlank + self.emission + self.rightFlank

    @property 
    def errorType(self):
        if self.isSnp:
            return 'SNP'
        elif self.isIndel:
            return 'INDEL'
    @property 
    def isSnp(self):
        """Is the error a SNP"""
        return len(self.true) == len(self.emission)
    @property 
    def isIndel(self):
        """Is the error an INDEL"""
        return len(self.true) != len(self.emission)

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

    def getRefRead(self,positions):
        """
        Function to return the reference sequence a read is aligned to
        """
        return self.ref[positions[0]:positions[-1]+1]

    def readNext(self):
        ## If there are errors left in the read 
        self.read = self.samfile.next()
        self.readCounter['Total'] += 1
        if self.read.is_unmapped:
            self.readCounter['UnMapped'] += 1
        else:
            self.readCounter['Mapped'] += 1
            # print read.qname
            self.checkRead(self.read)

    def next(self):
        ## if the error list is empty get the next read
        while not self.errorList:
            self.readNext()
        else:
            return self.errorList.pop(0)

        # self.readCounter['totalErrorBases'] = len(self.errorList)

    def checkRead(self,read):
        """
        A function which take a read and generates some error objects
        """
        NM = read.opt('NM')
        if NM == 0:
            self.readCounter['perfectAlignments'] += 1
            # for base in self.alphabet:
            #     self.res['errorMode'][base][base] += read.seq.count(base)
        else:
            self.readCounter['mismatchedAlignments'] += 1
            ## Check what type of mismatch it is.
            ## if the cigarsting only has M then only SNP errors
            if (read.rlen == read.cigar[0][1]) and (read.cigar[0][0] == 0):
                refRead = list(str(self.getRefRead(read.positions)))
                for i,tupl in enumerate(itertools.izip_longest(refRead,list(str(read.seq)))):
                    true,emission = tupl
                    
                    ## if it's an error, check the preceding  opt.maxOrder bases
                    if not true == emission:
                        ## Check preceding bases
                        self.errorList.append(error(true,emission,read,readPos=i))

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
            except:
                self.res['errorMode'][error.true][error.emission] = 1
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









