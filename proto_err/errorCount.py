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
    """ 
    Information about the errors in a read 

    Attributes
    ----------
    true : string
        reference bases
    emission : string 
        Read bases
    read : AlignedRead
        pysam AlignedRead object from which error was derived
    readPos : int
        position along read where error occured
    readPer : float
        position along read where error occured / length of read sequence
    alignedDist : int or None
        number of bases between position of start of alignment and where read was sampled None if sampled position n/a 
    alignedCorrectly : bool or None
        If the alignedDist is less then the read length True else False None is alignedDist is None

    Parameters
    ----------
    true: string
        reference base(s)

    emission: string
        read base(s)

    read: AlignedRead
         pysam AlignedRead object from which error was derived

    readPos: int
        position along read where error occured / length of read sequence



    See Also
    --------
    counter

    Examples
    --------
    >>> from pysam import AlignedRead
    >>> from import proto_err.errorCount import error
    >>> read = AlignedRead()
    >>> read.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    >>> read.qual = '++))++)+*)******)))+)**+*+++)**)*+)'
    >>> error = proto_err.errorCount.error('T','A',read,0)
    >>> error.isSnp
    True
    >>> error.isIndel
    False
    >>> error.qual
    11
    >>> error.before(2)
    NN 
    >>> error.after(2)
    GC 
    >>> error.qscore(0) ## equivalent to error.qual
    11
    >>> error.qscore(2)
    9

    """
    

    def __init__(self,true,emission,read,readPos):
        self.true = true
        self.emission = emission
        self.read = read
        self.readPos = readPos # position on read where error starts 
        self.readPer = float(readPos) / float(len(read.seq))
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
        """Return the quality score at a base +i i from the error start position
        error.qscore(0) is equivalent to error.qual """
        return asciiToInt(self.read.qqual[self.readPos+i])

    @property 
    def errorType(self):
        """The type of read error SNP, Insertion or Deletion"""
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
        """Return whether or not the error is a SNP"""
        return len(self.true) == len(self.emission)
    @property 
    def isIndel(self):
        """Return whether or not the error is a INDEL"""
        return len(self.true) != len(self.emission)
    @property 
    def isInsertion(self):
        """Return whether or not the error is an insertion"""
        return len(self.true) < len(self.emission)
    @property 
    def isDeletion(self):
        """Return whether or not the error is a deletion"""
        return len(self.true) > len(self.emission)
    @property
    def qual(self):
        """Return quality score of the base where the error occured"""
        return asciiToInt(self.read.qqual[self.readPos])


class errorReader():
    """ 
    Iterable over errors in aligned reads

    Attributes
    ----------
    readCounter : dict
        A dictonary of counts of read alignment

    Parameters
    ----------
    samfile: string
        path to samfile

    ref: string
        reference sequence


    See Also
    --------
    error, counter

    Examples
    --------
    >>> from proto_err.errorCount import errorReader
    >>> reader = errorReader(opt.samfile,ref)
    >>> for error in reader:
    >>>     print error
    SNP error(T to C)
    SNP error(A to T)
    Deletion error(CGCA to )
    ...

    >>> print reader.readCounter
    {'UnMapped': 454, 'mismatchedAlignments': 617, 'Total': 1225, 'perfectAlignments': 154, 'Mapped': 771}

    """

    def __init__(self, samfile,ref):
        self.__samfile = pysam.Samfile( samfile ).fetch()
        self.__ref = ref
        self.__alphabet = getAlphabet()
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
    def __refRead(self):
        """
        Get the reference sequence the read is aligned to
        """
        refRead = [str(self.__ref[pos]) for pos in self.__read.positions]
        return "".join(refRead)

    def __readNext(self):
        """
        Iterates to the next read in samfile
        """
        ## If there are errors left in the read 
        self.__read = self.__samfile.next()
        self.readCounter['Total'] += 1
        if self.__read.is_unmapped:
            self.readCounter['UnMapped'] += 1
        else:
            self.readCounter['Mapped'] += 1
            self.__checkRead()

    def next(self):
        """
        Returns the next error in the samfile aligned read
        """
        ## if the error list is empty get the next read
        while not self.errorList:
            self.__readNext()
        else:
            return self.errorList.pop(0)

        # self.readCounter['totalErrorBases'] = len(self.errorList)

    def __checkRead(self):
        """
        Checks current read for errors
        """
        if self.__read.opt('NM') == 0:
            self.readCounter['perfectAlignments'] += 1
        else:
            self.readCounter['mismatchedAlignments'] += 1
        self.__currentRefReadList = list(self.__refRead)
        self.__currentReadList = list(str(self.__read.seq))
        self.__readPos = 0
        for tup in self.__read.cigar:
            cigarInt = tup[0]
            numBases = tup[1]
            if cigarInt == 0:
                ## Match or mismatch
                self.__checkSNPs(N=numBases)
                self.__readPos += numBases
            elif cigarInt == 1:
                ## Insertion to the reference
                self.__checkInsertion(N=numBases)
                self.__readPos += numBases
            elif cigarInt == 2:
                ## Deletion from the reference
                self.__checkDeletion(N=numBases)
            elif cigarInt == 3:
                ## skipped region from the reference
                self.__checkSkipped(N=numBases)
            elif cigarInt == 4:
                ##  soft clipping (clipped sequences present in SEQ)
                self.__checkSoftClipped(N=numBases)
            elif cigarInt == 5:
                ##  hard clipping (clipped sequences NOT present in SEQ)
                self.__checkHardClipped(N=numBases)
            elif cigarInt == 6:
                ## padding (silent deletion from padded reference)
                self.__checkPadding(N=numBases)
            elif cigarInt == 7:
                ## sequence match
                self.__checkSeqMatch(N=numBases)
            elif cigarInt == 8:
                ## sequence mismatch
                self.__checkSeqMismatch(N=numBases)
        assert self.__currentRefReadList == []
        assert self.__currentReadList == []
    def __checkSNPs(self,N):
        """
        Checks read segment for SNP errors. called when cigarstring = M:N 
        """

        readSeg = popLong(self.__currentReadList,0,N)
        refSeg = popLong(self.__currentRefReadList,0,N)
        assert len(readSeg) == len(refSeg)
        for i,tupl in enumerate(itertools.izip(refSeg,readSeg)):
            true,emission = tupl
            ## if it's an error, check the preceding  opt.maxOrder bases
            if not true == emission:
                ## Check preceding bases
                self.errorList.append(error(true,emission,self.__read,readPos=i+self.__readPos))   
    def __checkInsertion(self,N):
        """
        Checks Insertion read segment for errors. called when cigarstring = I:N 
        """
        insSeg = popLong(self.__currentReadList,0,N)
        self.errorList.append(error(true='',emission="".join(insSeg),read=self.__read,readPos=self.__readPos))
    def __checkDeletion(self,N):
        """
        Checks Deletionread segment for  errors. called when cigarstring = D:N 
        """
        i = self.__read.positions[self.__readPos-1] + 1
        j =  self.__read.positions[self.__readPos]        
        delSeg = str(self.__ref[i:j])
        self.errorList.append(error(true=delSeg,emission="",read=self.__read,readPos=self.__readPos))
    def __checkSkipped(self,N):
        """
        Checks skipped read segment for errors. called when cigarstring = N:N 
        """
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSoftClipped(self,N):
        """
        Checks SoftClipped read segment for errors. called when cigarstring = S:N 
        """
        del self.__currentReadList[0:N]
    def __checkHardClipped(self,N):
        """
        Checks HardClipped read segment for errors. called when cigarstring = H:N 
        """
        pass
    def __checkPadding(self,N):
        """
        Checks Padding read segment for errors. called when cigarstring = P:N 
        """
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSeqMismatch(self,N):
        """
        Checks SeqMismatch read segment for errors. called when cigarstring = =:N 
        """
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSkipped(self,N):
        """
        Checks Skipped read segment for errors. called when cigarstring = X:N 
        """
        logging.error("Haven't written a handler for this case yet")
        0/0





class counter(): 
    """ 
    Takes a list of errors or samfile and does some kmer counting

    Attributes
    ----------
    ref : string
        The reference sequence
    errorList : list 
        A list of error objects (Optional, can take samfile instead)
    samfile : string
        Path to samfile (Optional, can take a list of errors instead)
    opt : dict
        Options passed from OptionParser

    Parameters
    ----------
    errorList : list
        A list of error objects on which counting is done 
    res : dict
        A dictonary containg resulting counts
    numErrors : int
        Number of errors in errorList


    See Also
    --------
    error, errorReader

    Examples
    --------
    >>> errorCounter = counter(ref,samfile=opt.samfile)
    >>> errorCounter.setup(opt)
    >>> errorCounter.countRefKmer()
    """

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
        self.__countErrorKmerRun = False
        if opt:
            self.setup(opt)


    def setup(self,opt):
        """
        Pass options from OptionParser
        """
        self.opt = opt

    def __kmerFreq(self,seq,kmer):
        """
        Calculate the number of times a kmer appears in a sequence
        seq : SeqRecord Object
        kmer : kmer string
        """
        count = seq.count(kmer)
        return count,float(count)/len(seq)

    def countRefKmer(self,maxKmerLength=None):
        """
        Count all kmers in the reference of length kmerLen or below


        Parameters
        ----------
        maxKmerLength : int
            This is a type.
            Maximum kmer length


        Returns
        ----------
        dict
             dictonary of counts in the form

            {'AAT':123,'AA':123,...}

            emmissionDic[trueBase][emmitedBase] = count of trueBase -> emmitedBase transition

            This also stored in the counter object in counter.res['RefCounts'] 

         See Also
        --------
        countErrorKmer

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> errorCounter.countRefKmer(maxKmerLength = 3) 
            >>> ## counts all possible kmers of length 1,2,3 in reference
            >>> {'AAT':123,'AA':123,...}
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
        Count all kmers in the list of errors of length kmerLen or below before 
        and after and error.


        Parameters
        ----------
        maxKmerLength : int
            This is a type.
            Maximum kmer length


        Returns
        ----------
        emmissionDic,beforeAfterCount
            dict,dict
            emmissionDic counts the number of times a transition occurs

            e.g. {'A':{'T':123,'C':123,...,},'T':{'A':123,...},...}

            emmissionDic[trueBase][emmitedBase] = count of trueBase -> emmitedBase transition



            beforeAfterCount counts the number of times a kmer appears before or after an error

            beforeAfterCount['before'][trueBase][emmitedBase][kmer] = count of kmer occurance before trueBase -> emmitedBase transition

            beforeAfterCount['after'][trueBase][emmitedBase][kmer] = count of kmer occurance after trueBase -> emmitedBase transition 



            These dictonaries are also stored in the counter object in counter.res['errorMode'],counter.res['kmerCounts']

         See Also
        --------
        countRefKmer

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> errorCounter.countErrorKmer(maxKmerLength = 3) 
            >>> ## counts all possible kmers of length 1,2,3 before and after an error. 
            >>> ({'A': {'C': 113, 'T': 105, 'G': 84}, ...}, {'after':{'A': {'C': {'CTT': 2, 'GCA': 1,...},..}} ,'before:{'A': {'C': {'CTT': 2, 'GCA': 1,...},..}}'})
        """
        self.__countErrorKmerRun = True
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

    def __readDiff(self,read,ref):
        return difflib.ndiff(read,ref)

    def countAlignedBases(self):
        """
        Counts the number of aligned bases in counter

        Returns
        ----------
        int
            Count of aligned bases

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> errorCounter.countAlignedBases() 
            >>> 12335
        """
        count = 0
        for e in self.errorList:
            read = e.read
            for t in read.cigar:
                count += t[1]
        return count

    def getCount(self,truth=None,emission=None,kmer='',after=False):
        """
        Gets the count for a given {truth,emmision,kmer}


        Parameters
        ----------
        truth : string
            truth base(s)
        emission : string
            emmited base(s)
        kmer : string
            preceding or following kmer
        after : bool
            Default False. If True kmer is before, if False kmer is after

        Returns
        ----------
        int
            Count of errors matching Parameter query

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> #### Count all the errors    
            >>>print errorCounter.getCount()
            1234
            >>> ####   Count all the errors preceded by kmer 'A'
            >>> print errorCounter.getCount(kmer='AA')
            123
            >>> ####  Count all the errors followed by kmer 'AA'
            >>> print errorCounter.getCount(kmer='AA',after=True)
            12
            >>> ####     Count all the errors for truth 'A' preceded by an A
            >>> print errorCounter.getCount(truth='A',kmer='A')
            123
            >>> ####     Count all the errors A->T preceded by A
            >>> print errorCounter.getCount(truth='A',emission='T', kmer='A')
            12
            >>> ####  Count all the errors where the emmited base is 'A' preceded by any kmer
            >>> print errorCounter.getCount(emission='A')
            2
        """

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






           


# class aggregator(): 
#     """Takes a counter object and does some aggreation"""
#     def __init__(self,counter):
#         self.counts = counter.res
#         self.precedingKmers = AutoVivification()

#     def countPrecedingKmers(self):
#         """Count all the kmers preceding any error"""
#         for truth, emmitedDic in self.counts['kmerCounts']['before'].iteritems():
#             for emmited,kmerDic in emmitedDic.iteritems():
#                 for kmer,count in kmerDic.iteritems():
#                     try:
#                         self.precedingKmers[kmer] += count
#                     except:
#                         self.precedingKmers[kmer] = count
#     def precedingKmersCount(self,kmer):
#         if not len(self.precedingKmers):
#             self.countPrecedingKmers()
#         return self.precedingKmers[kmer]









