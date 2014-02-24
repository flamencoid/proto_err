#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time
from utils import *
import difflib
import itertools
from query import errordb
from plot import *
from error import error
import os
from collections import Counter as listCounter
import pandas as pd
from pprint import pprint
import scipy.stats

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

import re

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
        self.readCounter['TotalReads'] = 0
        self.readCounter['UnMapped'] = 0
        self.readCounter['Mapped'] = 0
        self.readCounter['perfectAlignments'] = 0
        self.readCounter['mismatchedAlignments'] = 0
        self.readCounter['totalAlignedBases'] = 0
        self.readCounter['totalBases'] = 0
        self.readCounter['HardClippedReads'] = 0
        for c in list('MIDNSHP=X'):
            self.readCounter[c] = 0
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
        
        self.readCounter['TotalReads'] += 1
        self.readCounter['totalBases'] += self.__read.rlen
        if self.__read.is_unmapped:
            self.readCounter['UnMapped'] += 1
        else:
            self.readCounter['Mapped'] += 1
            self.readCounter['totalAlignedBases'] += self.__read.alen
            if 'H' in list(self.__read.cigarstring):
                self.readCounter['HardClippedReads'] += 1
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
        self.__readPosIndex = 0
        self.refPos = self.__read.positions[0]
        for tup in self.__read.cigar:
            cigarInt = tup[0]
            numBases = tup[1]
            if cigarInt == 0:
                ## Match or mismatch
                self.__checkSNPs(N=numBases)
            elif cigarInt == 1:                
                ## Insertion to the reference
                self.__checkInsertion(N=numBases)
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

        self.readCounter['M'] += N
        readSeg = popLong(self.__currentReadList,0,N)
        refSeg = popLong(self.__currentRefReadList,0,N)
        assert len(readSeg) == len(refSeg)
        for i,tupl in enumerate(itertools.izip(refSeg,readSeg)):
            true,emission = tupl
            ## if it's an error, check the preceding  opt.maxOrder bases
            if not true == emission:
                ## Check preceding bases
                if self.__read.is_reverse:
                    ## Need to take the reverse complement of the truth and emission to get the error which occured in the read
                    tmpTrue = str(SeqRecord(Seq(true),'','','').reverse_complement().seq)
                    tmpEmission = str(SeqRecord(Seq(emission),'','','').reverse_complement().seq)
                    self.errorList.append(error(tmpTrue,tmpEmission,self.__read,readPos=self.__read.rlen - (i+self.__readPos) - 1,refPos=self.refPos+i))

                else:
                    self.errorList.append(error(true,emission,self.__read,readPos=i+self.__readPos,refPos=self.refPos+i)) 
        self.__readPos += N
        self.__readPosIndex += N  
        self.refPos += N
    def __checkInsertion(self,N):
        """
        Checks Insertion read segment for errors. called when cigarstring = I:N 
        """
        self.readCounter['I'] += N
        insSeg = popLong(self.__currentReadList,0,N)
        if self.__read.is_reverse:
            tmpinsSeg = str(SeqRecord(Seq("".join(insSeg)),'','','').reverse_complement().seq)
            self.errorList.append(error(true='',emission=tmpinsSeg,read=self.__read,readPos=self.__read.rlen - self.__readPos -1,refPos=self.refPos))
        else:
            self.errorList.append(error(true='',emission="".join(insSeg),read=self.__read,readPos=self.__readPos,refPos=self.refPos)) # -1 as insertion occurs at previous base (if simulating)
        self.__readPos += N

    def __checkDeletion(self,N):
        """
        Checks Deletion read segment for  errors. called when cigarstring = D:N 
        """
        self.readCounter['D'] += N
        i = self.__read.positions[self.__readPosIndex-1] + 1
        j =  self.__read.positions[self.__readPosIndex]        
        delSeg = str(self.__ref[i:j])
        if self.__read.is_reverse:
            tmpdelSeg = str(SeqRecord(Seq(delSeg),'','','').reverse_complement().seq)
            self.errorList.append(error(true=tmpdelSeg,emission="",read=self.__read,readPos=self.__read.rlen - self.__readPos -1,refPos=self.refPos))
        else:
            self.errorList.append(error(true=delSeg,emission="",read=self.__read,readPos=self.__readPos,refPos=self.refPos))
        self.refPos += N
        

    def __checkSkipped(self,N):
        """
        Checks skipped read segment for errors. called when cigarstring = N:N 
        """
        self.readCounter['N'] += N
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSoftClipped(self,N):
        """
        Checks SoftClipped read segment for errors. called when cigarstring = S:N 
        """
        self.readCounter['S'] += N
        del self.__currentReadList[0:N]
    def __checkHardClipped(self,N):
        """
        Checks HardClipped read segment for errors. called when cigarstring = H:N 
        """
        self.readCounter['H'] += N
        logging.warning("We shouldn't have HardClipped bases in Samfile")


        
    def __checkPadding(self,N):
        """
        Checks Padding read segment for errors. called when cigarstring = P:N 
        """
        self.readCounter['P'] += N
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSeqMismatch(self,N):
        """
        Checks SeqMismatch read segment for errors. called when cigarstring = =:N 
        """
        self.readCounter['='] += N
        logging.error("Haven't written a handler for this case yet")
        0/0
    def __checkSkipped(self,N):
        """
        Checks Skipped read segment for errors. called when cigarstring = X:N 
        """
        self.readCounter['X'] += N
        logging.error("Haven't written a handler for this case yet")
        0/0





class counter(): 
    """ 
    Takes a list of errors or samfile and does some kmer counting

    Parameters
    ----------
    ref : string
        The reference sequence
    errorList : list 
        A list of error objects (Optional, can take samfile instead)
    samfile : string
        Path to samfile (Optional, can take a list of errors instead)
    opt : dict
        Options passed from OptionParser

    Attributes
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

    def __init__(self,ref,opt,samfile,makeDB=False):
        self.logger     = logging.getLogger()
        self.samfile = samfile

        self.errorList = []
        reader = errorReader(samfile,ref)
        logging.info("Parsing samfile for errors")
        
        self.readCounter = reader.readCounter

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
        self.setup(opt)
        ## Connection to mongoDB, runs without for now.
        self.errordb = {}
        self.errordb['errors'] = errordb(database=opt.dbName,collection=opt.observedErrorDBName)
        self.errordb['simulatedErrors'] = errordb(database=opt.dbName,collection=opt.simulatedErrorDBName)
        self.errordb['metaData'] = errordb(database=opt.dbName,collection='metaData')

        for error in reader:
            self.errorList.append(error)
        logging.info("Found %i errors in samfile" % (len(self.errorList)))
        if makeDB:
            self.errordb['errors'].deleteAll()
            self.errorList = self.errordb['errors'].addErrors(self.errorList)
        else:
            # logging.info("Downloading pre-generated errors from database")
            # self.errorList = self.errordb['errors'].find_errors()
            self.logger.warning("""### Using pre-generated database.""")
            self.logger.warning("""### Initiate counter with makeDB=True to wipe and repopulate""")

        self.qscoresCount = None ## So we don't do counting twice
        self.kmerSamCount = None
        self.probKmerRefDic = {}
        self.freqKmerDic = {}
        self.alignedStringQual = ""
        self.alignedString = ""
        self.kmerQualCount = {}
        


        

        



    def setup(self,opt):
        """
        Pass options from OptionParser
        """
        self.opt = opt
        ## Check that all required opts exist
        if not self.opt.maxKmerLength:
            self.logger.info("opt.maxKmerLength not set, Defaulting to 3")
            self.opt.maxKmerLength = 3

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
        for klen in [i+1 for i in range(self.opt.maxKmerLength)]:
            # generate all kmers of length klen
            kmerList = kmerCombo(klen)
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
    def plotHist(self):
        """
        Plots histograms


        Parameters
        ----------
        none : none
            no desc


        Returns
        ----------
        outputs a png file to opt.imgDir

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> errorCounter.plotHist()
        """
        # if not self.__countErrorKmerRun:
        #     self.countErrorKmer()
        # tmpDic =  self.res['errorMode']
        observedDic = {}
        expectedDic = {}
        simulatedDic = {}
        multi = AutoVivification()
        for t in getAlphabet():
            for e in getAlphabet():
                if t != e:
                    observed = self.getCount(truth=t,emission=e)
                    expected = self.getExpectedCount(truth=t,emission=e)
                    simulated = self.getSimulatedCount(truth=t,emission=e)
                    observedDic[t + '->' + e] = observed
                    expectedDic[t + '->' + e] = expected
                    simulatedDic[t + '->' + e] = simulated
                    multi[t + '->' + e]['Expected'] = expected
                    multi[t + '->' + e]['Observed'] = observed
                    multi[t + '->' + e]['Simulated'] = simulated
        histPlotter(dic=observedDic,opt=self.opt,filename="SNP_observed_transition").plot()
        histPlotter(dic=expectedDic,opt=self.opt,filename="SNP_expected_transition").plot()
        histPlotter(dic=simulatedDic,opt=self.opt,filename="SNP_simulated_transition").plot()
        # multiHistPlotter(dic=multi,opt=self.opt,filename="SNP_observed_vs_expected_transition").plot()
        multiHistPlotter(dic=multi,opt=self.opt,filename="SNP_observed_simulated_expected_transition").plot()
        # # Count deletion kmers
        # # get maximum deletion length
        observedSize = [d['tlen'] for d in self.errordb['errors'].find( {'type' : 'Deletion'}, {'tlen':1} )]
        simulatedSize = [d['tlen'] for d in self.errordb['simulatedErrors'].find( {'type' : 'Deletion'}, {'tlen':1} )]
        if observedSize and simulatedSize:
            densityPlotterFromLists(dic={'observed':observedSize,'simulated':simulatedSize},
                                opt=self.opt,filename="deletion_size_hist_dens").plot(geom='dens')
            densityPlotterFromLists(dic={'observed':observedSize,'simulated':simulatedSize},
                                opt=self.opt,filename="deletion_size_bar").plot(geom='bar')

            maxLen = max(observedSize + simulatedSize) 
            # for order in [i+1 for i in range(maxLen)]:
            for order in [i+1 for i in range(5)]:
                dic = self.__countToDic(self.getCount(type='Deletion',tlenRange=order,returnList=True)[1],attribute='true')
                histPlotter(dic=dic,opt=self.opt,filename="deletedKmerCount/deleted_kmer_observed_order_%i" % (order)).plot()
                dic = self.__countToDic(self.getSimulatedCount(type='Deletion',tlenRange=order,returnList=True)[1],attribute='true')
                histPlotter(dic=dic,opt=self.opt,filename="deletedKmerCount/deleted_kmer_simulated_order_%i" % (order)).plot()

        ## Count insterted kmers 
        observedSize = [d['tlen'] for d in self.errordb['errors'].find( {'type' : 'Insertion'}, {'tlen':1} )]
        simulatedSize = [d['tlen'] for d in self.errordb['simulatedErrors'].find( {'type' : 'Insertion'}, {'tlen':1} )]
        if observedSize and simulatedSize:
            densityPlotterFromLists(dic={'observed':observedSize,'simulated':simulatedSize},
                                opt=self.opt,filename="insertion_size_hist_dens").plot(geom='dens')
            densityPlotterFromLists(dic={'observed':observedSize,'simulated':simulatedSize},
                                opt=self.opt,filename="insertion_size_bar").plot(geom='bar')

            maxLen = max(observedSize + simulatedSize) 
            # for order in [i+1 for i in range(maxLen)]:
            for order in [i+1 for i in range(5)]:
                dic = self.__countToDic(self.getCount(type='Insertion',tlenRange=order,returnList=True)[1],attribute='emmision')
                histPlotter(dic=dic,opt=self.opt,filename="insertedKmerCount/inserted_kmer_observed_order_%i" % (order)).plot()
                dic = self.__countToDic(self.getSimulatedCount(type='Insertion',tlenRange=order,returnList=True)[1],attribute='emmision')
                histPlotter(dic=dic,opt=self.opt,filename="insertedKmerCount/inserted_kmer_simulated_order_%i" % (order)).plot()


        ## Count kmers      
        for order in [i+1 for i in range(self.opt.maxKmerLength)]:
            dicBefore = AutoVivification()
            dicAfter = AutoVivification()
            for errorType in ['SNP','Insertion','Deletion']:
                for kmer in kmerCombo(r=order):
                    count = self.getCount(kmerBefore = kmer,type=errorType)
                    if count != 0 :
                        dicBefore[errorType][kmer] = count
                    count = self.getCount(kmerAfter = kmer,type=errorType)
                    if count != 0:
                        dicAfter[errorType][kmer] = count 
                histPlotter(dic=dicBefore[errorType],opt=self.opt,filename="kmerBeforeCount/kmer_before_%s_order_%s"%(errorType,order)).plot()
                histPlotter(dic=dicBefore[errorType],opt=self.opt,filename="kmerAfterCount/kmer_after_%s_order_%s"%(errorType,order)).plot() 

    def __countToDic(self,errorList,attribute):
        """Convience funciton to help with converting counts to dic for plotting

        Parameters
        ----------
        errorList : list
            list of error objects 

        Returns
        ----------
        dic
            dictonary of counts
        """
        dic = {}
        for error in errorList:
            if attribute == 'true':
                try:
                    dic[error.true] += 1
                except:
                    dic[error.true] = 1
            elif attribute == 'emmision':
                try:
                    dic[error.emission] += 1
                except:
                    dic[error.emission] = 1                
        return dic

    def getFreqQualGivenKmer(self,kmer,qual):
        """Get the frequency of a particular quality value for a context in the samfile"""
        # if not self.qscoresCount:
        #     logging.info("Counting qual score frequency")
        #     self.qscoresCount = listCounter([asciiToInt(i) for i in list("".join(read.qual for read in pysam.Samfile( self.samfile ).fetch()))])
        qscoresCount =  self.getKmerQualCount(kmer,qual)    
        return float(qscoresCount) / float(sum(self.kmerQualCount[kmer].values()))

    def getFreqKmer(self,kmer):
        """Get the frequency of a particular kmer in the samfile"""
        try:
            p = self.freqKmerDic[kmer]
        except:
            count = 0
            for read in pysam.Samfile( self.samfile ).fetch():
                count += str(read.seq).count(kmer)
            p = float(count) / float(self.readCounter['totalAlignedBases'])
            self.freqKmerDic[kmer] = p

        return p

    def getKmerQualCount(self,kmer,qual=1000):
        """Get the count of a particular kmer, quality kmer context
        i.e. """
        try:
            kmerQualCount = self.kmerQualCount[kmer][qual]
        except:
            logging.info("Calculating mean quality score for all contexts of length %i" % len(kmer))
            qualDic = AutoVivification()
            for read in samReader(samfile=self.samfile,ref=self.ref):
                ## Prepopulated with trimers
                kmerList =  kmerCombo(len(kmer))
                posAdd = int(float(len(kmer) - 1) / 2.0 )
                patternDic = dict(zip(kmerList,[re.compile(kmer) for kmer in kmerList]))
                for kmer,pattern in patternDic.iteritems():
                    ## Should this be greping in the reference rather than the read?
                    # results = pattern.finditer("".join(read.seq))
                    results = pattern.finditer("".join(read.refRead))
                    for result in results:
                        try:
                            qualDic[kmer].append(asciiToInt(read.qual[result.start(0) + posAdd]))
                        except:
                            qualDic[kmer] = [asciiToInt(read.qual[result.start(0) + posAdd])]
            for kmer in kmerCombo(len(kmer)):
                self.kmerQualCount[kmer] = listCounter(qualDic[kmer])
                try:
                    # Remove any counts which lie in deletions
                    del self.kmerQualCount[kmer]['_']
                except:
                    pass
            kmerQualCount = self.kmerQualCount[kmer][qual]    
        return  kmerQualCount       

    def getContextMeanQualScore(self,kmer):
        """Get the mean quality score associated with a given context (3mer only atm) returns the qual of the middle base"""
        assert len(kmer) in [1,3,5] ## Too lazy to code for all odd numbers
        try:
            ## If the count already exists in the dict use that
            qual = self.kmerQualCount[kmer].keys()
            counts = self.kmerQualCount[kmer].values()
        except KeyError:
            ## If not, calculated all the counts of that kmer length 
            ## (to save iterating through whole samfile for every count)
            self.getKmerQualCount(kmer)
            qual = self.kmerQualCount[kmer].keys()
            counts = self.kmerQualCount[kmer].values()
        sumQual = 0
        for qual,count in zip(qual,counts):
            sumQual+= (qual * count)

        meanQual = float(sumQual) / sum(counts)
        return meanQual
        



    def getExpectedCount(self,truth='',emission='',kmerBefore='',kmerAfter='',qual=None,type='SNP'):
        """
        Gets the expected count of a transition. 

        E[true,emmision] = (frequency of truth in reference) * (Observed SNP Count) * (Probability of emmision | SNP ). 

        i.e for a random reference E[true,emission] = 1/4 * numSnps * 1/3 

        Parameters
        ----------
        truth : string
            truth base(s)
        emission : string
            emmited base(s)

        Returns
        ----------
        int
            Expected count

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> print errorCounter.getExpectedCount(truth='A',emission='T')
            100.36
            >>> print errorCounter.getCount(truth='A',emission='T')
            111
        """


        # simulationMetaData =  self.errordb['metaData'].find_one({'type':'simulation'},{'snpFreq':1,'readMean':1,'numReads':1,'SnpIndelRatio':1})
        if qual is not None:
            context = kmerBefore+truth+kmerAfter
            pContext = self.probKmerRef(context) * self.getFreqQualGivenKmer(kmer=context,qual=qual)
            # print context,qual, self.getFreqQualGivenKmer(kmer=context,qual=qual)
        elif truth or emission or kmerBefore or kmerAfter:
            ## Use the frequency of the kmer in the reference or in the samfile???? 
            ## Has to be in the reference, as observed contexts are after errors
            pContext = self.probKmerRef(kmerBefore+truth+kmerAfter)

            # pContext = self.getFreqKmer(kmerBefore+truth+kmerAfter)
            # print self.probKmerRef(kmerBefore+truth+kmerAfter) , pContext
        else:
            pContext = 1
        # print simulationMetaData
        # ExpectedSNPCount = self.readCounter['totalAlignedBases'] * simulationMetaData['snpFreq'] * simulationMetaData['SnpIndelRatio'] #snpFreq is actually the errorFrequencey
        # ExpectedSNPCount = self.readCounter['totalBases'] * simulationMetaData['snpFreq'] * simulationMetaData['SnpIndelRatio'] #snpFreq is actually the errorFrequencey
        # probEmmission = float(1)/float(3)
        if qual is not None:
            probError = qscoreToProb(qual)
        else:
            # probSNPError = simulationMetaData['snpFreq'] * simulationMetaData['SnpIndelRatio']
            ## Using the observed snp error rate
            # probSNPError = float(self.getCount(type='SNP')) / float(self.readCounter['totalAlignedBases']) # I don't think this is the correct way to do this
            if type == "SNP"  and kmerBefore+truth+kmerAfter:
                ## If you're looking at a particular qual score
                probError = qscoreToProb(self.getContextMeanQualScore(kmer=kmerBefore+truth+kmerAfter))
            else:
                ## If you want the total expectation
                probError = float(self.getCount(type=type)) / float(self.readCounter['totalAlignedBases'])
        if emission and type == 'SNP':
            expectedCount = (self.readCounter['totalAlignedBases'] * pContext * probError) / 3.0 # assume equal probabilites of eny transition
        elif type == 'SNP' and not emission:
            expectedCount = pContext * probError * self.readCounter['totalAlignedBases']
        elif not type == 'SNP':
            expectedCount = (self.readCounter['totalAlignedBases'] * pContext * probError)
        else:
            raise(ValueError)



        return round(expectedCount,2)


    def getSimulatedCount(self,truth=None,emission=None,kmerBefore=None,kmerAfter=None,
                type=None,maxAlignedDist=None,readLength=None,readPosRange=[],readPerRange=[],
                qualRange=[],tlenRange=[],returnList=False):
        """
        Gets the count of simulated errors.

        Parameters
        ----------
        truth : string
            truth base(s)
        emission : string
            emmited base(s)
        kmerBefore : string
            preceding kmer
        kmerAfter : string
            following kmer
        type : string
            type of error SNP Insertion or Deletion
        maxAlignedDist : int
            maximum aligned distance e.g. maxAlignedDist = 10 returns all errors where the read mapped within 10 bp of where it was sampled from 
        readPosRange : list or tuple of ints
            list or tuple of two elements range of error position along read. e.g. [0,10] returns all errors in the first 10bp of read
        readPerRange : list or tuple of ints
            list or tuple of two elements range of error percentage along read. e.g. [0,10] returns all errors in the first 10percent of read
        qualRange : list or tuple of ints
            list or tuple of two elements range of quality of error e.g. [10,10] returns all errors with quality score 10
        returnList : bool
            default False. If true returns list of error documents as an additional argument

        Returns
        ----------
        returnList = False int
            Count of errors matching Parameter query
        returnList = True int,list
            Count of errors matching Parameter query, list of error documents

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> #### Count all the errors    
            >>>print errorCounter.getSimulatedCount()
            1234
            >>> ####   Count all the errors preceded by kmer 'A'
            >>> print errorCounter.getSimulatedCount(kmerBefore='AA')
            123
            >>> ####  Count all the errors followed by kmer 'AA'
            >>> print errorCounter.getSimulatedCount(kmerAfter='AA')
            12
            >>> ####     Count all the errors for truth 'A' preceded by an A
            >>> print errorCounter.getSimulatedCount(truth='A',kmerBefore'A')
            123
            >>> ####     Count all the errors A->T preceded by A
            >>> print errorCounter.getSimulatedCount(truth='A',emission='T', kmerBefore='A')
            12
            >>> ####  Count all the errors where the emmited base is 'A' preceded by any kmer
            >>> print errorCounter.getSimulatedCount(emission='A')
            2
        """
        return self.getCount(truth=truth,emission=emission,kmerBefore=kmerBefore,
                            kmerAfter=kmerAfter,type=type,maxAlignedDist=maxAlignedDist,
                            readPosRange=readPosRange,readPerRange=readPerRange,
                            qualRange=qualRange,tlenRange=tlenRange,readLength=readLength,returnList=returnList,
                            collection='simulatedErrors')


    def probKmerRef(self,kmer):
        """
        Gets the fraction of 

        Parameters
        ----------
        kmer : string
            kmer to count

        Returns
        ----------
        float
            freq in ref
        """
        try:
            p = self.probKmerRefDic[kmer]
        except:
            p = float(self.ref.count(kmer)) / len(self.ref)
            self.probKmerRefDic[kmer] = p
        return p


    def compareSimulationToResults(self,query={}):
        errorList = self.errordb['errors'].find_errors(query=query)
        simErrorList = self.errordb['simulatedErrors'].find_errors(query=query)
        errorDict  = {}
        for err in errorList:
            errorDict[err.refPos] = err
        simErrorDict = {}
        for err in simErrorList:
            simErrorDict[err.refPos] = err

        TP = [] ## Corretly returned error
        FP = [] ## Found but not in simulation
        FN = [] ## In simulationg but not found
        for refPos,err in errorDict.iteritems():
            if refPos in simErrorDict:
                simErr = simErrorDict.pop(refPos)
                if err == simErr:
                    TP.append(err)
                else:
                    FP.append(err)
                    FN.append(simErr)
            else:
                FP.append(err)
        ## All remaining sim errors are FN
        FN.extend(simErrorDict.values())
        errorDic = {'TP':TP,'FP':FP,'FN':FN}
        return errorDic

    def summary(self):
        "Log a summary of errors found and samfile counts"

        ## Initialise some pandas dataframes to sort output for writing to csv. 
        SNPRow = [self.getCount(type='SNP'),self.getSimulatedCount(type='SNP')]
        INSRow = [self.getCount(type='Insertion'),self.getSimulatedCount(type='Insertion')]
        DELRow = [self.getCount(type='Deletion'),self.getSimulatedCount(type='Deletion')]
        


        for key,value in self.readCounter.iteritems():
            self.logger.info("### Count of %s in samfile = %i" % (key,value))
        self.logger.info("### Total SNP errors observed = %i" % (SNPRow[0]) )
        self.logger.info("### Total SNP errors simulated = %i" % (SNPRow[1]))
        self.logger.info("### Total Insertion errors observed = %i" % (INSRow[0]))
        self.logger.info("### Total Insertion errors simulated = %i" % (INSRow[1]))
        self.logger.info("### Total Deletion errors observed = %i" % (DELRow[0]))
        self.logger.info("### Total Deletion errors simulated = %i" % (DELRow[1]))

        ## How many errors are from mismapped reads?
        self.logger.info("### Total errors from reads mapped mapped correctly =%i" % (self.getCount(mappedCorrectly=1)))
        self.logger.info("### Total errors from mismapped reads =%i" % (self.getCount(mappedCorrectly=0)))
        self.logger.info("### Percentage of errors from  mismapped reads =%f " % (round(100*(float(self.getCount(mappedCorrectly=0)) / float(self.getCount()) ),2)))

        ## How many errors were correctly recovered?

        ## Download errors and simulated errors and generate two dicts of the form
        ## {redPos : errorObject}
        compareDict = self.compareSimulationToResults()
        self.logger.info("### Total simulated errors correctly found (TP) = %i" % (len(compareDict['TP'])))
        self.logger.info("### Total errors found which were not simulated (FP)= %i" % (len(compareDict['FP'])))
        self.logger.info("### Total errors simulated which were not found (FN)= %i" % (len(compareDict['FN'])))
        compareDict = self.compareSimulationToResults(query={'type':'SNP'})
        SNPRow.extend([len(compareDict['TP']),len(compareDict['FP']),len(compareDict['FN']) ])
        self.logger.info("### Total simulated SNPs correctly found (TP) = %i" % (len(compareDict['TP'])))
        self.logger.info("### Total SNPs found which were not simulated (FP)= %i" % (len(compareDict['FP'])))
        self.logger.info("### Total SNPs simulated which were not found (FN)= %i" % (len(compareDict['FN'])))
        compareDict = self.compareSimulationToResults(query={'type':'Insertion'})
        INSRow.extend([len(compareDict['TP']),len(compareDict['FP']),len(compareDict['FN']) ])
        self.logger.info("### Total simulated Insertions correctly found (TP) = %i" % (len(compareDict['TP'])))
        self.logger.info("### Total Insertions found which were not simulated (FP)= %i" % (len(compareDict['FP'])))
        self.logger.info("### Total Insertions simulated which were not found (FN)= %i" % (len(compareDict['FN'])))
        compareDict = self.compareSimulationToResults(query={'type':'Deletion'})
        DELRow.extend([len(compareDict['TP']),len(compareDict['FP']),len(compareDict['FN']) ])
        self.logger.info("### Total simulated Deletions correctly found (TP) = %i" % (len(compareDict['TP'])))
        self.logger.info("### Total Deletions found which were not simulated (FP)= %i" % (len(compareDict['FP'])))
        self.logger.info("### Total Deletions simulated which were not found (FN)= %i" % (len(compareDict['FN'])))

        



        TotalRow = []
        for i in range(len(SNPRow)):
            TotalRow.append(sum([SNPRow[i],INSRow[i],DELRow[i]]))

        try:
            SNPRow.append(round(float(SNPRow[2]) / float(SNPRow[2] + SNPRow[3]),2))
        except:
            SNPRow.append('-')
        try:
            INSRow.append(round(float(INSRow[2]) / float(INSRow[2] + INSRow[3]),2))
        except:
           INSRow.append('-')
        try:
            DELRow.append(round(float(DELRow[2]) / float(DELRow[2] + DELRow[3]),2))
        except:
            DELRow.append('-')
        try:
            TotalRow.append(round(float(TotalRow[2]) / float(TotalRow[2] + TotalRow[3]),2))
        except:
            TotalRow.append('-')
        try:
            SNPRow.append(round(float(SNPRow[2]) / float(SNPRow[2] + SNPRow[4]),2))
        except:
            SNPRow.append('-')
        try:
            INSRow.append(round(float(INSRow[2]) / float(INSRow[2] + INSRow[4]),2))
        except:
            INSRow.append('-')
        try:
            DELRow.append(round(float(DELRow[2]) / float(DELRow[2] + DELRow[4]),2))
        except:
            DELRow.append('-')
        try:
            TotalRow.append(round(float(TotalRow[2]) / float(TotalRow[2] + TotalRow[4]),2))
        except:
            TotalRow.append('-')
        outfile = self.opt.outDir + 'errorsCount.dat'
        logging.info('Writing error counts to %s' % (outfile))
        self.errorsDF = pd.DataFrame.from_items([('Total', TotalRow),('SNP', SNPRow), ('INS', INSRow),('DEL', DELRow)],
                                orient='index', columns=['samCount', 'simCount', 'TP','FP','FN','Precision','Recall'])
        self.errorsDF.to_csv(path_or_buf=outfile,sep='\t')

        outfile = self.opt.outDir + 'readCounts.dat'
        logging.info('Writing read counts to %s' % (outfile))
        self.readCounts = pd.DataFrame.from_items([('Counts', self.readCounter.values())],
                                orient='index', columns=self.readCounter.keys())
        self.readCounts.to_csv(path_or_buf=outfile,sep='\t')

    def DELTransitionStats(self,kmerLength=3):
        """
        Generate some stats about the context in which insertions occur 
        We want to know if certain kmer contexts occur more frequently before insetions then would by chance
        """
        logging.info("Generating Context bias statistics for insertions")
        alphabet = ['A','T','C','G']
        outdic = AutoVivification()
        ## Just do trimers before for now.
        totalExpectedCount = self.getExpectedCount(type="Deletion")
        cnt = 0
        for context in kmerCombo(kmerLength):
            cnt += 1
            logging.info("Calculating stats for context %s : %i of %i" %(context,cnt,4**kmerLength))
            ## Stats for contexts without q scores
            outdic['DEL'][context]['samCount'] = self.getCount(kmerBefore=context,type='Deletion')
            outdic['DEL'][context]['simCount'] = self.getSimulatedCount(kmerBefore=context,type='Deletion')
            outdic['DEL'][context]['expectedCount'] = self.getExpectedCount(kmerBefore=context,type='Deletion')
            outdic['DEL'][context]['expectedOccurancyofContext'] = int(round(self.probKmerRef(context) * self.readCounter['totalAlignedBases']))
            outdic['DEL'][context]['observedOccurancyofContext'] = int(round(self.getFreqKmer(context) * self.readCounter['totalAlignedBases']))
            outdic['DEL'][context]['avgQual'] = round(self.getContextMeanQualScore(kmer=context),2)
            # print outdic['DEL'][context]['samCount'] , self.getCount(type='Deletion'), outdic['DEL'][context]['expectedCount'],totalExpectedCount
            outdic['DEL'][context]['pvalue'] = scipy.stats.binom_test(outdic['DEL'][context]['samCount'], 
                                                                    self.getCount(type='Deletion'), 
                                                                    outdic['DEL'][context]['expectedCount']/totalExpectedCount)
        print self.getCount(type='Deletion'),self.getSimulatedCount(type='Deletion'),totalExpectedCount

        logging.info("Generating readable output for Deletion bias Stats")
        outputList = []
        outputListHeader = ['avgQual','expectedOccurancyofContext','observedOccurancyofContext','simCount','samCount','expectedCount','pvalue']
        for context in kmerCombo(kmerLength):
            row = [context] + [outdic['DEL'][context][t] for t in outputListHeader]
            outputList.append(row)
        self.contextDELStats = pd.DataFrame(outputList, columns=['ContextTrue']+outputListHeader)

        ## Adjust the p-values for multiple testing
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(list(self.contextDELStats['pvalue'])), method = 'BH')
        self.contextDELStats['pvalue-adjust'] = p_adjust
        self.contextDELStats = self.contextDELStats.sort(['pvalue'],ascending=True)
        outfile = self.opt.outDir + 'contextDELStats_kmer%i.dat' % kmerLength
        logging.info("Writing context bias stats to %s" % (outfile))
        self.contextDELStats.to_csv(path_or_buf=outfile,sep='\t') 

    def INSTransitionStats(self,kmerLength=3):
        """
        Generate some stats about the context in which insertions occur 
        We want to know if certain kmer contexts occur more frequently before insetions then would by chance
        """
        logging.info("Generating Context bias statistics for insertions")
        alphabet = ['A','T','C','G']
        outdic = AutoVivification()
        ## Just do trimers before for now.
        totalExpectedCount = self.getExpectedCount(type="Insertion")
        cnt = 0
        for context in kmerCombo(kmerLength):
            cnt += 1
            logging.info("Calculating stats for context %s : %i of %i" %(context,cnt,4**kmerLength))
            ## Stats for contexts without q scores
            outdic['INS'][context]['samCount'] = self.getCount(kmerBefore=context,type='Insertion')
            outdic['INS'][context]['simCount'] = self.getSimulatedCount(kmerBefore=context,type='Insertion')
            outdic['INS'][context]['expectedCount'] = self.getExpectedCount(kmerBefore=context,type='Insertion')
            outdic['INS'][context]['expectedOccurancyofContext'] = int(round(self.probKmerRef(context) * self.readCounter['totalAlignedBases']))
            outdic['INS'][context]['observedOccurancyofContext'] = int(round(self.getFreqKmer(context) * self.readCounter['totalAlignedBases']))
            outdic['INS'][context]['avgQual'] = round(self.getContextMeanQualScore(kmer=context),2)
            # print outdic['INS'][context]['samCount'] , self.getCount(type='Insertion'), outdic['INS'][context]['expectedCount'],totalExpectedCount
            outdic['INS'][context]['pvalue'] = scipy.stats.binom_test(outdic['INS'][context]['samCount'], 
                                                                    self.getCount(type='Insertion'), 
                                                                    outdic['INS'][context]['expectedCount']/totalExpectedCount)
        print self.getCount(type='Insertion'),self.getSimulatedCount(type='Insertion'),totalExpectedCount

        logging.info("Generating readable output for Insertion bias Stats")
        outputList = []
        outputListHeader = ['avgQual','expectedOccurancyofContext','observedOccurancyofContext','simCount','samCount','expectedCount','pvalue']
        for context in kmerCombo(kmerLength):
            row = [context] + [outdic['INS'][context][t] for t in outputListHeader]
            outputList.append(row)
        self.contextINSStats = pd.DataFrame(outputList, columns=['ContextTrue']+outputListHeader)

        ## Adjust the p-values for multiple testing
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(list(self.contextINSStats['pvalue'])), method = 'BH')
        self.contextINSStats['pvalue-adjust'] = p_adjust
        self.contextINSStats = self.contextINSStats.sort(['pvalue'],ascending=True)
        outfile = self.opt.outDir + 'contextINSStats_kmer%i.dat' % kmerLength
        logging.info("Writing context bias stats to %s" % (outfile))
        self.contextINSStats.to_csv(path_or_buf=outfile,sep='\t') 

    def SNPTransitionStats(self):
        """
        generate some statistics of transition probabilites
        """
        logging.info("Generating Context bias statistics")
        alphabet = ['A','T','C','G']

        logging.info("Generating unique list of quality scores")
        qscores = unique([d['qual'] for d in self.errordb['errors'].find({'qual':{'$exists':1},'type':'SNP' } )])
        outdic = AutoVivification()
        outDicContextOnly = AutoVivification()
        samCount = 0
        simCount = 0 
        expectedCount = 0
        totalExpectedCount = self.getExpectedCount()
        cnt = 0
        for before in alphabet:
            for truth in alphabet:
                for after in alphabet:
                    context = before + truth + after
                    cnt += 1
                    logging.info("Calculating stats for context %s : %i of 64" %(context,cnt))

                    for emission in alphabet:
                        if not truth == emission:
                            ## Stats for contexts without q scores
                            outDicContextOnly['SNP'][context]['samCount'][emission] = self.getCount(truth=truth,emission=emission,
                                                                                    kmerBefore=before,kmerAfter=after)
                            outDicContextOnly['SNP'][context]['simCount'][emission] = self.getSimulatedCount(truth=truth,emission=emission,
                                                                                    kmerBefore=before,kmerAfter=after)
                            outDicContextOnly['SNP'][context]['expectedCount'][emission] = self.getExpectedCount(truth=truth,emission=emission,
                                                                                    kmerBefore=before,kmerAfter=after)

                            outDicContextOnly['SNP'][context]['expectedOccurancyofContext'][emission] = int(round(self.probKmerRef(context) * self.readCounter['totalAlignedBases']))
                            outDicContextOnly['SNP'][context]['observedOccurancyofContext'][emission] = int(round(self.getFreqKmer(context) * self.readCounter['totalAlignedBases']))
                            outDicContextOnly['SNP'][context]['avgQual'][emission] = round(self.getContextMeanQualScore(kmer=context),2)

                            outDicContextOnly['SNP'][context]['pvalue'][emission] = scipy.stats.binom_test(outDicContextOnly['SNP'][context]['samCount'][emission], self.getCount(type='SNP'), outDicContextOnly['SNP'][context]['expectedCount'][emission]/totalExpectedCount)
                            if (outDicContextOnly['SNP'][context]['samCount'][emission]) > 0:
                                outDicContextOnly['SNP'][context]['expectedDiff'][emission] = round((float( outDicContextOnly['SNP'][context]['samCount'][emission]) - float(outDicContextOnly['SNP'][context]['expectedCount'][emission])) ,2)
                            ## Stats for contexts with q scores
                            
                            for qscore in qscores:
                                outdic['SNP'][context][qscore]['samCount'][emission] = self.getCount(truth=truth,emission=emission,
                                                                                kmerBefore=before,kmerAfter=after,qualRange=[qscore,qscore])
                                outdic['SNP'][context][qscore]['simCount'][emission] = self.getSimulatedCount(truth=truth,emission=emission,
                                                                                kmerBefore=before,kmerAfter=after,qualRange=[qscore,qscore])
                                outdic['SNP'][context][qscore]['expectedCount'][emission] = self.getExpectedCount(truth=truth,emission=emission,
                                                                                kmerBefore=before,kmerAfter=after,qual=qscore)
                                outdic['SNP'][context][qscore]['pvalue'][emission] = float(scipy.stats.binom_test(outdic['SNP'][context][qscore]['samCount'][emission], totalExpectedCount, outdic['SNP'][context][qscore]['expectedCount'][emission]/totalExpectedCount) )
                                
                                samCount +=  outdic['SNP'][context][qscore]['samCount'][emission]
                                simCount += outdic['SNP'][context][qscore]['simCount'][emission]
                                expectedCount += outdic['SNP'][context][qscore]['expectedCount'][emission]
                                # if (outdic['SNP'][context][qscore]['samCount'][emission]) > 0:
                                    # outdic['SNP'][context][qscore]['expectedDiff'][emission] = round((float( outdic['SNP'][context][qscore]['samCount'][emission]) - float(outdic['SNP'][context][qscore]['expectedCount'][emission])) ,2)
                    logging.info("p-values for context %s %s" %(context,outDicContextOnly['SNP'][context]['pvalue'].values()))

        print samCount,simCount,expectedCount
        print self.getCount(type='SNP'),self.getSimulatedCount(type='SNP'),totalExpectedCount

        ## Create friendly output

        ##                      simCount     ExptCount   SamCount    pvalue  pvalue-corrected
        ## Context ContextOut
        ## AAA  ATC 123 122     124 0.05    0.2
        ## ...
        logging.info("Generating readable output")
        outputList = []
        outputWithQscoresList = []
        outputListHeader = ['avgQual','expectedOccurancyofContext','observedOccurancyofContext','simCount','samCount','expectedCount','pvalue']
        for before in alphabet:
            for truth in alphabet:
                for after in alphabet:
                    for emission in alphabet:
                        if not truth == emission:
                            context = before + truth + after
                            contextOut = before + emission + after
                            row = [context,contextOut] + [outDicContextOnly['SNP'][context][t][emission] for t in outputListHeader]
                            outputList.append(row)
                            for qscore in qscores:
                                row = [context,contextOut,qscore] + [outdic['SNP'][context][qscore][t][emission] for t in ['simCount','samCount','expectedCount','pvalue']]
                                outputWithQscoresList.append(row)
        self.contextStats = pd.DataFrame(outputList, columns=['ContextTrue','ContextEmit']+outputListHeader)
        self.contextQualStats = pd.DataFrame(outputWithQscoresList, columns=['ContextTrue','ContextEmit','qscore','simCount','samCount','expectedCount','pvalue'])

        ## Adjust the p-values for multiple testing
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(list(self.contextStats['pvalue'])), method = 'BH')
        p_adjustQ = stats.p_adjust(FloatVector(list(self.contextQualStats['pvalue'])), method = 'BH')
        
        self.contextStats['pvalue-adjust'] = p_adjust
        self.contextQualStats['pvalue-adjust'] = p_adjustQ



        self.contextStats = self.contextStats.sort(['pvalue'],ascending=True)
        self.contextQualStats = self.contextQualStats.sort(['pvalue'],ascending=True)
        
        outfile = self.opt.outDir + 'contextStats.dat'
        logging.info("Writing context bias stats to %s" % (outfile))
        self.contextStats.to_csv(path_or_buf=outfile,sep='\t') 

        
        outfile = self.opt.outDir + 'contextQualityScoreStats.dat'
        logging.info("Writing context and qual score bias stats to %s" % (outfile))
        self.contextQualStats.to_csv(path_or_buf=outfile,sep='\t') 

                        
    









    def getCount(self,truth=None,emission=None,kmerBefore=None,kmerAfter=None,
                type=None,maxAlignedDist=None,mappedCorrectly=None,readLength=None,readPosRange=[],readPerRange=[],
                qualRange=[],tlenRange=[],returnList=False,collection='errors'):
        """
        Gets the count for a given {truth,emmision,kmer}


        Parameters
        ----------
        truth : string
            truth base(s)
        emission : string
            emmited base(s)
        kmerBefore : string
            preceding kmer
        kmerAfter : string
            following kmer
        type : string
            type of error SNP Insertion or Deletion
        maxAlignedDist : int
            maximum aligned distance e.g. maxAlignedDist = 10 returns all errors where the read mapped within 10 bp of where it was sampled from 
        mappedCorrectly : bool
            Mapped distance within read length
        readPosRange : list or tuple of ints
            list or tuple of two elements range of error position along read. e.g. [0,10] returns all errors in the first 10bp of read
        readPerRange : list or tuple of ints
            list or tuple of two elements range of error percentage along read. e.g. [0,10] returns all errors in the first 10percent of read
        qualRange : list or tuple of ints
            list or tuple of two elements range of quality of error e.g. [10,10] returns all errors with quality score 10
        tlenRange : int or list or tuple
            int list or tuple of two elements range of quality of error e.g. [10,10] returns all errors with quality score 10 equivalent to 10
        returnList : bool
            default False. If true returns list of error documents as an additional argument


        Returns
        ----------
        returnList = False int
            Count of errors matching Parameter query
        returnList = True int,list
            Count of errors matching Parameter query, list of error documents

        Examples
        --------
            >>> from errorCount import counter
            >>> errorCounter = counter(ref,samfile)
            >>> #### Count all the errors    
            >>>print errorCounter.getCount()
            1234
            >>> ####   Count all the errors preceded by kmer 'A'
            >>> print errorCounter.getCount(kmerBefore='AA')
            123
            >>> ####  Count all the errors followed by kmer 'AA'
            >>> print errorCounter.getCount(kmerAfter='AA')
            12
            >>> ####     Count all the errors for truth 'A' preceded by an A
            >>> print errorCounter.getCount(truth='A',kmerBefore'A')
            123
            >>> ####     Count all the errors A->T preceded by A
            >>> print errorCounter.getCount(truth='A',emission='T', kmerBefore='A')
            12
            >>> ####  Count all the errors where the emmited base is 'A' preceded by any kmer
            >>> print errorCounter.getCount(emission='A')
            2
        """

        ########### Monogo DB Verison
        query = {}
        if truth:
            query['true'] = truth
        if emission:
            query['emission'] = emission
        if kmerBefore:
            query['leftFlank'] = {'$regex':kmerBefore+'$','$options': 'i'}
        if kmerAfter:
            query['rightFlank'] = {'$regex':'^'+kmerAfter, '$options': 'i'}
        if type:
            query['type'] = type
        if maxAlignedDist:
            query['alignedDist'] = {'$lte':maxAlignedDist}
        if mappedCorrectly is not None:
            query['mappedCorrectly'] = int(mappedCorrectly)
        if readPosRange:
            query['readPos'] = {'$lte':readPosRange[1],'$gte':readPosRange[0]}
        if readPerRange:
            query['readPer'] = {'$lte':readPerRange[1],'$gte':readPerRange[0]}
        if readLength:
            query['readLength'] = readLength
        if qualRange:
            query['qual'] = {'$lte':qualRange[1],'$gte':qualRange[0]}
        if tlenRange:
            if tlenRange is list or tlenRange is tuple:
                query['tlen'] = {'$lte':tlenRange[1],'$gte':tlenRange[0]}
            else:
                query['tlen'] = tlenRange
        
        mongoPointer = self.errordb[collection].find(query)
        if returnList:
            errorList = self.errordb[collection].find_errors(query)
            return mongoPointer.count(),errorList
        else:
            return mongoPointer.count()


class db_summary():
    def __init__(self,opt):
        ## Connection to mongoDB, runs without for now.
        self.errordb = {}
        self.errordb['errors'] = errordb(database=opt.dbName,collection=opt.observedErrorDBName)
        self.errordb['simulatedErrors'] = errordb(database=opt.dbName,collection=opt.simulatedErrorDBName)
        self.errordb['metaData'] = errordb(database=opt.dbName,collection='metaData')
        self.opt = opt
        self.qscores = False
        self.errorQualList = False

    def getAllQualScores(self):
        "Returns a list of quality scores of all the bases"
        logging.info("Extracting Qscores from samfile")
        return [asciiToInt(i) for i in list("".join(read.qual for read in pysam.Samfile( self.opt.samfile ).fetch()))]

    def errorDistribution(self):
        ed =  [doc['readPer'] for doc in self.errordb['errors'].find(query={},filt={'readPer'})]
        if ed:
            densityPlotterFromLists(dic={'errorDistribution':ed},opt=self.opt,filename='errorDistribution_dens').plot()
            densityPlotterFromLists(dic={'errorDistribution':ed},opt=self.opt,filename='errorDistribution_hist',binwidth=0.1).plot(geom='hist')

        simEd = [doc['readPer'] for doc in self.errordb['simulatedErrors'].find(query={},filt={'readPer'})]
        if simEd:
            densityPlotterFromLists(dic={'errorDistribution':simEd},opt=self.opt,filename='simulatedErrorDistribution_dens').plot()
            densityPlotterFromLists(dic={'errorDistribution':simEd},opt=self.opt,filename='simulatedErrorDistribution_hist',binwidth=0.1).plot(geom='hist')

    def errorQualDistribution(self):
        if not self.errorQualList:
            self.errorQualList =  [doc['qual'] for doc in self.errordb['errors'].find(query={'type' : {'$ne': 'Deletion'}},filt={'qual'})]
        if ed:
            densityPlotterFromLists(dic={'errorQualDistribution':self.errorQualList},opt=self.opt,filename='errorQualDistribution_dens').plot()
            densityPlotterFromLists(dic={'errorQualDistribution':self.errorQualList},opt=self.opt,filename='errorQualDistribution_hist').plot(geom='hist')

        simEd = [doc['qual'] for doc in self.errordb['simulatedErrors'].find(query={},filt={'qual'})]
        if simEd:
            densityPlotterFromLists(dic={'simulatederrorQualDistribution':simEd},opt=self.opt,filename='simulatederrorQualDistribution_dens').plot()
            densityPlotterFromLists(dic={'simulatederrorQualDistribution':simEd},opt=self.opt,filename='simulatederrorQualDistribution_hist').plot(geom='hist')        

    def qualDistribution(self):
        if not self.qscores:
            self.qscores = self.getAllQualScores()

        densityPlotterFromLists(dic={'QualDistribution':self.qscores},opt=self.opt,filename='QualDistribution_dens').plot()
        densityPlotterFromLists(dic={'QualDistribution':self.qscores},opt=self.opt,filename='QualDistribution_hist').plot(geom='hist')

    def qScoreCalibrationTest(self,variantType='SNP'):
        "plots assigned qscore versus emperical qscore to assess how well calibrated the quals are"
        logging.info("Assess qscore calibration")
        if not self.qscores:
            self.qscores = self.getAllQualScores()
        qScoreCounts = listCounter(self.qscores)

        self.errorQualList =  [doc['qual'] for doc in self.errordb['errors'].find(query={'type' : variantType},filt={'qual'})]
        errorQscoreCounts = listCounter(self.errorQualList)
        empiricalQscore = []
        samQscore =  [qual for qual,count in qScoreCounts.iteritems() if count > 100]
        
        ## Calculate empirical qscore. 
        ## prob(qual) = number of errors with that qual / number of bases with that qual        
        for qual in samQscore:
            empiricalQscore.append(probToQscore(float(errorQscoreCounts[qual]) / float(qScoreCounts[qual])))

        scatterPloter(x=samQscore,y=empiricalQscore,
                    xlab='Reported',ylab='Empirical',
                    filename='qscoreCalibration_%s' % (variantType),opt=self.opt).plot()

    @property 
    def numErrors(self):
        "Get the number of errors in the database"
        return self.errordb['errors'].find().count()



class samReader():
    """ Read the same file and return human readable alignments"""
    def __init__(self,samfile,ref):
        self.samfile = pysam.Samfile( samfile ).fetch()
        self.ref = ref
    def __iter__(self):
        if self is not None:
            return self

    def __readNext(self):
        """
        Iterates to the next read in samfile
        """
        ## If there are errors left in the read 
        self.alignedRead = self.samfile.next()
        if self.alignedRead.is_unmapped:
            self.__readNext()
        else:
            self.__parseRead()


    def next(self):
        """
        Returns the next read in the samfile aligned read
        """
        self.__readNext()
        return self

            
    @property   
    def ID(self):
        return int(self.alignedRead.qname.split('id=')[1])

    def __parseRead(self):
        """
        parses the read and returns readable alignment
        """
        self.currentRead = []
        self.currentRefRead = []
        self.currentQual = []

        self.__currentRefReadList = list(self.__refRead)
        self.__currentReadList = list(str(self.alignedRead.seq))
        self.__currentQualList = list(self.alignedRead.qual)
        self.__readPos = 0
        self.__readPosIndex = 0
        for tup in self.alignedRead.cigar:
            cigarInt = tup[0]
            numBases = tup[1]
            if cigarInt == 0:
                ## Match or mismatch
                self.__checkSNPs(N=numBases)
            elif cigarInt == 1:                
                ## Insertion to the reference
                self.__checkInsertion(N=numBases)
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

        assert len(self.currentRefRead) == len(self.currentRead)
        assert len(self.currentQual) == len(self.currentRead)
    @property
    def __refRead(self):
        """
        Get the reference sequence the read is aligned to
        """
        refRead = [str(self.ref[pos]) for pos in self.alignedRead.positions]
        return "".join(refRead)
    

    def __checkSNPs(self,N):
        """
        Checks read segment for SNP errors. called when cigarstring = M:N 
        """

        readSeg = popLong(self.__currentReadList,0,N)
        refSeg = popLong(self.__currentRefReadList,0,N)
        qualSeg = popLong(self.__currentQualList,0,N)

        self.currentRead.extend(readSeg)
        self.currentRefRead.extend(refSeg)
        self.currentQual.extend(qualSeg)

        self.__readPos += N
        self.__readPosIndex += N  

    def __checkInsertion(self,N):
        """
        Checks Insertion read segment for errors. called when cigarstring = I:N 
        """

        ## add _ to reference
        insSeg = popLong(self.__currentReadList,0,N)
        insQualSeg = popLong(self.__currentQualList,0,N)

        self.currentRead.extend(insSeg)
        self.currentQual.extend(insQualSeg)
        self.currentRefRead.extend('_'*len(insSeg))
        self.__readPos += N
    def __checkDeletion(self,N):
        """
        Checks Deletion read segment for  errors. called when cigarstring = D:N 
        """
        ## add _ to read
        i = self.alignedRead.positions[self.__readPosIndex-1] + 1
        j =  self.alignedRead.positions[self.__readPosIndex]        
        delSeg = str(self.ref[i:j])
        self.currentRead.extend('_'*len(delSeg))
        self.currentQual.extend('_'*len(delSeg))
        self.currentRefRead.extend(delSeg)

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
        clippedSeg = popLong(self.__currentReadList,0,N)
        self.currentRead.extend(clippedSeg)
        self.currentRefRead.extend('*'*len(clippedSeg))

        
    def __checkHardClipped(self,N):
        """
        Checks HardClipped read segment for errors. called when cigarstring = H:N 
        """
        logging.warning("We shouldn't have HardClipped bases in Samfile")
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

    @property
    def read(self):
        return "".join(self.currentRead)
    @property
    def refRead(self):
        return "".join(self.currentRefRead) 
    @property
    def qual(self):

        return "".join(self.currentQual) 







