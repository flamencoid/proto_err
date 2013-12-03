#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time
from utils import *
import difflib
import itertools


## Compare aligned reads to reference and calculate stats
class comparison(): 

    def __init__(self,ref):
        # the results dictionary
        self.res = {}
        # set toplevel name 
        self.VERSION  = "v1.0"
        self.metaLabels = {'name':'proto_err','timestamp': 
                            time.strftime("%a %b %d %X %Y"),\
                           'runInfo':[],'version':self.VERSION,\
                           'metrics':[]}
        self.counts  = [['Total','Total number of reads',0],
                         ['Mapped','Total number of aligned reads',0],
                         ['UnMapped','Total number of non-aligned reads',0],
                         ['NM','Total number of SNP errors',0],
                         ['perfectAlignments','Total number of perfectly aligned reads',0],
                         ['mismatchedAlignments','Total number of aligned reads with mismatches',0]]
        self.ref = ref
        self.logger     = logging.getLogger()


    def setup(self,opt):
        """
        Function to set up dictionary structure for outputed stats
        """
        self.opt = opt
        self.res['Counts'] = {}
        self.res['kmerCounts'] = AutoVivification()
        self.res['RefCounts'] =  AutoVivification()
        self.errorList = []
        for metric in self.counts:
            self.res['Counts'][metric[0]] = metric[2]

        ## Create the dictonary of error modes
        ## res[A][T] = counts of A -> T SNPs

        self.res['errorMode'] = dict(zip(getAlphabet(),[dict(zip(getAlphabet(),[0]*4)) for i in getAlphabet()]))

    def kmerFreq(self,seq,kmer):
        """
        Function to calculate the number of times a kmer appears in a sequence
        """
        count = seq.count(kmer)
        return count,float(count)/len(seq)


    def countKmers(self,kmerLen):
        """
        Function to count all kmers in long sequence (reference) of length kmerLen or below
        """
        alpabet = getAlphabet()
        for klen in [i+1 for i in range(kmerLen)]:
            # generate all kmers of length klen
            kmerList = kmerCombo(alpabet,klen)
            self.res['RefCounts'] = dict(zip(kmerList,[0]*len(kmerList)))
            for kmer in kmerList:
                self.res['RefCounts'][str(kmer)] = self.ref.count(kmer)

    def precedingKmers(self):
        """
        Function which takes a list of errors and counts the kmers before  
        and after and error
        """
        for error in self.errorList:
            for j in [k +1 for k in range(self.opt.maxOrder)]:
                try:
                    self.res['kmerCounts']['before'][error.true][error.emission][error.before(j)] += 1
                except:
                    self.res['kmerCounts']['before'][error.true][error.emission][error.before(j)] = 1
                try:
                    self.res['kmerCounts']['after'][error.true][error.emission][error.after(j)] += 1
                except:
                    self.res['kmerCounts']['after'][error.true][error.emission][error.after(j)]= 1


    def compareReads(self,samfile,reffile):
        """
        Function which take a samfile iterates through the aligned reads and 
        generates stats
        """
        samfile = pysam.Samfile( samfile )
        for read in samfile.fetch():
            self.res['Counts']['Total'] += 1
            if read.is_unmapped:
                self.res['Counts']['UnMapped'] += 1
            else:
                self.res['Counts']['Mapped'] += 1
                self.checkRead(read)

    def readDiff(self,read,ref):
        return difflib.ndiff(read,ref)

    def getRefRead(self,positions):
        """
        Function to return the reference sequence a read is aligned to
        """
        return self.ref[positions[0]:positions[-1]+1]

    def checkRead(self,read):
        """
        A function which take a read and returns some error stats
        """
        NM = read.opt('NM')
        self.res['Counts']['NM'] += NM
        if NM == 0:
            self.res['Counts']['perfectAlignments'] += 1
            for base in getAlphabet():
                self.res['errorMode'][base][base] += read.seq.count(base)
        else:
            self.res['Counts']['mismatchedAlignments'] += 1
            ## Check what type of mismatch it is.
            ## if the cigarsting only has M then only SNP errors
            if (read.rlen == read.cigar[0][1]):
                refRead = list(str(self.getRefRead(read.positions)))

                for i,tupl in enumerate(itertools.izip_longest(refRead,list(str(read.seq)))):
                    true,emission = tupl
                    self.res['errorMode'][true][emission] += 1
                    ## if it's an error, check the preceding  opt.maxOrder bases
                    if not true == emission:
                        ## Check preceding bases
                        self.errorList.append(error(true,emission,read.seq[i-self.opt.maxOrder:i],read.seq[i+1:self.opt.maxOrder+i+1]))


class error():
    """ Information about the errors in a read """
    def __init__(self,true,emission,leftFlank='',rightFlank=''):
        self.true = true
        self.emission = emission
        self.leftFlank = leftFlank 
        self.rightFlank = rightFlank

    def before(self,j):
        """Return the preceding j bases,return N when bases missing"""
        b = self.leftFlank[-j:]
        while len(b) < j:
            b = 'N' +b
        return b
    def after(self,j):
        """Return the following j bases,return N when bases missing"""
        a = self.rightFlank[0:j]
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
    def flankLength(self):
        """Lengths of the flanks either side"""
        return  (len(self.leftFlank),len(self.rightFlank))
    @property 
    def isSnp(self):
        """Is the error a SNP"""
        return len(self.true) == len(self.emission)
    @property 
    def isIndel(self):
        """Is the error an INDEL"""
        return len(self.true) != len(self.emission)

class alignedRead():
    """
    Class for the aligned read from a samfile
    """

    pass









