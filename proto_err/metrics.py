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
        self.res['Counts'] = {}
        self.res['kmerCounts'] = {}
        self.res['RefCounts'] = {}
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
            self.res['RefCounts'][klen] = dict(zip(kmerList,[0]*len(kmerList)))
            for kmer in kmerList:
                self.res['RefCounts'][klen][str(kmer)] = self.ref.count(kmer)

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
                # print "\n".join(difflib.ndiff([str(refRead)], [str(read.seq)]))
                for t,e in itertools.izip_longest(refRead,list(str(read.seq))):
                    self.res['errorMode'][t][e] += 1




                # for s in self.readDiff(read.seq,refRead):
                #     if s[0] == ' ':
                #         self.res['errorMode'][s[2]][s[2]] += 1
                #     else:
                #         print s,s.next()
                        

                # print read
                # print read.positions[0],read.positions[-1]
                # print len(refRead)
                # print read.rlen
                # print read.tlen

                # for pos,emmittedBase in enumerate(read.seq):
                #     trueBase = self.ref[read.positions[pos]]
                #     # print str(trueBase),str(emmittedBase)
                #     self.res['errorMode'][str(trueBase)][str(emmittedBase)] += 1

        

class alignedRead():
    """
    Class for the aligned read from a samfile
    """

    pass









