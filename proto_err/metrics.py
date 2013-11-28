#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time
from utils import *


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
        for metric in self.counts:
            self.res['Counts'][metric[0]] = metric[2]

        ## Create the dictonary of error modes
        ## res[A][T] = counts of A -> T SNPs
        self.res['errorMode'] = dict(zip(getAlphabet(),[dict(zip(getAlphabet(),[0]*4))]*4))

    def kmerFreq(self,seq,kmer):
        """
        Function to calculate the number of times a kmer appears in a sequence
        """
        count = seq.count(kmer)
        return count,float(count)/len(seq)
    def compareReads(self,samfile,reffile):
        """
        Function which take a samfile iterates through the aligned reads and 
        generates stats
        """
        samfile = pysam.Samfile( samfile )
        ref = getRef(reffile)
        for read in samfile.fetch():
            self.res['Counts']['Total'] += 1
            if read.is_unmapped:
                self.res['Counts']['UnMapped'] += 1
            else:
                self.res['Counts']['Mapped'] += 1
                self.checkRead(read)

    def checkRead(self,read):
        """
        A function which take a read and returns some error stats
        """
        NM = read.opt('NM')
        self.res['Counts']['NM'] += NM
        if NM == 0:
            self.res['Counts']['perfectAlignments'] += 1
        else:
            self.res['Counts']['mismatchedAlignments'] += 1
            ## Check what type of mismatch it is.
            ## if the positions on reference are continous and length of 
            ## read is the same then there's only SNP errors??
            if (read.rlen == len(read.positions)) and (len(read.positions)==
                 read.positions[-1] - read.positions[0] + 1):
                # print read.tlen
                pass

        

class alignedRead():
    """
    Class for the aligned read from a samfile
    """

    pass









