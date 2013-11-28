#!/usr/bin/env python
import logging
import pysam
from pysam import AlignedRead
from fastaIO import getRef
import time


## Compare aligned reads to reference and calculate stats

class comparison(): 

    def __init__(self):
        # the results dictionary
        self.res = {}
        # set toplevel name 
        self.VERSION  = "v1.0"
        self.metaLabels = {'name':'proto_err','timestamp': 
        					time.strftime("%a %b %d %X %Y"),\
                           'runInfo':[],'version':self.VERSION,\
                           'metrics':[]}
        self.metrics  = [['Total','Total number of reads',0],
                         ['Mapped','Total number of aligned reads',0],
                         ['UnMapped','Total number of non-aligned reads',0]]
       
        self.logger     = logging.getLogger()
    def setup(self):
    	"""
    	Function to set up dictionary structure for outputed stats
    	Dic structure will be 
    	"""
    	for metric in self.metrics:
    		self.res[metric[0]] = metric[2]


    def kmerFreq(self,seq,kmer):
    	"""
    	Function to calculate the number of times a kmer appears in a sequence
    	"""
    	count = seq.count(kmer)
    	return count,float(count)/len(seq)
    def compareReads(self,samfile,reffile):
    	"""
    	Function which take a samfile iterates through the aligned reads and 
    	generates stats about error bias
    	"""
    	samfile = pysam.Samfile( samfile )
    	ref = getRef(reffile)
    	for read in samfile.fetch():
    		self.res['Total'] += 1
    		if read.is_unmapped:
    			self.res['UnMapped'] += 1
    		else:
    			self.res['Mapped'] += 1
    			self.checkRead(read)

    def checkRead(self,read):
    	"""
    	A funciton which take a read and returns some error stats
    	"""
    	pass

class alignedRead():
	"""
	Class for the aligned read from a samfile
	"""

	pass









