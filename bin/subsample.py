#!/usr/bin/env python

## python script that will take a fasta file and return a fasta file 
## of subsampled reads
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from fastaIO import getRef,writeFasta
from error import singleSNP,subsample,identity
from optparse import OptionParser
from metrics import comparison 
import align
import pysam

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",default="../data/ref.fa")
(opt, args) = parser.parse_args()

## Hardcod some numbers
opt.snpFreq = 0.1
opt.maxOrder = 3

opt.readFilename = opt.refFilename[:-3] + '.subsampled.fa' 
ref = getRef(opt.refFilename)
logging.info("Subsampling reads from reference")
seqList = subsample(ref,opt,readError=singleSNP,numReads=1000,readRange = [1000,20000])
logging.info("Writing Fasta file of subsampled reads")
writeFasta(filename = opt.readFilename,seqList = seqList)
# ## Index to the reference
align.refIndex(file=opt.refFilename)
# ## Align reads to the reference
samfileName = opt.readFilename + '.sam'
aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName)

logging.info("Doing read comparision")
compare = comparison(ref)
compare.setup(opt)
logging.info("countKmers")
compare.countKmers(opt.maxOrder)
logging.info("countKmers")
logging.info("compareReads")

compare.compareReads(samfile=samfileName,reffile=opt.readFilename)

# print compare.res
for ee in compare.errorList:
	print ee.true,ee.emission,ee.seq






     # print read.cigarstring
     # print read.NM

# pysam.sort( "-S",samfileName, "output" )
# print pysam.idxstats("-S",samfileName)

# '__class__', '__delattr__', '__doc__', '__format__', '__getattribute__', 
# '__hash__', '__init__', '__new__', '__reduce__', '__reduce_ex__', '__repr__',
 # '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'aend', 'alen', 
 # 'aligned_pairs', 'bin', 'cigar', 'cigarstring', 'compare', 'fancy_str', 'flag', 
 # 'is_duplicate', 'is_paired', 'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 
 # 'is_reverse', 'is_secondary', 'is_unmapped', 'isize', 'mapq', 'mate_is_reverse', 
 # 'mate_is_unmapped', 'mpos', 'mrnm', 'opt', 'overlap', 'pnext', 'pos', 'positions', 
 # 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'rlen', 'rname', 
 # 'rnext', 'seq', 'tags', 'tid', 'tlen']









