#!/usr/bin/env python
## python script that subsamples from a referecnce, induces errors in the reads 
## and realigns to the reference outputting a sam file
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from fastaIO import getRef,writeFasta,writeFastq
from simulation import subsample
from optparse import OptionParser
import align
import pysam

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
parser.add_option("--numReads", dest="numReads",help="""Number of reads to 
					sample from referecnce (Optional defaults to 1000)""",
					default=1000,type='int')
parser.add_option("--meanReadLength", dest="readMean",help="""Read length is 
					sampled from a normal of this mean (Optional defaults to 
						1000)""",default=1000,type='int')
parser.add_option("--errorFreqMean", dest="snpFreq",help="""Probablity of an error 
					at a base occurs is sampled from a normal with mean 
					of this Probablity (Optional defaults to 
						.1)""",default=0.1,type='float')
parser.add_option("--errorFreqSD", dest="snpFreqSd",help="""Probablity of an error 
					at a base occurs is sampled from a normal with SD 
					of this value (Optional defaults to 
						errorFreqMean/10)""",default=None,type='float')
parser.add_option("--strandBias", dest="strandBias",help="""Ratio of reads 
					sampled from forward or reverse strand. 1 samples only from 
					forward, 0 only from reverse.(Optional defaults to .5)""",
					default=0.5,type='float')
parser.add_option("--ReadLengthSD", dest="readSd",help="""Read length is 
					sampled from a normal with this SD (Optional defaults to 
						meanReadLength/3)""",default=None,type='int')
parser.add_option("--SnpIndelRatio", dest="SnpIndelRatio",help="""Ratio of SNP 
					errors to INDEL errors (Optional defaults to 0.5)
					value of 1 gives ONLY SNP errors, 0 only INDELs""",
					default=0.5,type='float')
parser.add_option("--meanIndelSize", dest="indelMean",help="""INDEL size is 
					sampled from a normal of this mean (Optional defaults to 
						5)""",
					default=5,type='int')
parser.add_option("--IndelSizeSD", dest="indelSd",help="""INDEL size is 
					sampled from a normal with this SD (Optional defaults to 
						meanIndelSize/2)""",default=None,type='int')
(opt, args) = parser.parse_args()
## Need a mean and a SD
if not opt.readSd :
	opt.readSd = int(float(opt.readMean)/3)
if not opt.indelSd:
	opt.indelSd = float(opt.indelMean)/2
if not opt.snpFreqSd:
	opt.snpFreqSd = float(opt.snpFreq)/10

opt.simID = 'tmp'

opt.readFilename = opt.refFilename[:-3] + '.subsampled.fq' 
ref = getRef(opt.refFilename)
logging.info("Subsampling reads from reference")
seqList = subsample(ref,opt)
logging.info("Writing Fasta file of subsampled reads")
writeFastq(filename = opt.readFilename,seqList = seqList)
# ## Index to the reference
logging.info("Indexing reference")
align.refIndex(file=opt.refFilename)
# ## Align reads to the reference
logging.info("Aligning reads to reference")
samfileName = opt.readFilename + '.sam'
aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName)






