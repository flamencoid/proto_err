#!/usr/bin/env python
## python script that aligns to the reference outputting a sam file
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../src'))
from optparse import OptionParser
import align

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
parser.add_option("-f","--reads",dest="fastqfile",help="fasta of reads",
						default='')
(opt, args) = parser.parse_args()
## Make run ID mandatory
if not opt.fastqfile:
	logging.error("Please specify a run ID with -f '''fasta''' ")
	raise ValueError("-f option is mandatory")

## Index to the reference
opt.readFilename = opt.fastqfile

logging.info("Indexing reference")
# align.refIndex(file=opt.refFilename)
# ## Align reads to the reference
logging.info("Aligning reads to reference")
samfileName = opt.refFilename[:-3] +'.sam'
aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName,algorithm='bwa-mem')

