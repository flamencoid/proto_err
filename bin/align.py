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
parser.add_option("-a","--algorithm",dest="algorithm",help="alignment algorithm (last(default)/bwa-mem/stapmy)",
						default='last')
# parser.add_option("-o","--outfile",dest="samfileName",help="outputsamfilename",
# 						default='')

(opt, args) = parser.parse_args()
## Make run ID mandatory
if not opt.fastqfile:
	logging.error("Please specify a run ID with -f '''fasta''' ")
	raise ValueError("-f option is mandatory")
# if not opt.samfileName:
# 	logging.error("Please specify an output filename with -o '''--outfile''' ")
# 	raise ValueError("-o option is mandatory")
## Index to the reference
opt.readFilename = opt.fastqfile

logging.info("Indexing reference")

align.refIndex(reference=opt.refFilename,algorithm=opt.algorithm )
# ## Align reads to the reference
logging.info("Aligning reads to reference")

if opt.algorithm  == "last":
	samfileName = "/".join(opt.fastqfile.split('/')[:-2])+"/sam/" + opt.fastqfile.split('/')[-1].split('.')[0] +"_%s.txt" % opt.algorithm 
else:
	samfileName = "/".join(opt.fastqfile.split('/')[:-2])+"/sam/" + opt.fastqfile.split('/')[-1].split('.')[0] +"_%s.sam" % opt.algorithm 

aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName,algorithm=opt.algorithm )

