#!/usr/bin/env python
## python script that aligns to the reference outputting a sam file
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from optparse import OptionParser
import align

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
parser.add_option("-i","--id",dest="simID",help="simulation identifier",
						default='')
(opt, args) = parser.parse_args()
## Make run ID mandatory
if not opt.simID:
	logging.error("Please specify a run ID with -i '''id''' ")
	raise ValueError("-i option is mandatory")

## Index to the reference
opt.readFilename = opt.refFilename[:-3] +'.subsampled.'+ opt.simID+'.fq' 

logging.info("Indexing reference")
align.refIndex(file=opt.refFilename)
# ## Align reads to the reference
logging.info("Aligning reads to reference")
samfileName = opt.refFilename[:-3] +'.subsampled.'+ opt.simID + '.sam'
aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName)
