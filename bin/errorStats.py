#!/usr/bin/env python
## A script to take a sam file of aligned reads and do some stats on errors

import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from optparse import OptionParser
from errorCount import errorReader,counter,aggregator
from fastaIO import getRef


parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
parser.add_option("-s", "--samfile", dest="samfile",help="Samfile of aligned reads",
						default="../data/ref.subsampled.fq.sam")
(opt, args) = parser.parse_args()

## Hardcode some options
opt.maxKmerLength = 3 

ref = getRef(opt.refFilename)
# logging.info("Finding and aggregating errors
reader = errorReader(opt.samfile,ref)
for error in reader:
	if error.isIndel:
		print error
print reader.readCounter

logging.info("Doing some kmer counting")
errorCounter = counter(ref,samfile=opt.samfile)
errorCounter.setup(opt)
errorCounter.countRefKmer()
errorCounter.countErrorKmer()
print errorCounter.res['qualCounter']



# errorAggregation = aggregator(errorCounter)
# print errorAggregation.precedingKmersCount('TTC')

# logging.info("Doing read comparision")
# compare = comparison(ref)
# compare.setup(opt)
# logging.info("countKmers")
# compare.countKmers(opt.maxOrder)
# logging.info("countKmers")
# logging.info("compareReads")

# compare.compareReads(samfile=opt.samfileName,reffile=opt.readFilename)
# compare.precedingKmers()
# print compare.res
# for ee in compare.errorList:
# 	print ee.true,ee.emission,ee.seq


