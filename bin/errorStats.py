#!/usr/bin/env python
## A script to take a sam file of aligned reads and do some stats on errors

import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from optparse import OptionParser
from errorCount import errorReader,counter
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
	print error
print reader.readCounter

logging.info("Doing some kmer counting")
errorCounter = counter(ref,opt,samfile=opt.samfile)
errorCounter.countRefKmer()
errorCounter.countErrorKmer(1)

# for post in errorCounter.errordb.find():
# 	print post
print errorCounter.getCount(), len(errorCounter.errorList)
print errorCounter.getCount(kmerBefore='A')
print errorCounter.getCount(maxAlignedDist=10)
print errorCounter.getCount(readPosRange=[0,10])
count,errorList = errorCounter.getCount(qualRange=[10,10],returnList=True)
for error in errorList:
	print error
print count,len(errorList)
# errorCounter.res['qualCounter']



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


