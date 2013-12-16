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
parser.add_option("--outDir", dest="outDir",help="Path to output directory (Optional)",
						default="../results")
(opt, args) = parser.parse_args()

## Hardcode some options
opt.maxKmerLength = 3 
opt.outDir = opt.outDir + '/'
opt.imgDir = opt.outDir+'img/'
opt.jsonDir = opt.outDir +'json/'
if not os.path.exists(opt.outDir):
    os.makedirs(opt.outDir)
if not os.path.exists(opt.imgDir):
    os.makedirs(opt.imgDir)
if not os.path.exists(opt.jsonDir):
    os.makedirs(opt.jsonDir)

ref = getRef(opt.refFilename)

logging.info("Doing some kmer counting")
errorCounter = counter(ref,opt,samfile=opt.samfile,makeDB=True)
for name,count in errorCounter.readCounter.iteritems():
	logging.info('### Count of %s == %i' % (name,count))
errorCounter.summary()

# errorCounter.countRefKmer()
# errorCounter.countErrorKmer(1)
errorCounter.plotHist()
print errorCounter.getSimulatedCount(truth='A',emission='T')
print errorCounter.getCount(truth='A',emission='T')


# for post in errorCounter.errordb.find():
# 	print post
# print errorCounter.getCount(), len(errorCounter.errorList)
# print errorCounter.getCount(kmerBefore='A')


# print errorCounter.getCount(maxAlignedDist=10)
# print errorCounter.getCount(readPosRange=[0,10])
# count,errorList = errorCounter.getCount(qualRange=[10,10],returnList=True)
# for error in errorList:
# 	print error
# print count,len(errorList)
# errorCounter.res['qualCounter']




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


