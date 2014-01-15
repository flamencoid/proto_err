#!/usr/bin/env python
## A script to take a sam file of aligned reads and do some stats on errors

import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from optparse import OptionParser
from errorCount import errorReader,counter,summary
from fastaIO import getRef
from query import errordb
import time
start = time.clock()

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
# parser.add_option("-s", "--samfile", dest="samfile",help="Samfile of aligned reads",
# 						default="../data/ref.subsampled.fq.sam")
parser.add_option("--outDir", dest="outDir",help="Path to output directory (Optional)",
						default="../results")
parser.add_option("-f","--forceUseCurrentDB",action="store_false", dest="force",help="Don't Remake database (Optional: default wipes and repopulates database)",
						default=True)
parser.add_option("-i","--id",dest="simID",help="simulation identifier",
						default='')
(opt, args) = parser.parse_args()

## Hardcode some options
opt.maxKmerLength = 3 
opt.outDir = opt.outDir + '/' + opt.simID + '/'
opt.imgDir = opt.outDir+'img/'
opt.jsonDir = opt.outDir +'json/'
if not os.path.exists(opt.outDir):
    os.makedirs(opt.outDir)
if not os.path.exists(opt.imgDir):
    os.makedirs(opt.imgDir)
if not os.path.exists(opt.jsonDir):
    os.makedirs(opt.jsonDir)
## Make run ID mandatory
if not opt.simID:
	logging.error("Please specify a run ID with -i '''id''' ")
	raise ValueError("-i option is mandatory")
opt.dbName = 'proto_err_' + opt.simID.replace('.','_')
opt.samfile = "../data/ref.subsampled."+opt.simID+".sam"
opt.simulatedErrorDBName = 'simulatedErrors'
opt.observedErrorDBName = 'errors'
ref = getRef(opt.refFilename)

if opt.force:
	errorCounter = counter(ref,opt,samfile=opt.samfile,makeDB=True)
else:
	errorCounter = counter(ref,opt,samfile=opt.samfile,makeDB=False)

summ = summary(opt)
summ.errorDistribution()
summ.qualDistribution()
summ.qScoreCalibrationTest()

# errorCounter.summary()
# errorCounter.plotHist()

end = time.clock()
logging.info("errorStats.py took %i seconds to run" % (end-start))


