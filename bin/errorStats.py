#!/usr/bin/env python
## A script to take a sam file of aligned reads and do some stats on errors

import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../src'))
from optparse import OptionParser
from errorCount import errorReader,counter,db_summary,samReader
from Bio import SeqIO
from query import errordb
import time
#from report import Reporter
from fastaIO import getRef
start = time.clock()

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file")
# parser.add_option("-s", "--samfile", dest="samfile",help="Samfile of aligned reads",
# 						default="../data/ref.subsampled.fq.sam")
parser.add_option("--outDir", dest="outDir",help="Path to output directory (Optional)",
						default="../results")
parser.add_option("-f","--forceUseCurrentDB",action="store_false", dest="force",help="Don't Remake database (Optional: default wipes and repopulates database)",
						default=True)
parser.add_option("-s","--samfile",dest="samfile",help="samfile of alignedReads",
						default='')
(opt, args) = parser.parse_args()

## Hardcode some options
opt.maxKmerLength = 3 
opt.runID = opt.samfile.split('/')[-1].split('.')[0]
opt.outDir = opt.outDir + '/' + opt.runID + '/'
opt.imgDir = opt.outDir+'img/'
opt.jsonDir = opt.outDir +'json/'
if not os.path.exists(opt.outDir):
    os.makedirs(opt.outDir)
if not os.path.exists(opt.imgDir):
    os.makedirs(opt.imgDir)
if not os.path.exists(opt.jsonDir):
    os.makedirs(opt.jsonDir)
## Make run ID mandatory
if not opt.samfile:
	logging.error("Please specify a samfile of aligned reads '''samfile''' ")
	raise ValueError("-s option is mandatory")
opt.dbName = 'proto_err_' + opt.runID
opt.observedErrorDBName = 'errors'
opt.observedReadDBName = 'alignedReads'
ref = getRef(opt.refFilename)
observedReadsDB = errordb(database=opt.dbName,collection=opt.observedReadDBName )


# Generate read documents for uploading to db
# if opt.force:
# 	observedReadsDB.deleteAll()
# 	for read in samReader(samfile=opt.samfile,ref=ref):
# 		post = {'id':read.ID,'read':read.read,'ref':read.refRead,
# 		'qual':read.qual,'cigar':read.alignedRead.cigarstring}
# 		observedReadsDB.insert(post)


if opt.force:
	errorCounter = counter(ref,opt,samfile=opt.samfile,makeDB=True)
else:
	errorCounter = counter(ref,opt,samfile=opt.samfile,makeDB=False)
errorCounter.INSTransitionStats()
errorCounter.DELTransitionStats()

errorCounter.SNPTransitionStats()
errorCounter.summary()
errorCounter.plotHist()
# # ## Do some meta and summary statistics
summ = db_summary(opt)
summ.errorDistribution()
summ.qualDistribution()
summ.qScoreCalibrationTest('SNP')
summ.qScoreCalibrationTest('Insertion')
summ.qScoreCalibrationTest('Deletion')

#report = Reporter(opt=opt,counter=errorCounter,outfileDir= opt.outDir ,latexTemplate='../data/template.tex')
#report.generatePdfReport()

end = time.clock()
logging.info("errorStats.py took %f seconds to run" % (end-start))


