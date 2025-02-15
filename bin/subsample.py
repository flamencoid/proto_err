#!/usr/bin/env python
## python script that subsamples from a referecnce, induces errors in the reads 
## and realigns to the reference outputting a sam file
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../src'))
from fastaIO import getRef,writeFasta,writeFastq
from simulation import complexError
from optparse import OptionParser
import pysam
from query import errordb
from utils import *
import numpy as np
import math
import random 
from Bio.SeqRecord import SeqRecord
import csv
import re
import time
## set seed

start = time.clock()
parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/refs/tb.ref.fasta")
parser.add_option("-i","--id",dest="simID",help="simulation identifier",
						default='')
parser.add_option("--numReads", dest="numReads",help="""Number of reads to 
					sample from referecnce (Optional defaults to 1000)""",
					default=100,type='int')
parser.add_option("--meanReadLength", dest="readMean",help="""Read length is 
					sampled from a normal of this mean (Optional defaults to 
						5000)""",default=100,type='int')
parser.add_option("--errorFreqMean", dest="errFreq",help="""Probablity of an error 
					at a base occurs is sampled from a normal with mean 
					of this Probablity (Optional defaults to 
						.1)""",default=0.01,type='float')
parser.add_option("--errorFreqSD", dest="errFreqSd",help="""Probablity of an error 
					at a base occurs is sampled from a normal with SD 
					of this value (Optional defaults to 
						errorFreqMean/5)""",default=None,type='float')
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
					default=3,type='int')
parser.add_option("--IndelSizeSD", dest="indelSd",help="""INDEL size is 
					sampled from a normal with this SD (Optional defaults to 
						meanIndelSize/2)""",default=None,type='int')
parser.add_option("--errorBiasFile", dest="errorBiasFile",help="""File 
					containing information about error biases 
					(Optional) default all errors occur with probability 
					--errorFreqMean""",default=None)
(opt, args) = parser.parse_args()
## Need a mean and a SD
if not opt.readSd :
	opt.readSd = int(float(opt.readMean)/3)
if not opt.indelSd:
	opt.indelSd = float(opt.indelMean)/2
if not opt.errFreqSd:
	opt.errFreqSd = float(opt.errFreq)/5
if opt.errorBiasFile:
	errorBias = AutoVivification()
	with open(opt.errorBiasFile,'rb') as errorFile:
		pattern = re.compile("#") # pattern to skip comments
		reader = csv.reader(errorFile,delimiter=',')
		for row in reader:
			## If the row isn't a comment
			if not pattern.search(row[0]):
				if row[2]:
					# If the bias includes a specific base
					patternKey = re.compile(row[1]+row[2]+row[3])
				else:
					## This is incase of the regular expression matching the 
					## start or end of a line
					patternKey = re.compile(row[1]+'.'+row[3])
				errorBias[patternKey] = (row[0],len(row[1]),float(row[4]),float(row[5]))
else:
	errorBias = None

## Make run ID mandatory
if not opt.simID:
	logging.error("Please specify a run ID with -i '''id''' ")
	raise ValueError("-i option is mandatory")

opt.dbName = 'proto_err_' + opt.simID.replace('.','_')	
opt.simulatedErrorDBName = 'simulatedErrors'
opt.simulatedReadDBName = 'simulatedReads'
# errordb(database=opt.dbName).addMetaData(opt=opt,t='simulation',errorBias=errorBias)
errordb(database=opt.dbName).addMetaData(opt=opt,t='simulation')

def subsample(ref,opt,errorBias=None,errorSimulator=complexError):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	refLength =  len(ref)
	seqList = []
	simulatedErrorDB = errordb(database=opt.dbName,collection=opt.simulatedErrorDBName)
	simulatedReadsDB = errordb(database=opt.dbName,collection=opt.simulatedReadDBName)
	counter = AutoVivification()

	## Add options to errorDB
	optionsErrorDB = errordb(database=opt.dbName,collection='metaData')

	
	simulatedErrorDB.deleteAll()
	simulatedReadsDB.deleteAll()
	for i in range(opt.numReads):
		logging.info("### Read %i of  %i" % (i+1,opt.numReads) )
		seqLength = abs(int(math.ceil(np.random.normal(opt.readMean,opt.readSd))))
		start = random.randrange(refLength)
		## randomly subsample from reference
		recordId = 'st=%s&id=%i' % (str(start),i)
		seq = ref[start:start+seqLength]
		record=SeqRecord(seq,recordId,'','')
		## Take the read from the reverse stand opt.strandBias% of the time
		opt.is_reverse = False
		if random.random() > opt.strandBias:
			record = record.reverse_complement()
			opt.is_reverse = True
		## Randomly generate errors
		opt.refPos = start
		simulatedErrors = errorSimulator(record,opt,id = recordId,errorBias=errorBias)
		errs = simulatedErrors.error()
		logging.info("### generated %i errors in a read of length %i" % (len(errs),seqLength))
		simulatedErrorDB.addErrors(errs)
		simulatedReadsDB.insert(post = {'id':i,'read':"".join(simulatedErrors.read),
										'ref':"".join(simulatedErrors.ref),
										'qscore':"".join([str(i) for i in simulatedErrors.qscore('int') ])})
		record = simulatedErrors.record
		seqList.append(record)
	return seqList

opt.readFilename = "../data/ref" +'.subsampled.'+ opt.simID+'.fq' 
ref = getRef(opt.refFilename)
logging.info("Subsampling reads from reference")
seqList = subsample(ref,opt,errorBias=errorBias)

logging.info("Writing Fasta file of subsampled reads")
writeFastq(filename = opt.readFilename,seqList = seqList)

end = time.clock()
logging.info("subsampled.py took %i seconds to run" % (end-start))






