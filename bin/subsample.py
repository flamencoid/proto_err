#!/usr/bin/env python
## python script that subsamples from a referecnce, induces errors in the reads 
## and realigns to the reference outputting a sam file
import sys
import os
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
sys.path.insert(0, os.path.abspath('../proto_err'))
from fastaIO import getRef,writeFasta,writeFastq
from simulation import complexError
from optparse import OptionParser
import align
import pysam
from query import errordb
from utils import *
import numpy as np
import math
import random 
from Bio.SeqRecord import SeqRecord
import csv
import re

parser = OptionParser()
parser.add_option("-r", "--ref", dest="refFilename",help="fasta input ref file",
					default="../data/ref.fa")
parser.add_option("-i","--id",dest="simID",help="simulation identifier",
						default='')
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
					default=5,type='int')
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
if not opt.snpFreqSd:
	opt.snpFreqSd = float(opt.snpFreq)/5
if opt.errorBiasFile:
	errorBias = AutoVivification()
	with open(opt.errorBiasFile,'rb') as errorFile:
		pattern = re.compile("#")
		reader = csv.reader(errorFile,delimiter='\t')
		for row in reader:
			if not pattern.search(row[0]):
				if row[1]:
					patternKey = re.compile(row[0]+row[1]+row[2])
				else:
					## This is incase of the regular expression matching the 
					## start or end of a line
					patternKey = re.compile(row[0]+'.'+row[2])
				errorBias[patternKey] = (len(row[0]),float(row[3]),float(row[4]))
else:
	errorBias = None

opt.dbName = 'proto_err_' + opt.simID	
opt.simulatedErrorDBName = 'simulatedErrors'

errordb(database=opt.dbName).addMetaData(opt=opt,t='simulation',errorBias=errorBias)

def subsample(ref,opt,errorBias=None,errorSimulator=complexError):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	refLength =  len(ref)
	seqList = []
	simulatedErrorDB = errordb(database=opt.dbName,collection=opt.simulatedErrorDBName)
	counter = AutoVivification()

	## Add options to errorDB
	optionsErrorDB = errordb(database=opt.dbName,collection='metaData')

	
	simulatedErrorDB.deleteAll()
	for i in range(opt.numReads):
		seqLength = abs(int(math.ceil(np.random.normal(opt.readMean,opt.readSd))))
		start = random.randrange(refLength)
		## randomly subsample from reference
		recordId = 'st=%s' % (str(start))
		seq = ref[start:start+seqLength]
		record=SeqRecord(seq,recordId,'','')
		## Take the read from the reverse stand opt.strandBias% of the time
		if random.random() > opt.strandBias:
			record = record.reverse_complement()
		## Randomly generate errors
		simulatedErrors = errorSimulator(record,opt,id = recordId,errorBias=errorBias)
		errs = simulatedErrors.error()
		logging.info("### generated %i errors in a read of length %i" % (len(errs),seqLength))
		simulatedErrorDB.addErrors(errs)
		record = simulatedErrors.record
		seqList.append(record)
	return seqList

opt.readFilename = opt.refFilename[:-3] + '.subsampled.fq' 
ref = getRef(opt.refFilename)
logging.info("Subsampling reads from reference")
seqList = subsample(ref,opt,errorBias=errorBias)



logging.info("Writing Fasta file of subsampled reads")
writeFastq(filename = opt.readFilename,seqList = seqList)
# ## Index to the reference
logging.info("Indexing reference")
align.refIndex(file=opt.refFilename)
# ## Align reads to the reference
logging.info("Aligning reads to reference")
samfileName = opt.readFilename + '.sam'
aligned = align.align(reference=opt.refFilename, read_file=opt.readFilename,stdout=samfileName)







