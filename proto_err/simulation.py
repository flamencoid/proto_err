#!/usr/bin/env python
from utils import *
import random 
import numpy as np
import math
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import copy
# A python module for simulating errors in reads.
class simulateError():
	"""Class of objects to simulate errors in a SeqRecord"""
	def __init__(self,record,opt,id):
		self.alphabet = ['A','T','C','G']
		self.seq = record.seq
		self.id = id
		self.opt = opt
		self.errorProb = [0]*len(self.seq)

	def snp(self,pos,rl):
		"""function to induce a SNP"""
		seq = list(self.seq)
		seq[pos] = rl
		self.seq = Seq("".join(seq))
		self.id += 's%s,%s' % (str(pos),str(rl))

	def ins(self,pos,rl):
		"""function to induce a insertion"""
		seq = list(self.seq)
		seq = seq[:pos+1] + list(rl) + seq[pos+1:]
		## Also need to insert the equivalent qscores
		self.errorProb = self.errorProb[:pos+1] + [self.errorProb[pos] for _ in list(rl)]  +self.errorProb[pos+1:]
		self.seq = Seq("".join(seq))
		self.id += 'i%s,%s' % (str(pos),str(rl))

	def deletion(self,pos,dlen):
		"""function to induce a deletion"""
		seq = list(self.seq)
		seq = seq[:pos] + seq[pos+dlen:]
		## Also need to delete the equivalent qscores
		self.errorProb = self.errorProb[:pos] + self.errorProb[pos+dlen:]
		self.seq = Seq("".join(seq))
		self.id += 'd%s,%s' % (str(pos),str(dlen))

	def indel(self,pos):
		"""function to induce an INDEL"""
		indelSize = int(round(random.gauss(self.opt.indelMean,self.opt.indelSd)))
		if random.random() < 0.5:
			self.ins(pos,rl="".join([random.choice(self.alphabet) for _ in range(indelSize)]))
		else:
			self.deletion(pos,dlen=indelSize)

	@property 
	def record(self):
		"""Return a Bio.SeqRecord"""
		return SeqRecord(self.seq,self.id,letter_annotations={'phred_quality':self.qscore})

	@property
	def qscore(self):
		return [int(-10 * math.log10(p+0.000001)) for p in self.errorProb]

class singleSNP(simulateError):
	def __init__(self,record,opt,id):
		simulateError.__init__(self,record,opt,id)
		self.errorProb = [random.gauss(opt.snpFreq, 0.01) for _ in range(len(self.seq))] 

	def error(self):
		"""Function to induce and error"""
		## iterate through the probability list
		for pos,prob in enumerate(self.errorProb):
			if random.random() < prob:
				letter = self.seq[pos]
				alphabet = copy(self.alphabet)
				reducedAlphabet = alphabet
				reducedAlphabet.remove(letter)
				replaceLetter = random.choice(reducedAlphabet)
				self.snp(pos,replaceLetter)

class complexError(simulateError):
	def __init__(self,record,opt,id,baseErrorProb=None):
		simulateError.__init__(self,record,opt,id)
		if baseErrorProb:
			self.errorProb = baseErrorProb
			assert len(self.errorProb) == len(self.seq)
		else:
			self.errorProb = [random.gauss(opt.snpFreq, 0.01) for _ in self.seq]

	def error(self):
		"""Function to induce and error"""
		## iterate through the probability list
		for pos,prob in enumerate(self.errorProb):
			if random.random() < prob and pos < len(self.errorProb):
				letter = self.seq[pos]
				alphabet = copy(self.alphabet)
				reducedAlphabet = alphabet
				reducedAlphabet.remove(letter)
				replaceLetter = random.choice(reducedAlphabet)
				if random.random() < self.opt.SnpIndelRatio:
					self.snp(pos,replaceLetter)
				else:
					self.indel(pos)
			

def subsample(ref,opt,errorSimulator=complexError):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	refLength =  len(ref)
	seqList = []
	for i in range(opt.numReads):
		seqLength = abs(int(math.ceil(np.random.normal(opt.readMean,opt.readSd))))
		start = random.randrange(refLength)
		## randomly subsample from reference
		recordId = 'st=%s&l=%s' % (str(start),str(seqLength))
		seq = ref[start:start+seqLength]
		record=SeqRecord(seq,recordId,'','')
		## Take the read from the reverse stand x% of the time
		if random.random() > opt.strandBias:
			record = record.reverse_complement()
		## Randomly generate errors
		simulatedErrors = errorSimulator(record,opt,id = recordId)
		simulatedErrors.error()
		record = simulatedErrors.record
		seqList.append(record)
	return seqList