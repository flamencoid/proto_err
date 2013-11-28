#!/usr/bin/env python
from utils import *
import random 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# A python module for simulating errors in reads.
def identity(seqObj):
	return Seq("".join(list(seqObj)))

def SNP(seq,pos,rl):
	seq = list(seq)
	seq[pos] = rl
	return Seq("".join(seq))

def singleSNP(seqObj,freq):
	"""
	Function to randomly add SNP errors
	"""
	seq = seqObj
	alphabet = getAlphabet()
	if random.random() < freq:
		pos = random.randint(0,len(seq))
		letter =seq[pos]
		reducedAlphabet = alphabet
		reducedAlphabet.remove(letter)
		replaceLetter = random.choice(reducedAlphabet)
		seq = SNP(seq=seqObj,pos=pos,rl=replaceLetter)
	return seq

def subsample(ref,readError=identity,numReads=10000,readRange = [1000,20000],errorFreq=0.1):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	refLength =  len(ref)
	seqList = []
	for i in range(numReads):
		seqLength = random.randrange(readRange[0],readRange[1])
		start = random.randrange(refLength)
		## randomly subsample from reference
		seq = ref[start:start+seqLength]
		seq = readError(seqObj=seq,freq=errorFreq)
		record=SeqRecord(seq,'fragment_%i' % (i+1),'','')
		seqList.append(record)
	return seqList