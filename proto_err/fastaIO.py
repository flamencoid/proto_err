from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import *
import random 

def identity(seqObj):
	return Seq("".join(list(seqObj)))

def singleSNP(seqObj,freq):
	"""
	Function to randomly add SNP errors
	"""
	seq = list(seqObj)
	alphabet = getAlphabet()
	if random.random() < freq:
		pos = random.randint(0,len(seq))
		letter =seq[pos]
		reducedAlphabet = alphabet
		reducedAlphabet.remove(letter)
		replaceLetter = random.choice(reducedAlphabet)
		seq[pos] = replaceLetter
	return Seq("".join(seq))

def getRef(filename):
	for seq_record in SeqIO.parse(filename, "fasta"):
	    ref = seq_record.seq
	return ref

def subsample(filename,readError=identity,numReads=10,readRange = [100,200],errorFreq=0.1):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	ref = getRef(filename)
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

def writeFasta(filename,seqList):
	"""
	Function to take list of Seq objects and write a fasta file
	"""
	output_handle = open(filename, "w")
	SeqIO.write(seqList,output_handle,'fasta')
	return 1


		







