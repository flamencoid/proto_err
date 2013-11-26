from Bio import SeqIO
from Bio.Seq import Seq
from random import randrange

def subsample(filename,numReads=1000,readRange = [1000,20000]):
	"""
	Function to take a fasta file subsample reads and generate a list of 
	subsampled reads
	"""
	print "reading " + filename
	for seq_record in SeqIO.parse(filename, "fasta"):
	    ref = seq_record
	refLength =  len(ref)
	seqList = []
	for i in range(numReads):
		seqLength = randrange(readRange[0],readRange[1])
		start = randrange(refLength)
		## randomly subsample from reference
		seqList.append(ref[start:start+seqLength])
	return seqList

def writeFasta(filename,seqList):
	output_handle = open(filename, "w")
	SeqIO.write(seqList,output_handle,'fasta')
	return 1


		







