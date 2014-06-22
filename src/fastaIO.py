#!/usr/bin/env python
## Handling fasta files
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def getRef(filename):
	for seq_record in SeqIO.parse(filename, "fasta"):
	    ref = seq_record.seq
	return ref

def writeFasta(filename,seqList):
	"""
	Function to take list of Seq objects and write a fasta file
	"""
	output_handle = open(filename, "w")
	SeqIO.write(seqList,output_handle,'fasta')
	return 1
def writeFastq(filename,seqList):
	"""
	Function to take list of Seq objects and write a fasta file
	"""
	output_handle = open(filename, "w")
	SeqIO.write(seqList,output_handle,'fastq')
	return 1


		







