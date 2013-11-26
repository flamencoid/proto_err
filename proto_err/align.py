import logging
from Bio.Sequencing import Applications
from bwamemWrapper import BwaMemAlignCommandline
import os

def refIndex(file):
	"""
	Function to generate BWA index
	"""
	if os.path.exists(file + '.bwt'):
		pass
	else:
		logging.info("Creating BW index of reference")
		index_cmd = Applications.BwaIndexCommandline(infile=file, algorithm="bwtsw")
		index_cmd()
	return 1

def align(reference, read_file, stdout,algorithm='bwa-mem'):
	if algorithm=='bwa-mem':
		logging.info("Aligning reads to reference with bwa-mem")
		alignCmd = BwaMemAlignCommandline( reference=reference, read_file=read_file)
	return alignCmd(stdout=stdout)


