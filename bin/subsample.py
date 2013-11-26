#!/usr/bin/env python

## python script that will take a fasta file and return a fasta file 
## of subsampled reads
import sys
import os
sys.path.insert(0, os.path.abspath('../proto_err'))

from fastaIO import subsample,writeFasta
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="fasta input file", metavar="FILE")
(opt, args) = parser.parse_args()

writeFasta(filename = opt.filename[:-3] + '.subsampled.fa' ,seqList = subsample(filename=opt.filename))
