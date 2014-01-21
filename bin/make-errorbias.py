#!/usr/bin/env python
## Script for making errorbias files
import csv
import sys
baseprob =  float(sys.argv[1])
bias = float(sys.argv[2])
biasprob = min([baseprob*bias,1.0])

row = ['G','A','T',biasprob,round(0.1 * biasprob,2)]

outfile = '../data/errorbias_%s_%s.dat' % (baseprob,bias)
print outfile
with open(outfile, 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(row)


