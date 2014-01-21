.. proto_err documentation master file, created by
   sphinx-quickstart on Mon Dec  9 11:48:02 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to proto_err's documentation!

Intro
================
Requirements
================
Install requirements with

**pip install -r requirements.txt**

**sh nonPiprequirements.sh**

nonPiprequirements installs bwa, R and required R packages. Comment out unneeded lines before running. 

Tests 
================
To run tests run:



**nosetests**

from the top directory. 

Sample Usage
================

Simulating reads with errors
-----------------------------
**python subsample.py -r referenceFile.fa**

subsample.py supports many optional arguments. Run **subsample.py -h** to see options.

Context specific errors
^^^^^^^^^^^^^^^^^^^^^^^^^
In order to change the error frequencey based on a given nucleotide context run.

**subsample.py --errorBiasFile path/to/errorBias.dat**

errorBias.dat is a tab delimited text file with 4 columns

e.g.

kmerBefore  errorBase   kmerAfter   probability

AA \\t A \\t \\t 0.3

changes the probability of an error to occur at the 3rd base of in an AAA nucleotide context. 
the context is searched via regex and can be included in the errorBias file. 

e.g.

A.T \\t \\t \\t 0.2

changes the probability of an error occuring after ATT|AGT|ACT|AAT to 20%. 


Counting errors in a samfile
-----------------------------
To iterate through errors in samfile
::
   from errorCount import errorReader
   for error in errorReader(ref,samfile):
       print error

   SNP error(T to A)
   Deletion error(GC to )

If you want to count the errors in a samfile
::
   from optparse import Values
   from errorCount import counter
   from fastaIO import getRef
   ## Alternatively pass these from OptionParser see https://github.com/flamencoid/proto_err/blob/master/bin/errorStats.py
   opt = Values()
   opt.maxKmerLength = 3 
   opt.outDir = 'results/'
   opt.imgDir = opt.outDir+'img/'
   opt.jsonDir = opt.outDir +' json/'
   ##
   ref = getRef(refFilename)
   ## Initalise a counter object
   errorCounter = counter(ref,opt,samfile=opt.samfile)
   errorCounter.plotHist() ## Plot predefined histograms
   errorCounter.getCount(kmerBefore='A') ## Count errors preceded by 'A'
   errorCounter.getCount(kmerBefore='A',type='SNP') ## Count SNP errors preceded by 'A'

See :class:`~proto_err.errorCount` for more methods