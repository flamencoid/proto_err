Counting errors in Samfile
===========================

To iterate through errors in samfile::
   from errorCount import errorReader
   for error in errorReader(ref,samfile):
       print error

   SNP error(T to A)
   Deletion error(GC to )

If you want to count the errors in a samfile::
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

.. autoclass:: proto_err.errorCount.errorReader
