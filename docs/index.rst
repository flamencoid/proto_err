.. proto_err documentation master file, created by
   sphinx-quickstart on Mon Dec  9 11:48:02 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to proto_err's documentation!
=====================================



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





Api
================

:mod:`proto_err` -- Main package
****************************

.. automodule:: proto_err

:mod:`proto_err.simulation` -- Simulating Error Reads
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. autoclass:: proto_err.simulation.simulateError
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.errorCount` -- Error Counting
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.errorCount
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.plot` -- Plotting
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.plot
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.query` -- Working with MongoDB
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.query
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.fastaIO` -- Working with Fasta Files
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.fastaIO
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.align` -- Aligning Reads
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.align
   :members:
   :undoc-members:
   :show-inheritance:

:mod:`proto_err.utils` -- Convinience functions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. automodule:: proto_err.utils
   :members:
   :undoc-members:
   :show-inheritance:

	

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

