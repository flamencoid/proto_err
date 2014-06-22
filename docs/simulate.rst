Simulating Errors
===================

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


.. autoclass:: proto_err.simulation.complexError