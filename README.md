[Proto_err](http://proto-err.readthedocs.org/)
=========

[See Detailed Documenation](http://proto-err.readthedocs.org/)
##Requirements##
Install requirements with

*sudo pip install -r requirements.txt*
### Tests ###
To run tests run:

*nosetests*

from the top directory. 
## Simulating Reads ##

To read in a reference file subsample and realign with random read errors and output a samfile

**python subsample.py -r referenceFile.fa**

or simply,

**./subsample.py** to run with default reference

Run **subsample.py -h** for additional options
```python
Options:
  -h, --help            show this help message and exit
  -r REFFILENAME, --ref=REFFILENAME
                        fasta input ref file
  --numReads=NUMREADS   Number of reads to
                        sample from referecnce (Optional defaults to 1000)
  --meanReadLength=READMEAN
                        Read length is
                        sampled from a normal of this mean (Optional defaults
                        to
                        1000)
  --errorFreqMean=SNPFREQ
                        Probablity of an error
                        at a base occurs is sampled from a normal with mean
                        of this Probablity (Optional defaults to
                        .1)
  --errorFreqSD=SNPFREQSD
                        Probablity of an error
                        at a base occurs is sampled from a normal with SD
                        of this value (Optional defaults to
                        errorFreqMean/10)
  --strandBias=STRANDBIAS
                        Ratio of reads
                        sampled from forward or reverse strand. 1 samples only
                        from                                          forward,
                        0 only from reverse.(Optional defaults to .5)
  --ReadLengthSD=READSD
                        Read length is
                        sampled from a normal with this SD (Optional defaults
                        to
                        meanReadLength/3)
  --SnpIndelRatio=SNPINDELRATIO
                        Ratio of SNP
                        errors to INDEL errors (Optional defaults to 0.5)
                        value of 1 gives ONLY SNP errors, 0 only INDELs
  --meanIndelSize=INDELMEAN
                        INDEL size is
                        sampled from a normal of this mean (Optional defaults
                        to                                                  5)
  --IndelSizeSD=INDELSD
                        INDEL size is
                        sampled from a normal with this SD (Optional defaults
                        to
                        meanIndelSize/2)
```

# Error Reading # 

**python errorStats.py -r referenceFile.fa -s samfile.sam**

## Error Reader ## 
To iterate through errors in samfile
```python
from metrics import errorReader
for error in errorReader:
    print error

SNP error(T to A)
Deletion error(GC to )
...
```

## Error Class##

[See Error Class](https://github.com/flamencoid/proto_err/blob/master/docs/_build/text/index.txt)

## Counting  ##

The counter object will take a list of errors or a samfle and a reference.
errorCounter = counter(ref,samfile=samfile)
or
errorCounter = counter(ref,errorList=errorList)

To count errorTypes
#### Count all the errors    
print errorCounter.getCount()

####   Count all the errors preceded by kmer 'A'
print errorCounter.getCount(kmer='AA')

####  Count all the errors followed by kmer 'AA'
print errorCounter.getCount(kmer='AA',after=True)

####     Count all the errors for truth 'A' preceded by an A
 print errorCounter.getCount(truth='A',kmer='A')

####     Count all the errors A->T preceded by A
 print errorCounter.getCount(truth='A',emission='T', kmer='A')

####  Count all the errors where the emmited base is 'A' preceded by any kmer
   print errorCounter.getCount(emission='A')

# To return dictonaries of all kmer counts within the reference of length 1,2,3
dic = errorCounter.countRefKmer(maxKmerLength = 3)
dic['ATC']
12345

dic = errorCounter.countErrorKmer(maxKmerLength = 3)
Will return a the emmison matrix (as a dict) and a dictonary of kmer counts before and after an error
emmisionMatrix,kmerCounts

The kmerCounts is of the form 
kmerCounts['before'][trueBase][EmittedBase][kmerString] and
kmerCounts['after'][trueBase][EmittedBase][kmerString]
e.g. to get the counts of the number of times 'TTT' appears before the A->T transisition
kmerCounts['before']['A']['T']['TTT']
12345







