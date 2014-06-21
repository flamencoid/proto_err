[Proto_err](http://proto-err.readthedocs.org/)
=========

[See Detailed Documenation](http://proto-err.readthedocs.org/)
##Requirements##
Install requirements with

*pip install -r requirements.txt*

*sh nonPiprequirements.sh*

nonPiprequirements installs bwa, R and required R packages. Comment out unneeded lines before running. 
### Tests ###
To run tests run:

**nosetests**

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




