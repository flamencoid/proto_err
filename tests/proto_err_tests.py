from nose.tools import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from proto_err.errorCount import *
from  proto_err.simulation import simulateError,complexError
from pysam import AlignedRead
from proto_err.fastaIO import *
from proto_err import align
from optparse import Values
from proto_err.query import *
from proto_err.plot import *
import random


def testErrorClassBasic():
    """Just to everything works in the most basis case"""
    a = AlignedRead()
    a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.qual = '++))++)+*)******)))+)**+*+++)**)*+)'
    errorObj = error('A','T',a,0)
    assert_equal(errorObj.true,'A')
    assert_equal(errorObj.emission,'T')
    assert_equal(errorObj.before(2),'NN')
    assert_equal(errorObj.after(2),'GC')
    assert_equal(errorObj.isSnp,True)
    assert_equal(errorObj.isIndel,False)
    assert_equal(errorObj.qual,11)
    assert_equal(errorObj.qscore(2),9)


def testComplexErrorSim():
    """Test the read error simulation"""
    N= 500
    randomLeftFlank = "".join([random.choice(getAlphabet()) for _ in range(N)])
    randomRightFlank = "".join([random.choice(getAlphabet()) for _ in range(N)])

    seq = randomLeftFlank+'AGTATACCTCGCATCGATCGATCG' +randomRightFlank# len 12 
    ref = "".join([random.choice(getAlphabet()) for _ in range(100000)])+seq + "".join([random.choice(getAlphabet()) for _ in range(100000)])
    refRecord=SeqRecord(Seq(ref),'Chromosome dna:chromosome chromosome:ASM19595v1:Chromosome:1:4411532:1','','')
    record=SeqRecord(Seq(seq),'st=%s'%(100000),'','')
    opt = Values()
    rep = [0.1]*N
    errorSim = complexError(record,opt,id = 'st=%s'%(100000),baseErrorProb = rep+[0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15]+rep)
    assert_equal(errorSim.errorProb,rep+[0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15]+rep)
    errorSim.snp(N,'T')
    assert_equal(errorSim.errorProb,rep+[0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15]+rep)
    assert_equal( str(errorSim.seq) , randomLeftFlank+'TGTATACCTCGCATCGATCGATCG'+randomRightFlank)
    errorSim.ins(N+3,'TTT')
    assert_equal(errorSim.errorProb,rep+[0.1,0.15,0.1,0.15,0.15,0.15,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15,0.1,0.15]+rep)
    assert_equal( str(errorSim.seq) , randomLeftFlank+'TGTATTTTACCTCGCATCGATCGATCG'+randomRightFlank)
    errorSim.deletion(N+17,5)
    assert_equal( str(errorSim.seq) , randomLeftFlank+'TGTATTTTACCTCGCATGATCG'+randomRightFlank)

    ## Write these to a fasta file
    refFilename = 'tests/ref.fa'
    readFilename = 'tests/read.fq'
    writeFasta(filename = refFilename,seqList = [refRecord])
    writeFastq(filename = readFilename,seqList = [errorSim.record])

    ## Align ref to read
    align.refIndex(file=refFilename)
    # ## Align reads to the reference
    samfileName = readFilename + '.sam'
    aligned = align.align(reference=refFilename, read_file=readFilename,stdout=samfileName)

    reader = errorReader(samfile=samfileName,ref=ref)
    errorList = []
    for error in reader:
        print error.read.cigar
        assert_equal(error.read.cigarstring,'M%sI3M%sD5M%s'%(N+4,10,N+5))
        errorList.append(error)
    assert_equal(len(errorList),3)  
    error1 = errorList[0]
    assert_equal(error1.isSnp,True)
    assert_equal(error1.errorType,'SNP')
    assert_equal(error1.true,'A')
    assert_equal(error1.emission,'T')
    assert_equal(error1.after(3),'GTA')

    error2 = errorList[1]
    assert_equal(error2.isIndel,True)
    assert_equal(error2.isInsertion,True)
    assert_equal(error2.errorType,'Insertion')
    assert_equal(error2.true,'')
    assert_equal(error2.emission,'TTT') 
    assert_equal(error2.after(3),'TAC')
    assert_equal(error2.before(3),'GTA')  

    error3 = errorList[2]
    assert_equal(error3.isIndel,True)
    assert_equal(error3.isDeletion,True)
    assert_equal(error3.errorType,'Deletion')
    assert_equal(error3.true,'CGATC')
    assert_equal(error3.emission,'')   
    assert_equal(error3.after(3),'GAT')
    assert_equal(error3.before(3),'CAT') 

    assert_equal(error1.alignedCorrectly,True)
    assert_equal(error1.alignedDist,0)

    ## Test counter
    opt.maxKmerLength = 4
    errorCounter = counter(ref,opt,errorList)
    assert_equal(errorCounter.probKmer('A') +  errorCounter.probKmer('T') +
                    errorCounter.probKmer('C')  +errorCounter.probKmer('G') ,1)
    count,errorQueryList = errorCounter.getCount(returnList=True)
    assert_equal(count,3)
    assert_equal(errorCounter.getCount(kmerBefore='GTA'),1)
    assert_equal(errorCounter.getCount(kmerBefore='CAT'),1)
    assert_equal(errorCounter.getCount(kmerAfter='TAC'),1)
    assert_equal(errorCounter.getCount(kmerAfter='GAT'),1)

def testbasicQuery():
    """Test DB stuff"""
    ## testing db stuff
    database = errordb()
    database.deleteAll()
    a = AlignedRead()
    a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.qual = '++))++)+*)******)))+)**+*+++)**)*+)'
    errorObj = error('A','T',a,0)
    errorObj = database.addError(error=errorObj)

def testPlotting():
    """Test the plotting module"""
    r = ro.r
    x = ro.IntVector(range(10))
    y = r.rnorm(10)

    
    grdevices.png(file="tests/img/test.png", width=512, height=512)
    r.layout(r.matrix(ro.IntVector([1,2,3,2]), nrow=2, ncol=2))
    r.plot(r.runif(10), y, xlab="runif", ylab="foo/bar", col="red")
    grdevices.dev_off()

    mtcars = datasets.__rdata__.fetch('mtcars')['mtcars']
    rnorm = stats.rnorm
    dataf_rnorm = ro.DataFrame({'value': rnorm(300, mean=0) + rnorm(100, mean=3),
                                  'other_value': rnorm(300, mean=0) + rnorm(100, mean=3),
                                  'mean': ro.IntVector([0, ]*300 + [3, ] * 100)})
    gp = ggplot2.ggplot(mtcars)

    pp = gp + \
         ggplot2.aes_string(x='wt', y='mpg') + \
         ggplot2.geom_point()
    grdevices.png(file="tests/img/test2.png", width=512, height=512)
    pp.plot()
    grdevices.dev_off()

def testPlottingMore():
    """Test the plotting module"""
    opt = Values()
    opt.imgDir = "tests/img/"
    opt.jsonDir = "tests/json/"
    testPlotter = plotter(opt)
    dic = {'TTT':10,'AAA':5}
    testHistPlotter = histPlotter(dic,opt,filename="testHist")
    testHistPlotter.plot()
    dic = {'AAA':{'Expected':1,'Observed':3},'TTT':{'Expected':10,'Observed':30}}
    testMultiHistPlotter = multiHistPlotter(dic,opt,filename="testMultiHist")
    testMultiHistPlotter.plot()

    





  




