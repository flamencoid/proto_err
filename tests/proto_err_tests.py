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
        assert_equal(error.read.cigarstring,'M%sI3M%sD5M%s'%(N+4,10,N+5))
        errorList.append(error)
    assert_equal(len(errorList),3)  
    error1 = errorList[0]
    assert_equal(error1.isSnp,True)
    assert_equal(error1.errorType,'SNP')
    assert_equal(error1.true,'A')
    assert_equal(error1.emission,'T')

    error2 = errorList[1]
    assert_equal(error2.isIndel,True)
    assert_equal(error2.isInsertion,True)
    assert_equal(error2.errorType,'Insertion')
    assert_equal(error2.true,'')
    assert_equal(error2.emission,'TTT')  

    error3 = errorList[2]
    assert_equal(error3.isIndel,True)
    assert_equal(error3.isDeletion,True)
    assert_equal(error3.errorType,'Deletion')
    assert_equal(error3.true,'CGATC')
    assert_equal(error3.emission,'')   

    assert_equal(error1.alignedCorrectly,True)
    assert_equal(error1.alignedDist,0)

    ## Test counter
    opt.maxKmerLength = 4
    errorCounter = counter(ref,opt,errorList)
    assert_equal(errorCounter.getCount(),12)

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
    print database.find_one()
    print errorObj.dbID

  




