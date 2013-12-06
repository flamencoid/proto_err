from nose.tools import *
from Bio.SeqRecord import SeqRecord
import proto_err.metrics
from  proto_err.error import simulateError

def testErrorClassBasic():
	"""Just to everything works in the most basis case"""
	error = proto_err.metrics.error('A','T')
	assert_equal(error.true,'A')
	assert_equal(error.emission,'T')
	assert_equal(error.leftFlank,'')
	assert_equal(error.rightFlank,'')
	assert_equal(error.before(2),'NN')
	assert_equal(error.after(2),'NN')
	assert_equal(error.trueSeq,'A')
	assert_equal(error.emissionSeq,'T')
	assert_equal(error.flankLength,(0,0))
	assert_equal(error.isSnp,True)
	assert_equal(error.isIndel,False)


# def testErrorClass():
#     true = 'A'
#     emission = 'T'
#     rightFlank = 'CC'
#     leftFlank = 'AT'
#     error = proto_err.metrics.error(true,emission,read)
#     assert_equal(error.flankLength,(2,2))
#     leftFlank = ''
#     error = proto_err.metrics.error(true,emission,read)
#     assert_equal(error.flankLength,(0,2))
#     assert_equal(error.trueSeq,'ACC')
#     assert_equal(error.leftFlank,'')
#     assert_equal(error.before(2),'NN')
#     assert_equal(error.before(8),'N'*8)
#     assert_equal(error.isSnp,True)
#     rightFlank = ''
#     error = proto_err.metrics.error(true,emission,read)
#     assert_equal(error.after(2),'NN')
#     emission = 'TT'
#     error = proto_err.metrics.error(true,emission,read)
#     assert_equal(error.isSnp,False)
#     assert_equal(error.isIndel,True)

def testErrorSim():
    """Test the read error simulation"""
    seq = 'ATCGATCGATCG'
    record=SeqRecord(seq,'fragment_id')
    opt = {}
    errorSim = simulateError(record,opt,id = 'fragment_id')
    errorSim.snp(0,'T')
    print errorSim.record.seq
    assert_equal( str(errorSim.seq) , 'TTCGATCGATCG')
    errorSim.ins(3,'TTT')
    assert_equal( str(errorSim.seq) , 'TTCGTTTATCGATCG')
    errorSim.deletion(4,5)
    assert_equal( str(errorSim.seq) , 'TTCGCGATCG')

# def testCounter():
#     print errorCounter.getCount()
#     print errorCounter.getCount(kmer='A')
#     print errorCounter.getCount(truth='A',kmer='A')
#     print errorCounter.getCount(truth='A',emission='T', kmer='A')
#     print errorCounter.getCount(emission='A')
#     print errorCounter.getCount(emission='A',kmer='A')



