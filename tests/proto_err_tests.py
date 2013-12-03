from nose.tools import *
import proto_err.metrics

def testErrorClass():
    true = 'A'
    emission = 'T'
    rightFlank = 'CC'
    leftFlank = 'AT'
    error = proto_err.metrics.error(true,emission,leftFlank,rightFlank)
    assert_equal(error.flankLength,(2,2))
    leftFlank = ''
    error = proto_err.metrics.error(true,emission,leftFlank,rightFlank)
    assert_equal(error.flankLength,(0,2))
    assert_equal(error.trueSeq,'ACC')
    assert_equal(error.leftFlank,'')
    assert_equal(error.before(2),'NN')
    assert_equal(error.before(8),'N'*8)
    assert_equal(error.isSnp,True)
    rightFlank = ''
    error = proto_err.metrics.error(true,emission,leftFlank,rightFlank)
    assert_equal(error.after(2),'NN')
    emission = 'TT'
    error = proto_err.metrics.error(true,emission,leftFlank,rightFlank)
    assert_equal(error.isSnp,False)
    assert_equal(error.isIndel,True)







