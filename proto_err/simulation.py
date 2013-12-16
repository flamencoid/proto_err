#!/usr/bin/env python
from utils import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random 
import math
from copy import copy
from error import error
from pysam import AlignedRead
# A python module for simulating errors in reads.
class simulateError():
    """Class of objects to simulate errors in a SeqRecord"""
    def __init__(self,record,opt,id):
        self.alphabet = ['A','T','C','G']
        self.seq = record.seq
        self.id = id
        self.opt = opt
        self.errorProb = [0]*len(self.seq)


    def snp(self,pos,rl):
        """function to induce a SNP"""
        seq = list(self.seq)
        truth = str(seq[pos])
        seq[pos] = rl
        self.seq = Seq("".join(seq))
        self.id += 's%s,%s' % (str(pos),str(rl))
        return error(true=truth,emission=rl,read=self.alignedRead,readPos=pos)

    def ins(self,pos,rl):
        """function to induce a insertion"""
        seq = list(self.seq)
        seq = seq[:pos+1] + list(rl) + seq[pos+1:]
        ## Also need to insert the equivalent qscores
        self.errorProb = self.errorProb[:pos+1] + [self.errorProb[pos] for _ in list(rl)]  +self.errorProb[pos+1:]
        self.seq = Seq("".join(seq))
        self.id += 'i%s,%s' % (str(pos),str(rl))
        return error(true='',emission=rl,read=self.alignedRead,readPos=pos)

    def deletion(self,pos,dlen):
        """function to induce a deletion"""
        seq = list(self.seq)
        newSeq = seq[:pos] + seq[pos+dlen:]
        deletedSeq = "".join(seq[pos:pos+dlen])
        ## Also need to delete the equivalent qscores
        self.errorProb = self.errorProb[:pos] + self.errorProb[pos+dlen:]
        self.seq = Seq("".join(newSeq))
        self.id += 'd%s,%s' % (str(pos),str(dlen))
        return error(true=deletedSeq,emission='',read=self.alignedRead,readPos=pos)

    def indel(self,pos):
        """function to induce an INDEL"""
        indelSize = 0
        while indelSize < 1:
            indelSize = int(round(random.gauss(self.opt.indelMean,self.opt.indelSd)))
        if random.random() < 0.5:
            return self.ins(pos,rl="".join([random.choice(self.alphabet) for _ in range(indelSize)]))
        else:
            return self.deletion(pos,dlen=indelSize)
    @property
    def alignedRead(self):
        """Returns an alignedRead object of a the subsampled read aligned perfectley to ref"""
        aRead = AlignedRead()
        aRead.seq=str(self.seq)
        aRead.qual = self.qscore('ascii')
        return aRead

    @property 
    def record(self):
        """Return a Bio.SeqRecord"""
        return SeqRecord(self.seq,self.id,letter_annotations={'phred_quality':self.qscore()})

    def qscore(self,t='int'):
        qscores =[int(-10 * math.log10(p+0.000001)) for p in self.errorProb]
        if t=='int':
            return qscores
        elif t=='ascii':
            return "".join([intToAscii(i) for i in qscores])
            


class singleSNP(simulateError):
    def __init__(self,record,opt,id):
        simulateError.__init__(self,record,opt,id)
        self.errorProb = [random.gauss(opt.snpFreq, 0.01) for _ in range(len(self.seq))] 

    def error(self):
        """Function to induce and error"""
        ## iterate through the probability list
        for pos,prob in enumerate(self.errorProb):
            if random.random() < prob:
                letter = self.seq[pos]
                alphabet = copy(self.alphabet)
                reducedAlphabet = alphabet
                reducedAlphabet.remove(letter)
                replaceLetter = random.choice(reducedAlphabet)
                self.snp(pos,replaceLetter)

class complexError(simulateError):
    def __init__(self,record,opt,id,baseErrorProb=None,errorBias=None):
        simulateError.__init__(self,record,opt,id)
        self.errorBias = errorBias
        if baseErrorProb:
            self.errorProb = baseErrorProb
        else:
            if self.errorBias:
                self.errorProb = [random.gauss(opt.snpFreq, opt.snpFreqSd) for _ in self.seq]
                ## Now check for regular expression matches
                for pattern,prob in self.errorBias.iteritems():
                    results = pattern.finditer(str(self.seq))
                    for result in results:
                    	# prob is a tuple (+pos to effected base,prob)
                    	self.errorProb[result.start(0)+prob[0]] = random.gauss(prob[1], opt.snpFreqSd) 

            else:
            	self.errorProb = [random.gauss(opt.snpFreq, opt.snpFreqSd) for _ in self.seq]

        assert len(self.errorProb) == len(self.seq)

    def error(self):
        """Function to induce errors"""
        ## iterate through the probability list
        errorList = []
        for pos,prob in enumerate(self.errorProb):
            if random.random() < prob and pos < len(self.errorProb):
                letter = self.seq[pos]
                alphabet = copy(self.alphabet)
                reducedAlphabet = alphabet
                reducedAlphabet.remove(letter)
                replaceLetter = random.choice(reducedAlphabet)
                if random.random() <= self.opt.SnpIndelRatio:
                    e = self.snp(pos,replaceLetter)
                    errorList.append(e)
                else:
                    e = self.indel(pos)
                    errorList.append( e)
        return errorList
            

