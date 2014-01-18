#!/usr/bin/env python
from utils import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random 
import math
from copy import copy
from error import error
from pysam import AlignedRead
import itertools
random.seed(1)
# A python module for simulating errors in reads.
class simulateError():
    """Class of objects to simulate errors in a SeqRecord"""
    def __init__(self,record,opt,id):
        self.alphabet = ['A','T','C','G']
        seq = str(record.seq)
        self.read = list(seq)
        self.ref = list(seq)
        self.id = id
        self.opt = opt
        self.errorProb = [0]*len(seq)


    def snp(self,pos,rl):
        """function to induce a SNP"""
        ## change the base in the read, not in ref
        old = self.read[pos]
        self.read[pos] = rl
        self.pos += 1

    def ins(self,pos,rl):
        """function to induce a insertion"""
        ## add bases in the read, add - * len(bases) in ref
        self.read = self.read[:pos+1] + list(rl) + self.read[pos+1:]
        self.ref = self.ref[:pos+1] + ['_']*len(rl) + self.ref[pos+1:]
        ## Also need to insert the equivalent qscores
        self.errorProb = self.errorProb[:pos+1] + [self.opt.snpFreq for _ in list(rl)]  +self.errorProb[pos+1:]
        # self.id += 'i%s,%s' % (str(pos),str(rl))
        self.pos += len(rl) + 1 # need to jump pos to avoid putting errors in inserted seq
        ## We haven't moved along the reference at all
        

    def deletion(self,pos,dlen):
        """function to induce a deletion"""
        ## add  - * len(del) in the self.read, don't edit but move along the reference
        self.read = self.read[:pos] + ['_']*dlen+self.read[pos+dlen:]
        ## Also need to delete the equivalent qscores
        self.errorProb = self.errorProb[:pos] + ['_']*dlen+self.errorProb[pos+dlen:]
        self.pos += dlen

        

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
        aRead.qname = self.id
        aRead.is_reverse = self.opt.is_reverse
        ## Add the start position of where the read should align to
        # aRead.positions = [] 
        return aRead
    @property 
    def seq(self):
        return Seq("".join(self.read).replace('_',''))

    @property 
    def record(self):
        """Return a Bio.SeqRecord"""
        return SeqRecord(self.seq,self.id,letter_annotations={'phred_quality':self.qscore(t='int',aligned=False)})
    @property 
    def SNPProb(self):
        "Convert errorProb into SNPProb list"
        out = []
        for p in self.errorProb:
            if not p == '_':
                out.append(self.opt.SnpIndelRatio * p)
            else:
                out.append('_')
        return out

    def qscore(self,t='int',aligned=True):
        "probability of a SNP only"
        qscores = []
        for p in self.SNPProb:
            if not p == '_':
                if p < 0:
                    p = 0
                qscores.append(int(-10 * math.log10(p+0.000001)))
            elif aligned:
                qscores.append('_')

        if t=='int':
            return qscores
        elif t=='ascii':
            return "".join([intToAscii(p) for p in qscores if not p == '_']  )
            


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
        self.pos = 0
        if baseErrorProb:
            self.errorProb = baseErrorProb
        else:
            if self.errorBias:
                self.errorProb = [random.gauss(opt.snpFreq, opt.snpFreqSd) for _ in self.read]
                ## Now check for regular expression matches
                for pattern,prob in self.errorBias.iteritems():
                    results = pattern.finditer(str(self.read))
                    for result in results:
                    	# prob is a tuple (+pos to effected base,prob,sd)
                    	self.errorProb[result.start(0)+prob[0]] = random.gauss(prob[1], prob[2]) 

            else:
            	self.errorProb = [random.gauss(opt.snpFreq, opt.snpFreqSd) for _ in self.read]

        assert len(self.errorProb) == len(self.read)

    def error(self):
        """Function to induce errors"""
        ## iterate through the probability list
        while self.pos < len(self.errorProb):
            assert len(self.read) == len(self.ref)
            assert len(self.read) == len(self.errorProb)
            r = random.random()
            prob = self.errorProb[self.pos]
            if r < prob:
                letter = self.ref[self.pos]
                alphabet = copy(self.alphabet)
                reducedAlphabet = alphabet
                reducedAlphabet.remove(letter)
                replaceLetter = random.choice(reducedAlphabet)
                if random.random() <= self.opt.SnpIndelRatio:
                    self.snp(self.pos,replaceLetter)
                else:
                    e = self.indel(self.pos)
            else:
                self.pos += 1

        ## Create error objects and upload to database
        errorList = []
        readLength = len(self.seq)
        il = zip(self.read,self.ref,self.qscore('int'))
        iterator = itertools.cycle(il)
        i = 0
        refPos = self.opt.refPos 
        readPos = 0
        base, truth, qscore = iterator.next()
        # print record.
        while i < len(il) and readPos < len(self.seq):
            if base == '_':
                deletedBases = []
                ## Deletion
                while base == '_':
                    deletedBases.append(truth)
                    i+=1
                    refPos +=1
                    base, truth, qscore = iterator.next()
                errorList.append(error(true="".join(deletedBases),
                                        emission="",read=self.alignedRead,
                                        readPos=readPos,readLength=readLength,
                                        refPos=refPos-len(deletedBases) ))
            elif truth == '_':
                ## Insertion
                insertedBases = []
                while truth == '_':
                    insertedBases.append(base)
                    i+=1
                    readPos +=1
                    base, truth, qscore = iterator.next()
                errorList.append(error(true="",
                                        emission="".join(insertedBases),
                                        read=self.alignedRead,
                                        readPos=readPos-len(insertedBases),
                                        readLength=readLength,refPos=refPos) )
            else: 
                ## M
                if base != truth:
                    ## SNP
                    errorList.append(error(true=truth,emission=base,
                                            read=self.alignedRead,readPos=readPos,
                                            readLength=readLength,refPos=refPos) )
                base, truth, qscore = iterator.next()
                i +=1
                refPos += 1
                readPos += 1
        return errorList
                    



            

        # for pos,prob in enumerate(self.errorProb):
        #     if  < prob and pos < len(self.errorProb):

                
            

