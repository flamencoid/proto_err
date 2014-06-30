from utils import *
class error():
    """ 
    Information about the errors in a read 

    Attributes
    ----------
    true : string
        reference bases
    emission : string 
        Read bases
    read : AlignedRead
        pysam AlignedRead object from which error was derived
    readPos : int
        position along read where error occured
    readPer : float
        position along read where error occured / length of read sequence
    alignedDist : int or None
        number of bases between position of start of alignment and where read was sampled None if sampled position n/a 
    alignedCorrectly : bool or None
        If the alignedDist is less then the read length True else False None is alignedDist is None
    errorType : string
        type of error (SNP,Insertion,Deletion)

    Parameters
    ----------
    true: string
        reference base(s)

    emission: string
        read base(s)

    read: AlignedRead
         pysam AlignedRead object from which error was derived

    readPos: int
        position along read where error occured / length of read sequence



    See Also
    --------
    counter

    Examples
    --------
    >>> from pysam import AlignedRead
    >>> from import proto_err.errorCount import error
    >>> read = AlignedRead()
    >>> read.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    >>> read.qual = '++))++)+*)******)))+)**+*+++)**)*+)'
    >>> error = proto_err.errorCount.error('T','A',read,0)
    >>> error.isSnp
    True
    >>> error.isIndel
    False
    >>> error.qual
    11
    >>> error.before(2)
    NN 
    >>> error.after(2)
    GC 
    >>> error.qscore(0) ## equivalent to error.qual
    11
    >>> error.qscore(2)
    9

    """
    

    def __init__(self,true,emission,read,readPos,refPos,readLength=None,):
        self.true = true
        self.emission = emission
        assert true != '' or emission != ''
        assert true != emission


        self.read = read
        self.readPos = readPos # position on read where error starts
        if readLength:
            self.readLength = readLength
        else:
            self.readLength = self.read.rlen
        assert self.readLength != None

        self.readID = self.read.qname
        # sampPos = int(self.read.qname.split('st=')[1].split('&')[0])

            
        self.readPer = float(readPos) / float(self.readLength)

        ## What are the coordinates of the error?
        self.refPos = refPos
        ## Was the read aligned correctly? 
        try:
            if self.read.cigar[0][0] in [4,5]:
                    clippedBases = self.read.cigar[0][1]
            else:
                clippedBases = 0
            # self.mappedDist =  abs(int(abs(sampPos - self.read.positions[0])) -clippedBases)
            self.mappedDist = None
        except:
            self.mappedDist = None
        if self.mappedDist is None:
            self.mappedCorrectly = None
        elif not self.mappedDist is None and (self.mappedDist < self.readLength):
            self.mappedCorrectly = 1
        else:
            self.mappedCorrectly = 0

        assert self.readPos >=0
        assert self.readPos < self.readLength


    def __eq__(self, other):
        "Checks if two errors are equivalent"
        if isinstance(other, self.__class__):
            refBase = self.true == other.true
            emitBase = self.emission == other.emission
            readPosBool = self.readPos == other.readPos
            coorBoole = self.refPos == other.refPos
            ## add check for reference coordinates
            return refBase and emitBase and readPosBool and coorBoole
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "%s error(%s to %s)" %(self.errorType,self.true,self.emission)

    def before(self,j):
        """Return the preceding j bases,return N when bases missing"""
        if self.isInsertion:
            ## For insertions, we need to move along the read otherwise the 
            ## "after" bases will include the inserted sequence
            i = self.readPos +1 
        else:
            i = self.readPos

        # i = self.readPos
        b = self.read.seq[i-j:i]
        while len(b) < j:
            b = 'N' +b
        return b
    def after(self,j):
        """Return the following j bases,return N when bases missing"""
        if self.isInsertion:
            ## For insertions, we need to move along the read otherwise the 
            ## "after" bases will include the inserted sequence
            i = self.readPos + len(self.emission) 
        elif self.isDeletion:
            i = self.readPos -1
        else:
            i = self.readPos
        
        a = self.read.seq[i+1:j+i+1]
        while len(a) < j:
            a = a + 'N'
        return a

    def qscore(self,i):
        """Return the quality score at a base +i i from the error start position
        error.qscore(0) is equivalent to error.qual """
        return asciiToInt(self.read.qqual[self.readPos+i])

    @property 
    def errorType(self):
        """The type of read error SNP, Insertion or Deletion"""
        if self.isSnp:
            return 'SNP'
        elif self.isInsertion:
            return 'Insertion'
        elif self.isDeletion:
            return 'Deletion'
        else:
            return 'Unknown Error Type'
    @property 
    def isSnp(self):
        """Return whether or not the error is a SNP"""
        return len(self.true) == len(self.emission)
    @property 
    def isIndel(self):
        """Return whether or not the error is a INDEL"""
        return len(self.true) != len(self.emission)
    @property 
    def isInsertion(self):
        """Return whether or not the error is an insertion"""
        return len(self.true) < len(self.emission)
    @property 
    def isDeletion(self):
        """Return whether or not the error is a deletion"""
        return len(self.true) > len(self.emission)
    @property
    def qual(self):
        """Return quality score of the base where the error occured"""
        if self.isDeletion:
            ## deletions deleted the base so there won't be a qual score at the deleted position
            ## return the previous position instead ?? 
            return asciiToInt(self.read.qqual[self.readPos-1])
        elif self.isInsertion:
            ## If it's an insertion return the mean quality across the indel + previous base
            try:
                qualList = [asciiToInt(self.read.qqual[i]) for i in range(self.readPos-1,self.readPos + self.tlen)]
            except IndexError:
                print "Index Error calculating insertion quality score"
                qualList = [asciiToInt(self.read.qqual[ self.readPos-1 ] )]
            meanqual = float(sum(qualList)) / float(len(qualList))
        else:
            return asciiToInt(self.read.qqual[self.readPos])
    @property 
    def tlen(self):
        """Return the length of the error"""
        return  int(max([len(self.emission),len(self.true)]))
    @property 
    def doc(self):
        """Return a pymongo document"""
        return {'true':self.true,'emission':self.emission,
                'readPos':self.readPos,'readPer':self.readPer,'mappedDist':self.mappedDist,
                'mappedCorrectly': self.mappedCorrectly,
                'leftFlank':self.before(10),'rightFlank':self.after(10),'type':self.errorType,
                'qual':self.qual,'tlen' :self.tlen,'readLength':self.readLength,'refPos': self.refPos,
                'readID':self.readID,'strand': 'reverse' if self.read.is_reverse else 'forward'}        