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
    

    def __init__(self,true,emission,read,readPos,readLength=None):
        self.true = true
        self.emission = emission
        assert true != '' or emission != ''
        assert true != emission
        self.read = read
        self.readPos = readPos # position on read where error starts
        if not self.readLength:
            self.readLength = len(self.read.seq)
        self.readPer = float(readPos) / float(self.readLength)

        ## If the read is on the reverse stand then 
        if self.read.is_reverse:
            self.true = reverse_complement(self.true)
            self.emission = reverse_complement(self.emission)


        ## Was the read aligned correctly? 
        try:
            ## messy way of extracting read length
            # s = list(self.read.qname.split('st=')[1])
            # curInt = '0'
            # tempList = []
            # while curInt.isdigit():
            #     curInt = s.pop(0)
            #     tempList.append(curInt)
            sampPos = int(self.read.qname.split('st=')[1])
            # sampPos = int("".join(tempList[:-1]))
            self.alignedDist =  int(abs(sampPos - self.read.positions[0]))
        except:
            self.alignedDist = None
        if self.alignedDist is None:
            self.alignedCorrectly = None
        elif not self.alignedDist is None and (self.alignedDist < len(self.read.seq)):
            self.alignedCorrectly = True
        else:
            self.alignedCorrectly = False

    def __str__(self):
        return "%s error(%s to %s)" %(self.errorType,self.true,self.emission)

    def before(self,j):
        """Return the preceding j bases,return N when bases missing"""
        i = self.readPos
        b = self.read.seq[i-j:i]
        while len(b) < j:
            b = 'N' +b
        return b
    def after(self,j):
        """Return the following j bases,return N when bases missing"""
        if self.isIndel:
            ## For insertions, we need to move along the read otherwise the 
            ## "after" bases will include the inserted sequence
            i = self.readPos + len(self.emission) -1 
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
            return None
        else:
            # print self.read.qqual
            return asciiToInt(self.read.qqual[self.readPos])
    @property 
    def tlen(self):
        """Return the length of the error"""
        return  int(max([len(self.emission),len(self.true)]))
    @property 
    def doc(self):
        """Return a pymongo document"""
        # return {'true':self.true,'emission':self.emission,'read':str(self.read.seq),
        #         'readPos':self.readPos,'readPer':self.readPer,'alignedDist':self.alignedDist,
        #         'leftFlank':self.before(10),'rightFlank':self.after(10),'type':self.errorType,
        #         'qual':self.qual,'tlen' :self.tlen}
        return {'true':self.true,'emission':self.emission,
                'readPos':self.readPos,'readPer':self.readPer,'alignedDist':self.alignedDist,
                'leftFlank':self.before(10),'rightFlank':self.after(10),'type':self.errorType,
                'qual':self.qual,'tlen' :self.tlen,'readLength':self.readLength}        