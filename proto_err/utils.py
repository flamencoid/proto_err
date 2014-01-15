import itertools
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import math
def getAlphabet():
    return ['A','T','C','G']

def kmerCombo(r):
    """
    Function to return all combinations of kmer length r in alphabet
    """
    a = getAlphabet()


    return ["".join(i) for i in itertools.product(a,repeat=r)]
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
def asciiToInt(s):
    return ord(s) - 33
def intToAscii(i):
    return  chr(i + 33 )
def popLong(l,i,j):
    """Pops along multiple index"""
    return [l.pop(0) for _ in l[i:j]]
def getDBError():
    return """ ########################################\n####### Is a pymongo instance running at 27017? ##### \n####### Install mongodb and run mongod from shell ####\n########################################"""

def getKeysFromValuesObject(opt):
    """Function to extract keys from a Values object"""
    return [key for key in dir(opt) if not key in ['__cmp__', '__doc__', '__init__', '__module__', '__repr__', '__str__', '_update', '_update_careful', '_update_loose', 'ensure_value', 'read_file', 'read_module']]
# def cigarIntToString(i):
#     dic = {0:'M',1:'I',''}
def reverse_complement(string):
    record  = SeqRecord(Seq(string),'','','')
    return str(record.reverse_complement().seq)

def probToQscore(p):
    return int(-10 * math.log10(p+0.000001))


