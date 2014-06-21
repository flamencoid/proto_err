import itertools
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import math
import csv
import re
import numpy as np
import logging
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
    if not s =='_':
        return ord(s) - 33
    else:
        return '_'
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
def qscoreToProb(q):
    return 10 ** (-float(q)/10)

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


def makeGeneToDrugTable(filepath = "../data/drugGene.csv"):
    geneToDrug = AutoVivification()
    with open(filepath,"U") as drugGeneFile:
        csvReader = csv.reader(drugGeneFile,delimiter=",")
        for row in csvReader:
            gene = row[1].strip()
            try:
                geneToDrug[row[1]].append(row[0].strip())
            except:
                geneToDrug[row[1]] = [row[0].strip()]
    return geneToDrug

def getSampleIDList(filepath):
    """Returns a list of all sampleID"""
    with open(filepath,'r') as infile:
        csvReader = csv.reader(infile)
        sample_id_list = []
        for i,row in enumerate(csvReader):
            sample_id = row[0].split('/')[-1].split('.kmer')[0]
            sample_id_list.append(sample_id)
    return sample_id_list

def sampleIDToTruthGenotype(filepath = "../data/FinalTableGenotypeLookup.csv"):
    """Get the genotype calls for each sample from the paper"""
    variantPattern = re.compile("[a-zA-Z]\d{1,4}[a-zA-Z]") 
    sampleIDtoTruthGenotypeDic = {}
    with open(filepath,"U") as drugGeneFile:
        csvReader = csv.reader(drugGeneFile,delimiter=",")
        for i,inrow in enumerate(csvReader):
            if i == 0:
                header = [gene.replace("_aa","") for gene in inrow]
            else:
                variantsOut = []
                for geneIndex,variants in enumerate(inrow):
                    for variant in variants.split(" "):
                            reMatch = variantPattern.search(variant)
                            if reMatch:
                                variantsOut.append(header[geneIndex]+"_"+reMatch.group(0))
                if inrow[0]:
                    sampleIDtoTruthGenotypeDic[inrow[0]] = variantsOut
    return sampleIDtoTruthGenotypeDic

def getSamplePhenotype(filepath="../data/FinalTableLookup.csv"):
    """Returns a dictonary mapping sample_id to drug resistance"""
    idToDrugResistant= AutoVivification()
    with open(filepath,"U") as drugGeneFile:
        csvReader = csv.reader(drugGeneFile,delimiter=",")
        header = csvReader.next()
        csvReader.next()
        for inrow in csvReader:
            drugs = {}
            ## Replace nt with ""
            row = []
            for rs in inrow:
            	if rs == "nt":
            		row.append("")
            	elif rs == "DR":
            		row.append("R")
            	elif rs == "DS":
            		row.append("S")
            	else:
            		row.append(rs)

            identifier =  row[0]
            for ii in range(1,len(row)):
                drugs[header[ii]] = row[ii]
            idToDrugResistant[identifier] = drugs
    return idToDrugResistant
drugList = ["Gentamicin","Penicillin","Trimethaprim","Erythromycin","Methicillin","Ciprofloxacin","Rifampicin","Tetracycline","Vancomycin","Mupirocin","Fusidic","Clindamycin"]
def makeSubToDrugMap():
    subToDrug = AutoVivification()
    with open('../data/subToDrug.txt','U') as subtodrugfile:
        csvReader = csv.reader(subtodrugfile,delimiter="\t")
        header = csvReader.next()
        for row in csvReader:
            drug =  row[0]
            gene = row[1]
            sublist = [i.replace('*','').strip() for i in row[2].split(',')]
            for sub in sublist:
                subToDrug[gene][sub] = drug
    return subToDrug


def cor_idToBlastIDTable(cortexConf):
    cor_idToid = AutoVivification()
    id_stringToCor = AutoVivification()
    with open(cortexConf,'r') as infile:
        csvReader = csv.reader(infile)
        print cortexConf
        for i,row in enumerate(csvReader):
            id_string = row[0].split('/')[-1].split('.kmer')[0]
            cor_idToid[i+1] = id_string
            id_stringToCor[id_string] = i+1
    return cor_idToid,id_stringToCor


def intersect(a, b):
     return list(set(a) & set(b))

def extractGeneNameFromFastaName(fastaname_string):
    return  fastaname_string.split('_col')[0]

def get_gene_presence_list(fastafile):
    gene_presence_list = []
    coverage_grep = re.compile(".*kmer_coverages")
    for panal_fasta_record in SeqIO.parse(fastafile, "fasta"):
        if coverage_grep.match(panal_fasta_record.name):
            gene = extractGeneNameFromFastaName(panal_fasta_record.name)
            gene_presence_list.append(gene)
    return gene_presence_list

def named_nparray_to_nparray(namedarray,dtype=np.float64):
    return namedarray.view(dtype).reshape(namedarray.shape + (-1,))

def getComboMuts():
    listOFComboMuts = {'grlB_E422D':'gyrA_S84L+grlA_V41G+grlA_S80F|grlA_S80F+grlA_I45M',
                    'grlA_I45M':'grlA_S80F+grlA_V41G',
                    "gyrA_S84L":'grlA_S80F|grlB_A116E|grlA_E84K',
                    "grlA_E84G":"gyrA_S84L+grlA_S80F|grlA_S80Y+gyrA_S84L"}

    return listOFComboMuts

