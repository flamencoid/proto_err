import itertools

def getAlphabet():
	return ['A','T','C','G']

def kmerCombo(a,r):
	"""
	Function to return all combinations of kmer length i in alphabet
	"""
	return ["".join(i) for i in itertools.product(a,repeat=r)]

