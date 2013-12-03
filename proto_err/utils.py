import itertools

def getAlphabet():
	return ['A','T','C','G']

def kmerCombo(a,r):
	"""
	Function to return all combinations of kmer length i in alphabet
	"""
	return ["".join(i) for i in itertools.product(a,repeat=r)]
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

