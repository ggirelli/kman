
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
'''

import oligo_melting as om

class Sequence(om.Sequence):
	"""docstring for Sequence"""

	def __init__(self, seq, t, name = None):
		super(Sequence, self).__init__(seq, t, name)

	def kmers(self, k):
		"""Extract k-mers from Sequence.
		Args:
			k {int} -- substring length
		Returns:
			generator -- kmer generator
		"""
		return self.kmerator(self.text, k, self.natype)

	@staticmethod
	def kmerator(seq, k, t, name = "ref"):
		"""Extract k-mers from seq.
		Args:
			seq {string} -- input sequence
			k {int} -- substring length
		Returns:
			generator -- kmer generator
		"""
		return (KMer(name, i, i+k, seq[i:i+k], t)
			for i in range(len(seq)-k+1))

class KMer(Sequence):
	"""docstring for KMer"""

	def __init__(self, chrom, start, end, seq, t = om.NATYPES.DNA):
		super(KMer, self).__init__(seq, t)
		assert start >= 0
		assert end >= 0
		self.__chrom = chrom
		self.__start = start
		self.__end = end

	@property
	def header(self):
		return "%s:%d-%d" % (self.__chrom, self.__start, self.__end)
	@property
	def seq(self):
		return self.text

	def as_fasta(self):
		return ">%s\n%s\n" % (self.header, self.seq)
	
	def is_ab_checked(self):
		"""Check if AB is fully respected.
		
		Checks if the KMer sequence respects the AB. I.e., it does not contain
		any foreign characters.
		
		Returns:
			bool -- whether AB is respected.
		"""
		for c in set(self.text):
			if not c in self.ab[0]:
				return False
		return True
