
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
'''

from enum import Enum
import oligo_melting as om

class Sequence(om.Sequence):
	"""docstring for Sequence"""

	def __init__(self, seq, t, name = None):
		super().__init__(seq, t, name)

	def kmers(self, k):
		"""Extract k-mers from Sequence.
		Args:
			k {int} -- substring length
		Returns:
			generator -- kmer generator
		"""
		return self.kmerator(self.text, k, self.natype, self.name)

	def batches(self, k, batchSize):
		"""Split the sequence in batches ready to be fed to the kmerator.
		
		This is useful when parallelizing, due to the fact that generators
		cannot be pickled and can cause problems when using a process-based
		back-end. Instead, pass the sequences to the processes and run kmerator
		there.
		
		Arguments:
			seq {string} -- sequence to be batched
			k {int} -- length of substrings for kmerator
			batchSize {int} -- number of kmers per batch
		"""
		return self.batcher(self.text, k, batchSize)

	def kmers_batched(self, k, batchSize = 1):
		"""Extract batches of k-mers from Sequence.
		
		Arguments:
			k {int} -- substring length

		Keyword Arguments:
			batchSize {number} -- number of k-mers per batch (default: {1})

		Returns:
			generator -- k-mer batch generator
		"""
		assert batchSize >= 1
		if batchSize == 1:
			return self.kmers(k)
		else:
			return self.kmerator_batched(self.text, k, self.natype,
				batchSize, self.name)

	@staticmethod
	def kmerator(seq, k, t, prefix = "ref", prefix_start = 0):
		"""Extract k-mers from seq.
		
		Arguments:
			seq {string} -- input sequence
			k {int} -- substring length
			t {om.NATYPES} -- nucleic acid type
		
		Keyword Arguments:
			prefix {str} -- reference record name (default: {"ref"})
			prefix_start {number} -- if this is a batch, current location for
			                         shifting (default: {0})
		"""
		return (KMer(prefix, i+prefix_start, i+prefix_start+k, seq[i:i+k], t)
				for i in range(len(seq)-k+1))

	@staticmethod
	def batcher(seq, k, batchSize):
		"""Split the sequence in batches ready to be fed to the kmerator.
		
		This is useful when parallelizing, due to the fact that generators
		cannot be pickled and can cause problems when using a process-based
		back-end. Instead, pass the sequences to the processes and run kmerator
		there.
		
		Arguments:
			seq {string} -- sequence to be batched
			k {int} -- length of substrings for kmerator
			batchSize {int} -- number of kmers per batch
		"""
		for i in range(0, len(seq)-k+1, batchSize):
			yield (seq[max(i, i-k+1):min(len(seq)-k+1, i+batchSize)], i)

	@staticmethod
	def kmerator_batched(seq, k, t, batchSize = 1, prefix = "ref"):
		"""Extract batches of k-mers from seq.
		
		Arguments:
			seq {string} -- input sequence
			k {int} -- substring length
			t {om.NATYPES} -- nucleic acid type
		
		Keyword Arguments:
			batchSize {number} -- number of kmers per batch (default: {1})
			prefix {str} -- reference record name (default: {"ref"})
		
		Returns:
			generator -- k-mer batch generator
		"""
		assert batchSize >= 1
		if batchSize == 1:
			return Sequence.kmerator(seq, k, t, prefix)
		else:
			for (seq2beKmered, i) in Sequence.batcher(seq, k, batchSize):
				yield Sequence.kmerator(seq2beKmered, k, t, prefix,
					prefix_start = i)

class KMer(Sequence):
	"""docstring for KMer"""

	class STRAND(Enum):
		"""docstring for STRAND"""
		PLUS = 1
		MINUS = 2

	def __init__(self, chrom, start, end, seq,
		t = om.NATYPES.DNA, strand = STRAND.PLUS):
		super().__init__(seq, t)
		assert start >= 0
		assert end >= 0
		assert strand in self.STRAND
		self.__chrom = chrom
		self.__start = start
		self.__end = end
		self.__strand = strand

	@property
	def header(self):
		return "%s:%d-%d" % (self.__chrom, self.__start, self.__end)
	@property
	def seq(self):
		return self.text

	@staticmethod
	def from_fasta(record, t = om.NATYPES.DNA):
		"""Reads a KMer from a Fasta record.
		
		Arguments:
			record {tuple} -- (header, seq)
		
		Keyword Arguments:
			t {om.NATYPES} -- nucleic acid type (default: {om.NATYPES.DNA})
		
		Returns:
			KMer
		"""
		chrom, coords = record[0].split(":")
		start, end = [int(n) for n in coords.split("-")]
		return KMer(chrom, start, end, record[1], t)
	
	@staticmethod
	def from_file(*args, **kwargs):
		return KMer.from_fasta(*args, **kwargs)

	def as_fasta(self):
		"""Fasta-like representation."""
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

class SequenceCount(Sequence):
	"""docstring for SequenceCount"""

	__headers = []

	def __init__(self, seq, headers, t = om.NATYPES.DNA):
		super().__init__(seq, t)
		self.__headers = headers

	@property
	def header(self):
		return self.__headers.copy()
	@property
	def seq(self):
		return self.text

	@staticmethod
	def from_text(line, t = om.NATYPES.DNA):
		"""Reads a KMer from a Fasta record.
		
		Arguments:
			record {tuple} -- (header, seq)
		
		Keyword Arguments:
			t {om.NATYPES} -- nucleic acid type (default: {om.NATYPES.DNA})
		
		Returns:
			KMer
		"""
		seq, headers = line.strip().split("\t")
		return SequenceCount(seq, headers.split(" "), t)

	@staticmethod
	def from_file(*args, **kwargs):
		return SequenceCount.from_text(*args, **kwargs)

	def as_text(self):
		return "%s\t%s\n" % (self.seq, " ".join(self.header))
