
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
from kman.seq import KMer, Sequence
import oligo_melting as om
import os
from tqdm import tqdm

class Batcher(object):
	"""docstring for Batcher"""

	DEFAULT_BATCH_SIZE = 1e6
	DEFAULT_BATCH_TYPE = KMer
	DEFAULT_NATYPE = om.NATYPES.DNA

	__size = DEFAULT_BATCH_SIZE
	__type = DEFAULT_BATCH_TYPE
	__natype = DEFAULT_NATYPE

	def __init__(self):
		super(Batcher, self).__init__()
		self.__batches = [Batch(self.size, self.type)]

	@property
	def size(self):
		return self.__size

	@property
	def type(self):
		return self.__type

	@property
	def natype(self):
		return self.__natype

	@property
	def batches(self):
		return self.__batches

	def new_batch(self):
		if self.batches[-1].is_full():
			self.__batches.append(Batch(self.size, self.type))

	def add_record(self, record):
		self.new_batch() # Add new batch if needed
		self.batches[-1].add(record)

class FastaBatcher(Batcher):
	"""docstring for FastaBatcher"""

	def __init__(self):
		super(FastaBatcher, self).__init__()

	def do(self, fasta, k):
		"""Start batching the fasta file.
		
		Batches a fasta file up to the specified number (self.size) of k-mers.
		
		Arguments:
			fasta {string} -- path to fasta file.
			k {int} -- length of k-mers
		"""
		assert os.path.isfile(fasta)
		assert k > 1

		batcher = RecordBatcher(self)
		with open(fasta, "r+") as FH:
			for record in SimpleFastaParser(FH):
				batcher.do(record, k)

class RecordBatcher(Batcher):
	"""docstring for RecordBatcher"""

	def __init__(self, parent = None):
		"""Initialize RecordBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			parent {Batcher} -- [description] (default: {None})
		"""
		super(RecordBatcher, self).__init__()
		if type(None) != type(parent):
			self.__size = parent.size
			self.__natype = parent.natype
			self.__batches = parent.batches

	def do(self, record, k):
		record_name = record[0].split(" ")[0]
		for kmer in tqdm(Sequence.kmerator(record[1], k, self.natype, record_name)):
			if kmer.is_ab_checked():
				self.add_record(kmer)

class Batch(object):
	"""docstring for Batch"""
	def __init__(self, size, t):
		super(Batch, self).__init__()
		assert size >= 1
		self.__size = int(size)
		self.__records = [None] * self.__size
		self.__i = 0
		self.__type = t

	@property
	def current_size(self):
		return self.__i
	@property
	def size(self):
		return self.__size
	@property
	def sorted(self):
		return sorted(self.__records, key = lambda x: x.seq)
	@property
	def type(self):
		return self.__type

	def add(self, record):
		assert not self.is_full(), "this batch is full."
		assert type(record) == self.type, "record must be %s, not %s" % (
			self.type, type(record))
		self.__records[self.__i] = record
		self.__i += 1

	def fastaGen(self):
		return (">%s\n%s\n" % (r.header, r.seq)
			for r in self.__records if not type(None) == type(r))

	def fasta(self):
		return("".join(list(self.fastaGen())))

	def is_full(self):
		return self.current_size == self.__size
