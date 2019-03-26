
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
import io
from kman.seq import KMer, Sequence
import oligo_melting as om
import os
import tempfile
from time import time
from tqdm import tqdm

class Batcher(object):
	"""docstring for Batcher"""

	DEFAULT_BATCH_SIZE = 1e6
	DEFAULT_BATCH_JOIN_SIZE = 10
	DEFAULT_BATCH_TYPE = KMer
	DEFAULT_NATYPE = om.NATYPES.DNA

	_tmp = None
	_batches = None
	__batch_join_size = DEFAULT_BATCH_JOIN_SIZE
	__size = DEFAULT_BATCH_SIZE
	__type = DEFAULT_BATCH_TYPE
	__natype = DEFAULT_NATYPE

	def __init__(self, size = None):
		super(Batcher, self).__init__()
		if type(None) == type(self.tmp):
			self._tmp = tempfile.TemporaryDirectory(prefix = "kmanBatch")
		if type(None) != type(size):
			assert size >= 1
			self.__size = size
		if type(None) == type(self._batches):
			self._batches = [Batch(self)]

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
		return self._batches
	@property
	def tmp(self):
		return self._tmp

	def new_batch(self):
		if self.batches[-1].is_full():
			self.batches[-1].write()
			self._batches.append(Batch(self))

	def add_record(self, record):
		self.new_batch() # Add new batch if needed
		self.batches[-1].add(record)

class FastaBatcher(Batcher):
	"""docstring for FastaBatcher"""

	def __init__(self, size = None):
		super(FastaBatcher, self).__init__(size)

	def do(self, fasta, k):
		"""Start batching the fasta file.
		
		Batches a fasta file up to the specified number (self.size) of k-mers.
		
		Arguments:
			fasta {string} -- path to fasta file.
			k {int} -- length of k-mers
		"""
		assert os.path.isfile(fasta)
		assert k > 1

		batcher = RecordBatcher(parent = self)
		with open(fasta, "r+") as FH:
			for record in SimpleFastaParser(FH):
				batcher.do(record, k)

class RecordBatcher(Batcher):
	"""docstring for RecordBatcher"""

	def __init__(self, size = None, parent = None):
		"""Initialize RecordBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			parent {Batcher} -- [description] (default: {None})
		"""
		if type(None) != type(parent):
			self.__size = parent.size
			size = parent.size
			self.__natype = parent.natype
			self._batches = parent.batches
			self._tmp = parent.tmp
		super(RecordBatcher, self).__init__(size)

	def do(self, record, k):
		record_name = record[0].split(" ")[0]
		for kmer in tqdm(Sequence.kmerator(
			record[1], k, self.natype, record_name)):
			if kmer.is_ab_checked():
				self.add_record(kmer)

class Batch(object):
	"""docstring for Batch"""
	
	__written = False

	def __init__(self, batcher):
		super(Batch, self).__init__()
		assert batcher.size >= 1
		self.__size = int(batcher.size)
		self.__records = [None] * self.__size
		self.__i = 0
		self.__type = batcher.type
		self.__tmp = tempfile.NamedTemporaryFile(mode = "w+",
			dir = batcher.tmp.name, prefix = str(hash(time())), suffix = ".fa")

	@property
	def is_written(self):
		return self.__written
	@property
	def current_size(self):
		return self.__i
	@property
	def size(self):
		return self.__size
	@property
	def type(self):
		return self.__type
	@property
	def tmp(self):
		return self.__tmp.name
	
	@property
	def sorted(self):
		if self.is_written:
			return sorted(self.record_gen(), key = lambda x: x[1])
		else:
			return sorted(self.record_gen(), key = lambda x: x.seq)

	def record_gen(self):
		if self.is_written:
			self.__tmp.seek(0)
			return (r for r in SimpleFastaParser(self.__tmp))
		else:
			return (r for r in self.__records)

	def add(self, record):
		"""Add a record to the current batch.
		
		Does not work if the batch is full. Also, the record type must match the
		batch type.
		
		Arguments:
			record -- record of the same type as self.type
		"""
		assert not self.is_written, "this batch has been stored locally."
		assert not self.is_full(), "this batch is full."
		assert type(record) == self.type, "record must be %s, not %s" % (
			self.type, type(record))
		self.__records[self.__i] = record
		self.__i += 1

	def to_write(self, f = "as_fasta"):
		"""Generator of writeable records.
		
		Generator function that converts batch records into writeable ones.
		
		Keyword Arguments:
			f {str} -- name of record in-built function for writeable format]
			           (default: {"as_fasta"})
		"""
		return (getattr(r, f)() for r in self.__records
			if not type(None) == type(r))

	def write(self, f = "as_fasta"):
		self.__tmp.write("".join(self.to_write(f)))
		self.__records = None
		self.__written = True

	def is_full(self):
		return self.current_size == self.__size
