
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ggc.args import check_threads
from kman.seq import KMer, Sequence
from joblib import Parallel, delayed
import oligo_melting as om
import os
import tempfile
import time
from tqdm import tqdm

class Batcher(object):
	"""docstring for Batcher"""

	DEFAULT_BATCH_SIZE = int(1e6)
	DEFAULT_BATCH_TYPE = KMer
	DEFAULT_NATYPE = om.NATYPES.DNA

	_tmp = None
	_batches = None
	__size = DEFAULT_BATCH_SIZE
	__type = DEFAULT_BATCH_TYPE
	__natype = DEFAULT_NATYPE

	def __init__(self, size = None):
		super(Batcher, self).__init__()
		if type(None) != type(size):
			assert size >= 1
			self.__size = int(size)
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
		if type(None) == type(self._tmp):
			self._tmp = tempfile.TemporaryDirectory(prefix = "kmanBatch")
		return self._tmp

	def new_batch(self):
		if self.batches[-1].is_full():
			self.batches[-1].write()
			self._batches.append(Batch(self))

	def add_record(self, record):
		self.new_batch() # Add new batch if needed
		self.batches[-1].add(record)

class BatcherThreading(Batcher):
	"""docstring for BatcherThreading"""

	__threads = 1

	def __init__(self, threads = 1, size = None):
		super(BatcherThreading, self).__init__(size)
		self.threads = threads

	@property
	def threads(self):
		return self.__threads
	@threads.setter
	def threads(self, t):
		self.__threads = check_threads(t)

	def feed_batches(self, batchList, replace = False):
		assert all([type(x) == self.type for x in batchList])
		if replace:
			self._batches = batchList
		else:
			for batch in tqdm(batchList):
				for record in batch.record_gen():
					self.add_record(record)

class FastaBatcher(BatcherThreading):
	"""docstring for FastaBatcher"""

	def __init__(self, threads = 1, size = None):
		"""Initialize FastaBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			threads {int} -- number of threads for parallelization
			                 (overridden by parent.threads)
			size {int} -- batch size (overridden by parent.size)
		"""
		super(FastaBatcher, self).__init__(threads, size)

	def do(self, fasta, k):
		"""Start batching the fasta file.
		
		Batches a fasta file up to the specified number (self.size) of k-mers.
		
		Arguments:
			fasta {string} -- path to fasta file
			k {int} -- length of k-mers
		"""
		assert os.path.isfile(fasta)
		assert k > 1

		batcher = RecordBatcher(parent = self)
		with open(fasta, "r+") as FH:
			for record in SimpleFastaParser(FH):
				batcher.do(record, k)

class RecordBatcher(BatcherThreading):
	"""docstring for RecordBatcher"""


	def __init__(self, threads = 1, size = None, parent = None):
		"""Initialize RecordBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			threads {int} -- number of threads for parallelization
			                 (overridden by parent.threads)
			size {int} -- batch size (overridden by parent.size)
			parent {Batcher} -- parent batcher to inherit attributes from
			                    (default: {None})
		"""
		if type(None) != type(parent):
			self.__size = parent.size
			size = parent.size
			self.threads = parent.threads
			threads = parent.threads
			self.__natype = parent.natype
			self._batches = parent.batches
			self._tmp = parent.tmp
		super(RecordBatcher, self).__init__(threads, size)

	def do(self, record, k):
		"""Start batching a fasta record.
		
		Requires a fasta record with header and sequence.
		
		Arguments:
			record {tuple} -- (header, sequence)
			k {int} -- length of k-mers
		"""
		record_name = record[0].split(" ")[0]
		if 1 == self.threads:
			kmerGen = Sequence.kmerator(record[1], k, self.natype, record_name)
			for kmer in tqdm(kmerGen):
				if kmer.is_ab_checked():
					self.add_record(kmer)
		else:
			#assert self.threads == 1, "parallelization not implemented yet."
			seqBatchGen = Sequence.batcher(record[1], k, self.size)

			batchList = Parallel(n_jobs = self.threads, verbose = 11
				)(delayed(build_batch_for_parallel
					)(seq, record_name, k, self, i) for (seq, i) in seqBatchGen)
			
			self.feed_batches(batchList)

def build_batch_for_parallel(seq, name, k, batcher, i = 0):
	batch = Batch(batcher)
	recordGen = Sequence.kmerator(seq, k, batcher.natype, name, i)
	batch.add_all((k for k in recordGen if k.is_ab_checked()))
	batch.write()
	return batch

class Batch(object):
	"""docstring for Batch"""
	
	__written = False
	__i = 0

	def __init__(self, batcher):
		super(Batch, self).__init__()
		assert batcher.size >= 1
		self.__size = int(batcher.size)
		self.__remaining = self.__size
		self.__records = [None] * self.__size
		self.__type = batcher.type
		with tempfile.NamedTemporaryFile(mode = "w+",
			dir = batcher.tmp.name,
			prefix = str(hash(time.time())),
			suffix = ".fa") as TH:
			self.__tmp = TH.name

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
	def remaining(self):
		return self.__remaining
	@property
	def type(self):
		return self.__type
	@property
	def tmp(self):
		return self.__tmp
	@property
	def info(self):
		return "%s\ntype: %s\nsize: %d\ni: %d\nremaining: %d" %(self.tmp,
			self.type, self.size, self.current_size, self.remaining)

	@property
	def sorted(self):
		if self.is_written:
			return sorted(self.record_gen(), key = lambda x: x[1])
		else:
			return sorted(self.record_gen(), key = lambda x: x.seq)

	def record_gen(self):
		if self.is_written:
			with open(self.tmp, "r+") as TH:
				TH.seek(0)
				return (r for r in SimpleFastaParser(TH))
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
		assert type(record) == self.type, "record must be %s, not %s." % (
			self.type, type(record))
		self.__records[self.__i] = record
		self.__i += 1
		self.__remaining -= 1

	def add_all(self, recordGen):
		for record in recordGen:
			self.add(record)

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
		with open(self.tmp, "w+") as TH:
			TH.write("".join(self.to_write(f)))
		self.__records = None
		self.__written = True

	def is_full(self):
		return 0 == self.remaining
