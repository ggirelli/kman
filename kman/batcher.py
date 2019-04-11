
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from enum import Enum
from ggc.args import check_threads
import gzip
import itertools
from kman.batch import Batch
from kman.seq import KMer, Sequence
from kman.io import SmartFastaParser
from joblib import Parallel, delayed
import oligo_melting as om
import os
import tempfile
from tqdm import tqdm

class BatcherBase(object):
	"""Basic batching system.
	
	Builds a collection of equally sized batches. In each Batch, a record is an
	instance of a class of the given size. Batch record type is consistent
	across the collection.
	
	Variables:
		DEFAULT_BATCH_SIZE {int} -- default batch size
		DEFAULT_BATCH_TYPE {type} -- default batched record type
		DEFAULT_NATYPE {om.NATYPES} -- default nucleic acid type
		_tmp {tempfile.TemporaryDirectory}
		_batches {list} -- list of Batch instances
		__size {int} -- batch size
		_type {type} -- batched record type
		__natype {om.NATYPES} -- nucleic acid type
	"""

	DEFAULT_BATCH_SIZE = int(1e6)
	DEFAULT_BATCH_TYPE = KMer
	DEFAULT_NATYPE = om.NATYPES.DNA

	_tmpH = None
	_tmp = None
	_batches = None
	__size = DEFAULT_BATCH_SIZE
	_type = DEFAULT_BATCH_TYPE
	__natype = DEFAULT_NATYPE

	def __init__(self, size = None, natype = None, tmp = None):
		"""Initializes BatcherBase.
		
		Keyword Arguments:
			size {int} -- batching size (default: {None})
			natype {om.NATYPES} -- nucleic acid type
			tmp {tempfile.TemporaryDirectory}
		"""
		super().__init__()
		if type(None) != type(size):
			assert size >= 1
			self.__size = int(size)
		if type(None) != type(natype):
			assert natype in om.NATYPES
			self.__natype = natype
		if tempfile.TemporaryDirectory == type(tmp):
			self._tmpH = tmp
			self._tmp = tmp.name
		elif type(None) != type(tmp):
			self._tmp = tmp
		if type(None) == type(self._batches):
			self._batches = [Batch.from_batcher(self)]

	@property
	def size(self):
		return self.__size
	@property
	def type(self):
		return self._type
	@property
	def natype(self):
		return self.__natype
	@property
	def collection(self):
		return self._batches
	@property
	def tmp(self):
		if type(None) == type(self._tmp):
			self._tmpH = tempfile.TemporaryDirectory(prefix = "kmanBatch")
			self._tmp = self._tmpH.name
		return self._tmp

	def new_batch(self):
		"""Add a new empty batch to the current collection."""
		if self.collection[-1].is_full():
			self.collection[-1].write()
			self._batches.append(Batch.from_batcher(self))

	def add_record(self, record):
		"""Add a record to the currentcollection.
		
		The record is added to the last empty batch in the collection. If no
		empty batches are left, a new empty batch is added.
		
		Arguments:
			record
		"""
		self.new_batch() # Add new batch if needed
		self.collection[-1].add(record)

	def write_all(self, f = "as_fasta", doSort = False):
		"""Write all batches to file.
		
		Keyword Arguments:
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
			doSort {bool} -- whether to sort when writing (default: {False})
		"""
		for bi in range(len(self.collection)):
			if 0 != self.collection[bi].current_size:
				self.collection[bi].write(f, doSort)

class BatcherThreading(BatcherBase):
	"""Parallelized batching system.
	
	Extends BatcherBase for parallelization.
	
	Extends:
		BatcherBase
	
	Variables:
		__threads {number} -- number of threads for parallelization
	"""

	class FEED_MODE(Enum):
		"""Feeding modes.
		
		Used with feed_collection() method.
		
		Extends:
			Enum
		
		Variables:
			REPLACE {number} -- replace current collection.
			FLOW {number} -- flow records into current collection.
			APPEND {number} -- append to current collection.
		"""
		REPLACE = 1
		FLOW = 2
		APPEND = 3

	__threads = 1

	def __init__(self, threads = 1, size = None, natype = None, tmp = None):
		"""Initialize BatcherThreading.
		
		Keyword Arguments:
			threads {number} -- number of threads for parallelization
			                    (default: {1})
			size {int} -- batching size (default: {None})
			natype {om.NATYPES} -- nucleic acid type
			tmp {tempfile.TemporaryDirectory}
		"""
		super().__init__(size, natype, tmp)
		self.threads = threads

	@property
	def threads(self):
		return self.__threads
	@threads.setter
	def threads(self, t):
		self.__threads = check_threads(t)

	def feed_collection(self, new_collection, mode = FEED_MODE.FLOW):
		"""Feed new batch collection to the current one.
		
		Different modes of feeding are available, see documentation of
		BatcherThreading.FEED_MODE for more details.
		
		Arguments:
			new_collection {list} -- list of Batches
		
		Keyword Arguments:
			mode {BatcherThreading.FEED_MODE} -- (default: {FEED_MODE.FLOW})
		"""
		assert all([b.type == self.type for b in new_collection])
		if mode == self.FEED_MODE.REPLACE:
			self._batches = new_collection
		elif mode == self.FEED_MODE.FLOW:
			for bi in tqdm(range(len(new_collection)), desc = "Flowing"):
				batch = new_collection.pop()
				for record in batch.record_gen():
					self.add_record(record)
				batch.reset()
		elif mode == self.FEED_MODE.APPEND:
			self._batches.extend(new_collection)

	@staticmethod
	def from_files(dirPath, threads, t = KMer, isFasta = True, reSort = False):
		"""Load batches from file.
		
		Each file in the provided directory should be a written Batch.
		
		Arguments:
			dirPath {str} -- path to batch directory
			threads {int} -- number of threads for parallelization
		
		Keyword Arguments:
			t {class} -- type of batch record (default: {KMer})
			isFasta {bool} -- whether batches are fasta files (default: {True})
		
		Returns:
			list -- list of Batches
		"""
		assert os.path.isdir(dirPath)
		threads = check_threads(threads)
		if 1 == threads:
			return [Batch.from_file(os.path.join(dirPath, fname), t, isFasta)
				for fname in tqdm(os.listdir(dirPath))]
		else:
			return Parallel(n_jobs = threads, verbose = 11)(
				delayed(Batch.from_file)(
					os.path.join(dirPath, fname), t, isFasta, reSort = reSort)
				for fname in os.listdir(dirPath))

class FastaBatcher(BatcherThreading):
	"""FASTA file k-mer batching.
	
	Divides k-mer from the records of a FASTA file into batches.
	
	Extends:
		BatcherThreading

	Variables:
		_doReverseComplement {bool} -- whether to batch also the reverse
		                               complement fo the sequences
	"""

	class MODE(Enum):
		KMERS = 1
		RECORDS = 2

	_doReverseComplement = False
	_mode = MODE.KMERS

	def __init__(self, threads = 1, size = None, natype = None, tmp = None):
		"""Initialize FastaBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			threads {int} -- number of threads for parallelization
			                 (default: {1})
			size {int} -- batch size (overridden by parent.size)
			natype {om.NATYPES} -- nucleic acid type
			tmp {tempfile.TemporaryDirectory}
		"""
		super().__init__(threads, size, natype, tmp)

	@property
	def mode(self):
		return self._mode
	@mode.setter
	def mode(self, m):
		assert m in self.MODE
		self._mode = m
	@property
	def doReverseComplement(self):
		return self._doReverseComplement
	@doReverseComplement.setter
	def doReverseComplement(self, rc):
		assert type(True) == type(rc)
		self._doReverseComplement = rc
	
	def __do_over_kmers(self, FH, k,
		feedMode = BatcherThreading.FEED_MODE.APPEND):
		"""Parallelize over kmers.
		
		Use RecordBatcher to parallelize when batching the kmers.
		
		Arguments:
			FH {io.TextIOWrapper} -- fasta file buffer
			k {int} -- k-mer length
			feedMode {BatcherThreading.FEED_MODE}
		"""
		batcher = FastaRecordBatcher.from_parent(self)
		for record in SmartFastaParser(FH).parse():
			batcher.do(record, k)
			for batch in batcher.collection:
				batch.unwrite()
		self.feed_collection(batcher.collection, feedMode)
		self.write_all(doSort = True)

	def __do_over_records(self, FH, k,
		feedMode = BatcherThreading.FEED_MODE.APPEND):
		"""Parallelize over FASTA records.
		
		Ran a non-parallelized RecordBatcher for each fasta record, in parallel.
		
		Arguments:
			FH {io.TextIOWrapper} -- fasta file buffer
			k {int} -- k-mer length
			feedMode {BatcherThreading.FEED_MODE}
		"""
		def do_record(size, natype, tmp, record, k):
			"""Batch a single record.
			
			Function to be passed to Parallel(delayed(*)).
			
			Arguments:
				fastaBatcher {FastaBatcher} -- batcher
				record {tuple} -- (header, sequence)
				k {int} -- k-mer length
			
			Returns:
				list -- list of Batches
			"""
			batcher = FastaRecordBatcher(1, size, natype, tmp)
			batcher.do(record, k, False)
			return batcher.collection

		batchCollections = Parallel(n_jobs = self.threads, verbose = 11
			)(delayed(do_record)(self.size, self.natype, self.tmp,
				record, k) for record in SmartFastaParser(FH).parse())

		self.feed_collection(list(itertools.chain(
			*batchCollections)), feedMode)
		self.write_all(doSort = True)

	def do(self, fasta, k, feedMode = BatcherThreading.FEED_MODE.APPEND):
		"""Start batching the fasta file.
		
		Batches a fasta file up to the specified number (self.size) of k-mers.
		
		Arguments:
			fasta {string} -- path to fasta file
			k {int} -- length of k-mers
			feedMode {BatcherThreading.FEED_MODE}
		"""
		assert os.path.isfile(fasta)
		assert k > 1

		if fasta.endswith(".gz"):
			FH = gzip.open(fasta, "rt") 
		else:
			FH = open(fasta, "r+")

		if self._mode == self.MODE.KMERS:
			self.__do_over_kmers(FH, k, feedMode)
		elif self._mode == self.MODE.RECORDS:
			self.__do_over_records(FH, k, feedMode)

		FH.close()

class FastaRecordBatcher(BatcherThreading):
	"""FASTA record batchin system.
	
	Divides k-mer from a single FASTA record into batches.
	
	Extends:
		BatcherThreading

	Variables:
		_doReverseComplement {bool} -- whether to batch also the reverse
		                               complement fo the sequences
	"""

	_doReverseComplement = False

	def __init__(self, threads = 1, size = None, natype = None, tmp = None):
		"""Initialize FastaRecordBatcher instance.
		
		Keyword Arguments:
			threads {int} -- number of threads for parallelization
			                 (default: {1})
			size {int} -- batch size (overridden by parent.size)
			natype {om.NATYPES} -- nucleic acid type
			tmp {tempfile.TemporaryDirectory}
		"""
		super().__init__(threads, size, natype, tmp)

	@property
	def doReverseComplement(self):
		return self._doReverseComplement
	@doReverseComplement.setter
	def doReverseComplement(self, rc):
		assert type(True) == type(rc)
		self._doReverseComplement = rc

	def do(self, record, k, verbose = True):
		"""Start batching a fasta record.
		
		Requires a fasta record with header and sequence.
		
		Arguments:
			record {tuple} -- (header, sequence)
			k {int} -- length of k-mers
		"""
		record_name = record[0].split(" ")[0]
		if verbose: print("Batching record '%s'..." % record_name)
		if 1 == self.threads:
			kmerGen = Sequence.kmerator(record[1], k, self.natype, record_name,
				rc = self.doReverseComplement)
			if verbose: kmerGen = tqdm(kmerGen)
			for kmer in kmerGen:
				if kmer.is_ab_checked():
					self.add_record(kmer)
		else:
			batches = Parallel(n_jobs = self.threads, verbose = 11
				)(delayed(FastaRecordBatcher.build_batch
					)(seq, record_name, k, self, i)
					for (seq, i) in Sequence.batcher(record[1], k, self.size))
			self.feed_collection(batches, self.FEED_MODE.APPEND)
		self.write_all()

	@staticmethod
	def build_batch(seq, name, k, batcher, i = 0):
		"""Builds a Batch.
		
		Arguments:
			seq {string} -- sequence
			name {string} -- header
			k {int} -- k for k-mering
			batcher {BatcherBase} -- parent
		
		Keyword Arguments:
			i {number} -- position offset (default: {0})
		
		Returns:
			Batch
		"""
		batch = Batch.from_batcher(batcher)
		recordGen = Sequence.kmerator(seq, k, batcher.natype, name, i,
			rc = batcher.doReverseComplement)
		batch.add_all((k for k in recordGen if k.is_ab_checked()))
		batch.write(doSort = True)
		return batch

	@staticmethod
	def from_parent(parent):
		"""Initialize FastaRecordBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			parent {Batcher} -- parent batcher to inherit attributes from
			                    (default: {None})
		"""
		batcher = FastaRecordBatcher(parent.threads, parent.size,
			parent.natype, parent.tmp)
		batcher._doReverseComplement = parent.doReverseComplement
		return batcher
