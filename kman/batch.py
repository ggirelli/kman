
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
from enum import Enum
from ggc.args import check_threads
import gzip
from kman.seq import KMer, Sequence
from joblib import Parallel, delayed
import oligo_melting as om
import os
import tempfile
import time
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

	_tmp = None
	_batches = None
	__size = DEFAULT_BATCH_SIZE
	_type = DEFAULT_BATCH_TYPE
	__natype = DEFAULT_NATYPE

	def __init__(self, size = None):
		"""Initializes BatcherBase.
		
		Keyword Arguments:
			size {int} -- batching size (default: {None})
		"""
		super().__init__()
		if type(None) != type(size):
			assert size >= 1
			self.__size = int(size)
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
			self._tmp = tempfile.TemporaryDirectory(prefix = "kmanBatch")
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

	def __init__(self, threads = 1, size = None):
		"""Initialize BatcherThreading.
		
		Keyword Arguments:
			threads {number} -- number of threads for parallelization
			                    (default: {1})
			size {int} -- batching size (default: {None})
		"""
		super().__init__(size)
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
			for bi in tqdm(range(len(new_collection))):
				batch = new_collection.pop()
				for record in batch.record_gen():
					self.add_record(record)
				batch.reset()
		elif mode == self.FEED_MODE.APPEND:
			self._batches.extend(new_collection)

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

	def __init__(self, threads = 1, size = None):
		"""Initialize FastaBatcher instance.
		
		A parent batcher class can be specified, whose attributes are inherited.
		
		Keyword Arguments:
			threads {int} -- number of threads for parallelization
			                 (overridden by parent.threads)
			size {int} -- batch size (overridden by parent.size)
		"""
		super().__init__(threads, size)

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
	
	def __do_over_kmers(self, FH, k):
		"""Parallelize over kmers.
		
		Use RecordBatcher to parallelize when batching the kmers.
		
		Arguments:
			FH {io.TextIOWrapper} -- fasta file buffer
			k {int} -- k-mer length
		"""
		batcher = FastaRecordBatcher(parent = self)
		for record in SimpleFastaParser(FH):
			batcher.do(record, k)
			if 1 != self.threads:
				self.feed_collection(batcher.collection,
					self.FEED_MODE.APPEND)

	def __do_over_records(self, FH, k):
		"""Parallelize over FASTA records.
		
		Ran a non-parallelized RecordBatcher for each fasta record, in parallel.
		
		Arguments:
			FH {io.TextIOWrapper} -- fasta file buffer
			k {int} -- k-mer length
		"""
		def do_record(fastaBatcher, record, k):
			"""Batch a single record.
			
			Function to be passed to Parallel(delayed(*)).
			
			Arguments:
				fastaBatcher {FastaBatcher} -- batcher
				record {tuple} -- (header, sequence)
				k {int} -- k-mer length
			
			Returns:
				list -- list of Batches
			"""
			batcher = FastaRecordBatcher(parent = self)
			batcher.threads = 1
			batcher.do(record, k, False)
			return batcher.collection

		batchCollections = Parallel(n_jobs = self.threads, verbose = 11
			)(delayed(do_record)(self, record, k)
			for record in SimpleFastaParser(FH))
		for collection in batchCollections:
			self.feed_collection(collection, self.FEED_MODE.APPEND)


	def do(self, fasta, k):
		"""Start batching the fasta file.
		
		Batches a fasta file up to the specified number (self.size) of k-mers.
		
		Arguments:
			fasta {string} -- path to fasta file
			k {int} -- length of k-mers
		"""
		assert os.path.isfile(fasta)
		assert k > 1

		if fasta.endswith(".gz"):
			FH = gzip.open(fasta, "rt") 
		else:
			FH = open(fasta, "r+")

		if self._mode == self.MODE.KMERS:
			self.__do_over_kmers(FH, k)
		elif self._mode == self.MODE.RECORDS:
			self.__do_over_records(FH, k)

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

	def __init__(self, threads = 1, size = None, parent = None):
		"""Initialize FastaRecordBatcher instance.
		
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
			self._doReverseComplement = parent.doReverseComplement
			self.__natype = parent.natype
			self._batches = parent.collection
			self._tmp = parent.tmp
		super().__init__(threads, size)

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
			self.feed_collection(batches, self.FEED_MODE.REPLACE)

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

class Batch(object):
	"""Batch container.
	
	Records in the Batch are accessible through the record_gen and sorted
	methods, which source either from memory or from written files. A Batch
	cannot be resized. After full size is reached, a new Batch should be created
	
	Variables:
		__written {bool} -- if the batch was written to file
		__i {number} -- current record location
		__tmp_dir {tempfile.TemporaryDirectory}
		__tmp {tempfile.TemporaryFile}
		isFasta {bool} -- whether the output should be in fasta format
		suffix {str} -- extension for the output temporary file
	"""
	
	__written = False
	__i = 0
	__tmp_dir = None
	__tmp = None
	isFasta = True
	suffix = ".fa"

	def __init__(self, t, tmpDir, size = 1):
		"""Initialize a Batch.
		
		Arguments:
			t {class} -- record type
			tmpDir {str} -- path to temporary directory
		
		Keyword Arguments:
			size {number} -- batch size (default: {1})
		"""
		super().__init__()
		assert size >= 1
		self.__size = int(size)
		self.__remaining = self.__size
		self.__records = [None] * self.__size
		self.__type = t
		self.__tmp_dir = tmpDir

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
		if type(None) == type(self.__tmp):
			with tempfile.NamedTemporaryFile(mode = "w+",
				dir = self.__tmp_dir,
				prefix = str(hash(time.time())),
				suffix = self.suffix) as TH:
				self.__tmp = TH.name
		return self.__tmp
	@property
	def info(self):
		"""Prints Batch information in a readable format."""
		info = "%s\ntype: %s\nsize: %d" % (self.tmp, self.type, self.size)
		info += "\ni: %d\nremaining: %d" %(self.current_size, self.remaining)
		info += "\nwritten: %r\n" % self.is_written
		return info

	@property
	def sorted(self, keyAttr = "seq"):
		"""Generator of sorted records.
		
		Keyword Arguments:
			keyAttr {str} -- attribute key for sorting (default: {"seq"})
		
		Returns:
			generator
		"""
		return sorted(self.record_gen(), key = lambda x: getattr(x, keyAttr))

	def record_gen(self):
		"""Generator of records.
		
		Yields:
			record
		"""
		if self.is_written:
			if self.tmp.endswith(".gz"):
				TH = gzip.open(self.tmp, "rt")
			else:
				TH = open(self.tmp, "r+")
			if self.isFasta:
				for record in SimpleFastaParser(TH):
					yield self.__type.from_file(record)
			else:
				for line in TH:
					yield self.__type.from_file(line)
			TH.close()
		else:
			for record in self.__records:
				if not type(None) == type(record):
					yield record

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
		"""Adds all records from a generator to the current Batch.
		
		Arguments:
			recordGen {generator} -- record generator
		"""
		for record in recordGen:
			self.add(record)

	def to_write(self, f = "as_fasta", doSort = False):
		"""Generator of writeable records.
		
		Generator function that converts batch records into writeable ones.
		
		Keyword Arguments:
			f {str} -- name of record in-built function for writeable format]
			           (default: {"as_fasta"})
		"""
		if doSort:
			return (getattr(r, f)() for r in self.sorted
				if not type(None) == type(r))
		else:
			return (getattr(r, f)() for r in self.record_gen()
				if not type(None) == type(r))

	def write(self, f = "as_fasta", doSort = False):
		"""Writes the batch to file.
		
		Keyword Arguments:
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
			doSort {bool} -- whether to sort when writing (default: {False})
		"""
		with open(self.tmp, "w+") as TH:
			TH.write("".join(self.to_write(f, doSort)))
		self.__records = None
		self.__written = True

	@staticmethod
	def from_file(path, t = KMer, isFasta = True):
		"""Generate a Batch from a file.
		
		Used to link an existing file to a Batch.
		
		Arguments:
			path {str} -- path to previously generated batch
		
		Keyword Arguments:
			t {class} -- record type (default: {KMer})
			isFasta {bool} -- whether the input is a fasta (default: {True})
		"""
		if isFasta and path.endswith(".gz"):
			with gzip.open(path, "rt") as FH:
				size = sum([1 for record in SimpleFastaParser(FH)])
		else:
			with open(path, "r+") as FH:
				size = sum([1 for line in FH])
		batch = Batch(t, os.path.dirname(path), size)
		batch.__tmp = FH.name
		batch.__i = size
		batch.__remaining = 0
		batch.__written = True
		return(batch)

	@staticmethod
	def from_batcher(batcher, size = 1):
		"""Initialize a Batch.
		
		Uses a Batcher to set default values.
		
		Arguments:
			batcher {BatcherBase} -- parent batcher.
		"""
		assert batcher.size >= 1
		size = max(int(batcher.size), size)
		return Batch(batcher.type, batcher.tmp.name, size)

	def reset(self):
		"""Reset the batch.
		
		Empties current records collection and any written file.
		"""
		os.remove(self.tmp)
		self.__written = False
		self.__i = 0
		self.__remaining = self.size
		self.__records = [None] * self.__size

	def is_full(self):
		"""Whether the Batch collection is full."""
		return 0 == self.remaining
