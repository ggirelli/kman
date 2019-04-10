
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
from kman.seq import KMer
from kman.io import SmartFastaParser
import os
import tempfile
import time

class Batch(object):
	"""Batch container.
	
	Records in the Batch are accessible through the record_gen and sorted
	methods, which source either from memory or from written files. A Batch
	cannot be resized. After full size is reached, a new Batch should be created
	
	Variables:
		_written {bool} -- if the batch was written to file
		_i {number} -- current record location
		_tmp_dir {tempfile.TemporaryDirectory}
		_tmp {tempfile.TemporaryFile}
		isFasta {bool} -- whether the output should be in fasta format
		suffix {str} -- extension for the output temporary file
	"""
	
	_written = False
	_i = 0
	_tmp_dir = None
	_tmp = None
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
		self._remaining = self.__size
		self.__type = t
		self._tmp_dir = tmpDir
		self.__records = [None] * self.__size

	@property
	def is_written(self):
		return self._written
	@property
	def current_size(self):
		return self._i
	@property
	def size(self):
		return self.__size
	@property
	def remaining(self):
		return self._remaining
	@property
	def type(self):
		return self.__type
	@property
	def tmp(self):
		if type(None) == type(self._tmp):
			with tempfile.NamedTemporaryFile(mode = "w+",
				dir = self._tmp_dir,
				prefix = str(hash(time.time())),
				suffix = self.suffix) as TH:
				self._tmp = TH.name
		return self._tmp
	@property
	def info(self):
		"""Prints Batch information in a readable format."""
		info = "%s\ntype: %s\nsize: %d" % (self.tmp, self.type, self.size)
		info += "\ni: %d\nremaining: %d" %(self.current_size, self.remaining)
		info += "\nwritten: %r\n" % self.is_written
		return info

	@property
	def sorted(self, smart = False, keyAttr = "seq"):
		"""Generator of sorted records.
		
		Keyword Arguments:
			keyAttr {str} -- attribute key for sorting (default: {"seq"})
		
		Returns:
			generator
		"""
		return sorted(self.record_gen(smart),
			key = lambda x: getattr(x, keyAttr))

	def _record_gen_from_file(self, smart = False):
		"""Generator of records, reading from file.
		
		Yields:
			record
		"""
		if self.tmp.endswith(".gz"):
			TH = gzip.open(self.tmp, "rt")
		else:
			TH = open(self.tmp, "r+")
		if self.isFasta:
			if smart:
				for record in SmartFastaParser(TH).parse():
					yield self.__type.from_file(record)
			else:
				for record in SimpleFastaParser(TH):
					yield self.__type.from_file(record)
		else:
			for line in TH:
				yield self.__type.from_file(line)
		if not TH.closed:
			TH.close()

	def record_gen(self, smart = False):
		"""Generator of records.
		
		Yields:
			record
		"""
		if self.is_written:
			for record in self._record_gen_from_file(smart):
				yield record
		else:
			for record in self.__records:
				if not type(None) == type(record):
					yield record

	def check_record(self, record):
		"""Check that record type matches the Batch."""
		assert_msg = "record must be %s, not %s." % (self.type, type(record))
		assert type(record) == self.type, assert_msg

	def add(self, record):
		"""Add a record to the current batch.
		
		Does not work if the batch is full. Also, the record type must match the
		batch type. Automatically selects whether to append in memory or to the
		temporary file (if self.is_appending).
		
		Arguments:
			record -- record of the same type as self.type
		"""
		assert not self.is_full(), "this batch is full."
		assert not self.is_written, "this batch has been stored locally."
		self.check_record(record)
		self.__records[self._i] = record
		self._i += 1
		self._remaining -= 1

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

	def write(self, f = "as_fasta", doSort = False, force = False):
		"""Writes the batch to file.
		
		Keyword Arguments:
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
			doSort {bool} -- whether to sort when writing (default: {False})
			force {bool} -- force overwriting, useful with doSort.
			                (default: {False})
		"""
		if not self.is_written or force:
			output = "".join(self.to_write(f, doSort))
			with open(self.tmp, "w+") as TH:
				TH.write(output)
			self.__records = None
			self._written = True

	@staticmethod
	def from_file(path, t = KMer, isFasta = True,
		smart = False, reSort = False):
		"""Generate a Batch from a file.
		
		Used to link an existing file to a Batch.
		
		Arguments:
			path {str} -- path to previously generated batch
		
		Keyword Arguments:
			t {class} -- record type (default: {KMer})
			isFasta {bool} -- whether the input is a fasta (default: {False})
			smart {bool} -- use smarter IO (open only when needed) parser when
			                available. Might cause higher overhead.
			reSort {bool} -- sort the written batch. (default: {False})
		"""

		if isFasta:
			if path.endswith(".gz"):
				FH = gzip.open(path, "rt")
			else:
				FH = open(path, "r+")
			
			if smart:
				size = max(2, sum(1 for record in SmartFastaParser(FH).parse()))
			else:
				size = max(2, sum(1 for record in SimpleFastaParser(FH)))
		else:
			FH = open(path, "r+")
			size = max(2, sum(1 for line in FH))

		batch = Batch(t, os.path.dirname(path), size)
		batch._tmp = FH.name
		batch._i = size
		batch._remaining = 0
		batch._written = True

		if reSort:
			batch.write(doSort = True, force = True)

		FH.close()
		return batch

	@staticmethod
	def from_batcher(batcher, size = 1):
		"""Initialize a Batch.
		
		Uses a Batcher to set default values.
		
		Arguments:
			batcher {BatcherBase} -- parent batcher.
			size {number} -- batch size. Uses the largest between provided here
			                 and batcher.size. (default: {1})
		"""
		assert batcher.size >= 1
		size = max(int(batcher.size), size)
		return Batch(batcher.type, batcher.tmp, size)

	def reset(self):
		"""Reset the batch.
		
		Empties current records collection and any written file.
		"""
		os.remove(self.tmp)
		self._written = False
		self._i = 0
		self._remaining = self.size
		self.__records = [None] * self.__size

	def is_full(self):
		"""Whether the Batch collection is full."""
		return 0 == self.remaining

class BatchAppendable(Batch):
	"""Batch container.
	
	Records in the Batch are accessible through the record_gen and sorted
	methods, which source either from memory or from written files. A Batch
	cannot be resized. After full size is reached, a new Batch should be created
	
	Variables:
		_written {bool} -- if the batch was written to file
		_i {number} -- current record location
		_tmp_dir {tempfile.TemporaryDirectory}
		_tmp {tempfile.TemporaryFile}
		isFasta {bool} -- whether the output should be in fasta format
		suffix {str} -- extension for the output temporary file
	"""
	
	def __init__(self, t, tmpDir, size = 1):
		"""Initialize a Batch.
		
		Arguments:
			t {class} -- record type
			tmpDir {str} -- path to temporary directory
		
		Keyword Arguments:
			size {number} -- batch size (default: {1})
		"""
		super().__init__(t, tmpDir, size)
		self._written = True

	@property
	def info(self):
		"""Prints Batch information in a readable format."""
		info = "%s\ntype: %s\nsize: %d" % (self.tmp, self.type, self.size)
		info += "\ni: %d\nremaining: %d" %(self.current_size, self.remaining)
		info += "\nappending: True\n"
		return info

	def record_gen(self, smart = False):
		"""Generator of records.
		
		Yields:
			record
		"""
		for record in self._record_gen_from_file(smart):
			yield record

	def add(self, record, f = "as_fasta"):
		"""Add a record to the current batch.
		
		Does not work if the batch is full. Also, the record type must match the
		batch type. Automatically selects whether to append in memory or to the
		temporary file (if self.is_appending).
		
		Arguments:
			record -- record of the same type as self.type
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
		"""
		assert not self.is_full(), "this batch is full."
		super().check_record(record)
		with open(self.tmp, "a+") as OH:
			OH.write(getattr(record, f)())
		self._i += 1
		self._remaining -= 1

	def add_all(self, recordGen, f = "as_fasta"):
		"""Adds all records from a generator to the current Batch.
		
		Arguments:
			recordGen {generator} -- record generator
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
		"""
		for record in recordGen:
			self.add(record, f)

	def write(self, f = "as_fasta", doSort = False):
		"""Writes the batch to file.

		Does something only if doSort is True.
		
		Keyword Arguments:
			f {str} -- name of method in records class for string-like
			           representation (default: {"as_fasta"})
			doSort {bool} -- whether to sort when writing (default: {False})
		"""
		if doSort:
			output = "".join(self.to_write(f, doSort))
			with open(self.tmp, "w+") as TH:
				TH.write(output)

	@staticmethod
	def from_file(path, t = KMer, isFasta = True, smart = False):
		"""Generate a Batch from a file.
		
		Used to link an existing file to a Batch.
		
		Arguments:
			path {str} -- path to previously generated batch
		
		Keyword Arguments:
			t {class} -- record type (default: {KMer})
			isFasta {bool} -- whether the input is a fasta (default: {False})
			smart {bool} -- use smarter IO (open only when needed) parser when
			                available. Might cause higher overhead.
		"""

		if isFasta:
			if path.endswith(".gz"):
				FH = gzip.open(path, "rt")
			else:
				FH = open(path, "r+")
			
			if smart:
				size = max(2, sum(1 for record in SmartFastaParser(FH).parse()))
			else:
				size = max(2, sum(1 for record in SimpleFastaParser(FH)))
		else:
			FH = open(path, "r+")
			size = max(2, sum(1 for line in FH))

		batch = BatchAppendable(t, os.path.dirname(path), size)
		batch._tmp = FH.name
		batch._i = size
		batch._remaining = 0
		batch._written = True

		FH.close()
		return batch

	@staticmethod
	def from_batcher(batcher, size = 1):
		"""Initialize a Batch.
		
		Uses a Batcher to set default values.
		
		Arguments:
			batcher {BatcherBase} -- parent batcher.
			size {number} -- batch size. Uses the largest between provided here
			                 and batcher.size. (default: {1})
		"""
		assert batcher.size >= 1
		size = max(int(batcher.size), size)
		return BatchAppendable(batcher.type, batcher.tmp, size)

	def reset(self):
		"""Reset the batch.
		
		Empties current records collection and any written file.
		"""
		os.remove(self.tmp)
		self._i = 0
		self._remaining = self.size
