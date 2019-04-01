
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
'''

from enum import Enum
from ggc.args import check_threads
import gzip
from heapq import merge
import io
from joblib import Parallel, delayed
from kman.batch import Batch, BatcherThreading
from kman.seq import SequenceCount
import numpy as np
import os
import re
import tempfile
import time
from tqdm import tqdm


class Crawler(object):
	"""docstring for Crawler"""

	doSort = False
	verbose = True

	def __init__(self):
		super().__init__()

	def count_records(self, batches):
		return sum([b.current_size for b in batches])

	def do_records(self, batches):
		assert all([type(b) == Batch for b in batches]), batches

		if self.doSort:
			generators = [((r.header, r.seq) for r in b.sorted)
				for b in batches]
		else:
			generators = [((r.header, r.seq) for r in b.record_gen())
				for b in batches]

		crawler = merge(*generators, key = lambda x: x[1])

		return crawler

	def do_batch(self, batches):
		crawler = self.do_records(batches)
		
		first_record = next(crawler)
		current_seq = first_record[1]
		current_headers = [first_record[0]]

		if self.verbose:
			crawler = tqdm(crawler, initial = 1,
				total = self.count_records(batches))

		for record in crawler:
			if current_seq == record[1]:
				current_headers.append(record[0])
			else:
				yield (current_headers, current_seq)
				current_seq = record[1]
				current_headers = [record[0]]

class KJoiner(object):
	"""docstring for KJoiner"""

	class MODES(Enum):
		"""Modes of joining.
		
		Extends:
			Enum
		
		Variables:
			UNIQUE {number} -- preserve only unique sequences
			SEQ_COUNT {number} -- count sequences
			VEC_COUNT {number} -- produce abundance vector
		"""
		UNIQUE = 1
		SEQ_COUNT = 2
		VEC_COUNT = 3
	DEFAULT_MODE = MODES.UNIQUE

	__mode = DEFAULT_MODE
	__join_function = None

	def __init__(self, mode = None):
		"""Initialize KJoiner.
		
		Keyword Arguments:
			mode {KJoiner.MODES} -- (default: {None})
		"""
		super().__init__()
		if type(mode) != type(None):
			assert mode in self.MODES
			self.__mode = mode
		self.__set_join_function()

	@property
	def mode(self):
		return self.__mode
	@mode.setter
	def mode(self, mode):
		assert mode in self.MODES
		self.__mode = mode
		self.__set_join_function()
	@property
	def join_function(self):
		return self.__join_function

	def __set_join_function(self):
		if self.mode == self.MODES.UNIQUE:
			self.__join_function = self.__join_unique
		elif self.mode == self.MODES.SEQ_COUNT:
			self.__join_function = self.__join_sequence_count
		elif self.mode == self.MODES.VEC_COUNT:
			self.__join_function = self.__join_vector_count

	def __join_unique(self, headers, seq, OH, **kwargs):
		"""Perform unique joining.
		
		Retains only unique records.
		
		Arguments:
			OH {io.TextIOWrapper} -- buffer to ouptut file
			headers {list} -- list of headers with seq
			seq {str} -- sequence
		"""
		if len(headers) == 1:
			batch = (headers[0], seq)
			OH.write(">%s\n%s\n" % batch)
			return(batch)

	def __join_sequence_count(self, headers, seq, OH, **kwargs):
		"""Perform sequence counting through joining.
		
		Counts sequence occurrences.
		
		Arguments:
			OH {io.TextIOWrapper} -- buffer to output file
			headers {list} -- list of headers with seq
			seq {str} -- sequence
		"""
		batch = (seq, len(headers))
		OH.write("%s\t%d\n" % batch)
		return(batch)

	def __join_vector_count(self, headers, seq, OH, vector, **kwargs):
		"""Generate abundance vectors through joining.
		
		Arguments:
			OH {io.TextIOWrapper} -- buffer to output file
			headers {list} -- list of headers with seq
			seq {str} -- sequence
			vector {AbundanceVector}
		"""
		regexp = re.compile(
			r'^(?P<name>[a-zA-Z0-9\\.]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$')
		hcount = len(headers)
		for header in headers:
			m = regexp.search(header)
			name, start, end = m.group("name", "start", "end")
			vector.add_count(name, "+", int(start), hcount)

	def join(self, batches, outpath, doSort = False):
		"""Join batches.
		
		Perform k-joining of batches.
		
		Arguments:
			batches {list} -- list of Batch instances
			outpath {str} -- path to output file
		
		Keyword Arguments:
			doSort {bool} -- whether batches need to be sorted
							 (default: {False})
		"""

		kwargs = {'OH' : outpath}
		if self.mode != self.MODES.VEC_COUNT:
			kwargs['OH'] = open(outpath, "w+")
		else:
			kwargs["vector"] = AbundanceVector()

		crawler = Crawler()
		for batch in crawler.do_batch(batches):
			self.join_function(*batch, **kwargs)

		if self.mode != self.MODES.VEC_COUNT:
			kwargs['OH'].close()
		else:
			kwargs['vector'].write_to(kwargs['OH'])

class SeqCountBatcher(BatcherThreading):
	"""docstring for SeqCountBatcher"""

	_type = SequenceCount
	__doSort = False

	def __init__(self, n_batches = 10, threads = 1, parent = None):
		if type(None) != type(parent):
			assert KJoinerThreading == type(parent)
			self._threads = parent.threads
			threads = parent.threads
			self._tmp = parent.tmp
		super().__init__(threads)
		self.n_batches = n_batches

	@property
	def doSort(self):
		return self.__doSort
	@doSort.setter
	def doSort(self, doSort):
		assert type(True) == type(doSort)
		self.__doSort = doSort

	def do(self, recordBatch):
		batchList = [recordBatch[i:min(len(recordBatch), i+self.n_batches)]
			for i in range(0, len(recordBatch), self.n_batches)]

		batches = Parallel(n_jobs = self.threads, verbose = 11
			)(delayed(SeqCountBatcher.build_batch
				)(batchedRecords, self) for batchedRecords in batchList)

		self.feed_collection(batches, self.FEED_MODES.REPLACE)

	@staticmethod
	def build_batch(recordBatchList, batcher):
		crawling = Crawler()
		crawling.doSort = batcher.doSort
		crawling.verbose = False

		batch = Batch(batcher, crawling.count_records(recordBatchList))
		batch.isFasta = False
		batch.add_all((SequenceCount(seq, headers)
			for (headers, seq) in crawling.do_batch(recordBatchList)))
		batch.suffix = ".txt"
		batch.write(f = "as_text", doSort = batcher.doSort)

		return batch

	def join(self):
		crawler = Crawler()

		for batch in crawler.do_batch(self.collection):
			print((batch, len(batch)))
			if 1 < len(batch):
				import sys; sys.exit()

class KJoinerThreading(KJoiner):
	"""docstring for KJoinerThreading"""

	_tmp = None
	_threads = 1
	__batch_size = 10

	def __init__(self, mode = None):
		super().__init__(mode)

	@property
	def threads(self):
		return self._threads
	@threads.setter
	def threads(self, t):
		self._threads = check_threads(t)
	@property
	def batch_size(self):
		return self.__batch_size
	@batch_size.setter
	def batch_size(self, batch_size):
		assert type(0) == type(batch_size)
		assert batch_size >= 2
		self.__batch_size = batch_size
	@property
	def tmp(self):
		if type(None) == type(self._tmp):
			self._tmp = tempfile.TemporaryDirectory(prefix = "kmanJoin")
		return self._tmp

	def __parallel_join(self, recordBatches, outpath, doSort = False):
		batcher = SeqCountBatcher(self.batch_size, parent = self)
		batcher.doSort = doSort
		batcher.do(recordBatches)
		batcher.join()
		time.sleep(1000)

	def join(self, batches, outpath, doSort = False):
		if 1 == self.threads:
			super().join(batches, outpath, doSort)
		else:
			self.__parallel_join(batches, outpath, doSort)

class AbundanceVector(object):
	"""docstring for AbundanceVector"""

	"""{ref:{strand:np.1darray}}"""
	__data = {}

	def __init__(self):
		super().__init__()

	@property
	def data(self):
		return self.__data

	def add_count(self, ref, strand, pos, count, replace = False):
		"""Add occurrence count to the vector.
		
		Arguments:
			ref {str} -- reference record name
			strand {str} -- strand type
			pos {int} -- position on ref:strand
			count {int} -- occurrence count
		
		Keyword Arguments:
			replace {bool} -- whether to allow for replacement of non-zero
							  counts (default: {False})
		"""
		self.add_ref(ref, strand, pos + 1)
		if not replace:
			assert_msg = "cannot update a non-zero count without replace."
			assert_msg += " (%s, %s, %d, %d)" % (ref, strand, pos, count)
			assert self.__data[ref][strand][pos] == 0, assert_msg
		self.__data[ref][strand][pos] = count

	def add_ref(self, ref, strand, size):
		"""Add/resize reference:strand vector.
		
		Adds a new reference:strand vector of given size  if absent, or resizes
		the current one is size is greater than the current one.
		
		Arguments:
			ref {str} -- reference record name
			strand {str} -- strand type
			size {int} -- new size
		"""
		if not ref in self.__data.keys():
			self.__data[ref] = {}
		if not strand in self.__data[ref].keys():
			self.__data[ref][strand] = np.zeros(size)
		elif size > self.__data[ref][strand].shape[0]:
			self.__data[ref][strand].resize(size)

	def write_to(self, dirpath):
		"""Write AbundanceVectors to a folder.

		The extension is removed from dirpath before proceeding.
		
		Arguments:
			dirpath {str} -- path to output directory.
		"""
		dirpath = os.path.splitext(dirpath)[0]
		assert not os.path.isfile(dirpath)
		print('Writing output in "%s"' % dirpath)
		os.makedirs(dirpath, exist_ok = True)
		for ref in self.__data.keys():
			for strand in self.__data[ref].keys():
				fname = "%s___%s.gz" % (ref, strand)
				OH = gzip.open(os.path.join(dirpath, fname), "wb")
				for count in tqdm(self.__data[ref][strand], desc = fname):
					OH.write(b"%d\n" % count)
				OH.close()
