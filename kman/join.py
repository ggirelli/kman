
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
from itertools import chain
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
		"""Count records across batches.
		
		Arguments:
			batches {list} -- list of Batches
		
		Returns:
			int -- number of records
		"""
		return sum([b.current_size for b in batches])

	def do_records(self, batches):
		"""Crawl through the batches.

		Produces a generator function that yields (r.header, r.seq) for each
		record r across batches. Batches are crawled through their r.record_gen
		if self.doSort==False, otherwise using r.sorted.
		
		Arguments:
			batches {list} -- list of Batches
		
		Returns:
			generator -- record generator
		"""
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
		
		yield (current_headers, current_seq)

class KJoiner(object):
	"""docstring for KJoiner"""

	class MODE(Enum):
		"""Modes of joining.
		
		Extends:
			Enum
		
		Variables:
			UNIQUE {number} -- preserve only unique sequences
			SEQ_COUNT {number} -- count sequences
			VEC_COUNT {number} -- produce abundance vector
			VEC_COUNT_MASKED {number} -- produce abundance vector w/o counting
										 abundance in same reference record
		"""
		UNIQUE = 1
		SEQ_COUNT = 2
		VEC_COUNT = 3
		VEC_COUNT_MASKED = 4
	DEFAULT_MODE = MODE.UNIQUE

	__mode = DEFAULT_MODE
	__join_function = None

	def __init__(self, mode = None):
		"""Initialize KJoiner.
		
		Keyword Arguments:
			mode {KJoiner.MODE} -- (default: {None})
		"""
		super().__init__()
		if type(mode) != type(None):
			assert mode in self.MODE
			self.__mode = mode
		self.__set_join_function()

	@property
	def mode(self):
		return self.__mode
	@mode.setter
	def mode(self, mode):
		assert mode in self.MODE
		self.__mode = mode
		self.__set_join_function()
	@property
	def join_function(self):
		return self.__join_function

	def __set_join_function(self):
		if self.mode == self.MODE.UNIQUE:
			self.__join_function = self.join_unique
		elif self.mode == self.MODE.SEQ_COUNT:
			self.__join_function = self.join_sequence_count
		elif self.mode == self.MODE.VEC_COUNT:
			self.__join_function = self.join_vector_count
		elif self.mode == self.MODE.VEC_COUNT_MASKED:
			self.__join_function = self.join_vector_count_masked

	@staticmethod
	def join_unique(headers, seq, OH, **kwargs):
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

	@staticmethod
	def join_sequence_count(headers, seq, OH, **kwargs):
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

	@staticmethod
	def join_vector_count(headers, seq, OH, vector, **kwargs):
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

	@staticmethod
	def join_vector_count_masked(headers, seq, OH, vector, **kwargs):
		"""Generate abundance vectors through joining.
		
		Arguments:
			OH {io.TextIOWrapper} -- buffer to output file
			headers {list} -- list of headers with seq
			seq {str} -- sequence
			vector {AbundanceVector}
		"""
		regexp = re.compile(
			r'^(?P<name>[a-zA-Z0-9\\.]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$')

		headers = [regexp.search(h).group("name", "start") for h in headers]
		if not 1 == len(headers):
			refList, refCounts = np.unique([h[0] for h in headers],
				return_counts = True)

			if not 1 == len(refList):
				for name, start in headers:
					hcount = refCounts[refList != name].sum()
					vector.add_count(name, "+", int(start), hcount)

	def _pre_join(self, outpath):
		kwargs = {'OH' : outpath}
		if not self.mode.name.startswith("VEC_"):
			kwargs['OH'] = open(outpath, "w+")
		else:
			kwargs["vector"] = AbundanceVector()
		return kwargs

	def _post_join(self, **kwargs):
		if not self.mode.name.startswith("VEC_"):
			kwargs['OH'].close()
		else:
			kwargs['vector'].write_to(kwargs['OH'])

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

		kwargs = self._pre_join((outpath))

		crawler = Crawler()
		for batch in crawler.do_batch(batches):
			self.join_function(*batch, **kwargs)

		self._post_join(**kwargs)

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

		self.feed_collection(batches, self.FEED_MODE.REPLACE)

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

	def join(self, fjoin, **kwargs):
		crawler = Crawler()
		for headers, seq in crawler.do_batch(self.collection):
			headers = list(chain(*headers))
			fjoin(headers, seq, **kwargs)

class KJoinerThreading(KJoiner):
	"""docstring for KJoinerThreading"""

	_tmp = None
	_threads = 1
	__batch_size = 10
	__doSort = False

	def __init__(self, mode = None):
		super().__init__(mode)

	@property
	def doSort(self):
		return self.__doSort
	@doSort.setter
	def doSort(self, doSort):
		assert type(True) == type(doSort)
		self.__doSort = doSort
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

	def __parallel_join(self, recordBatches, outpath):
		kwargs = self._pre_join(outpath)

		batcher = SeqCountBatcher(self.batch_size, parent = self)
		batcher.doSort = self.doSort
		batcher.do(recordBatches)
		batcher.join(self.join_function, **kwargs)
		
		self._post_join(**kwargs)

	def join(self, batches, outpath):
		if 1 == self.threads:
			super().join(batches, outpath, self.doSort)
		else:
			self.__parallel_join(batches, outpath)

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
