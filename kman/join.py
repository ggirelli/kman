
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
from kman.batch import Batch
import numpy as np
import os
import re
import time
from tqdm import tqdm

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
		assert all([type(b) == Batch for b in batches])
		if doSort:
			generators = [((r.header, r.seq) for r in b.sorted)
				for b in batches]
		else:
			generators = [((r.header, r.seq) for r in b.record_gen())
				for b in batches]

		nrecords = sum([b.current_size for b in batches])
		crawler = merge(*generators, key = lambda x: x[1])

		kwargs = {}
		OH = outpath
		if self.mode != self.MODES.VEC_COUNT:
			OH = open(outpath, "w+")
		else:
			kwargs["vector"] = AbundanceVector()

		first_record = next(crawler)
		current_seq = first_record[1]
		current_headers = [first_record[0]]
		for record in tqdm(crawler, total = nrecords):
			if current_seq == record[1]:
				current_headers.append(record[0])
			else:
				self.join_function(OH, current_headers, current_seq, **kwargs)
				current_seq = record[1]
				current_headers = [record[0]]

		if self.mode != self.MODES.VEC_COUNT:
			OH.close()
		else:
			kwargs['vector'].write_to(OH)

	def __join_unique(self, OH, headers, seq, **kwargs):
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

	def __join_sequence_count(self, OH, headers, seq, **kwargs):
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

	def __join_vector_count(self, OH, headers, seq, vector, **kwargs):
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

class KJoinerThreading(object):
	"""docstring for KJoinerThreading"""

	_threads = 1

	def __init__(self, mode = None):
		super().__init__(mode)

	@property
	def threads(self):
		return self.__threads
	@threads.setter
	def threads(self, t):
		self.__threads = check_threads(t)

	def __intermediate_counts(self, batches):
		# Join subset of batches into intermediate table with format:
		# sequence | space-separated-headers
		pass

	def __join_intermediates(self):
		# Join intermediates
		pass

	def join(self, batches, outpath, doSort = False):
		if 1 == self.threads:
			super().join(batches, outpath, doSort)
		else:
			#self.__intermediate_counts(batches)
			#self.__join_intermediates()
			#self.join_function()
			pass

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
