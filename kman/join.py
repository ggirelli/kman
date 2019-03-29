
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
'''

from enum import Enum
from heapq import merge
import io
from kman.batch import Batch
import numpy as np
import re
import time
from tqdm import tqdm

class KJoiner(object):
	"""docstring for KJoiner"""

	class MODES(Enum):
		UNIQUE = 1
		SEQ_COUNT = 2
		VEC_COUNT = 3
	DEFAULT_MODE = MODES.UNIQUE

	__mode = DEFAULT_MODE
	__join_function = None

	def __init__(self, mode = None):
		super(KJoiner, self).__init__()
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
		assert all([type(b) == Batch for b in batches])
		if doSort:
			generators = [((r.header, r.seq) for r in b.sorted)
				for b in batches]
		else:
			generators = [((r.header, r.seq) for r in b.record_gen())
				for b in batches]

		nrecords = sum([b.current_size for b in batches])
		crawler = merge(*generators, key = lambda x: x[1])

		OH = outpath
		if self.mode != self.MODES.VEC_COUNT:
			OH = open(outpath, "w+")

		first_record = next(crawler)
		current_seq = first_record[1]
		current_headers = [first_record[0]]
		for record in tqdm(crawler, total = nrecords):
			if current_seq == record[1]:
				current_headers.append(record[0])
			else:
				self.join_function(OH, current_headers, current_seq)
				current_seq = record[1]
				current_headers = [record[0]]

		if self.mode != self.MODES.VEC_COUNT:
			OH.close()

	def __join_unique(self, OH, headers, seq, **kwargs):
		if len(headers) == 1:
			batch = (headers[0], seq)
			OH.write(">%s\n%s\n" % batch)
			return(batch)

	def __join_sequence_count(self, OH, headers, seq, **kwargs):
		batch = (seq, len(headers))
		OH.write("%s\t%d\n" % batch)
		return(batch)

	def __join_vector_count(self, OH, headers, seq, **kwargs):
		regexp = re.compile(
			r'(?P<name>[a-zA-Z0-9\\.]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)')
		for header in headers:
			m = regexp.search(header)
			record_name, start, end = m.group("name", "start", "end")
			import sys; sys.exit()
		pass

class RecordStorage(io.FileIO):
	"""docstring for RecordStorage"""

	__n_lines = 0
	__line_locations = {}

	def __init__(self, name, delim = "\t", mode = "a+b", **kwargs):
		super(RecordStorage, self).__init__(name, mode, **kwargs)
		self.__delim = delim

		self.seek(0, 2)
		file_end_location = self.tell()
		self.seek(0)
		if 0 != file_end_location:
			while self.tell() != file_end_location:
				current_location = self.tell()
				line = next(self).decode("utf-8")
				current_key = int(line.strip().split(self.__delim)[0])
				self.__line_locations[current_key] = current_location
			self.seek(0)

	def add_count(self, pos, count):
		assert_msg = "a count for position %d has already been stored." % pos
		assert not pos in self.__line_locations.keys(), assert_msg
		self.seek(0, 2)
		self.__line_locations[pos] = self.tell()
		self.write(bytes("%d%s%d\n" % (pos, self.__delim, count), "utf-8"))

	def get_count(self, pos):
		self.seek(self.__line_locations[pos])
		line = next(self).decode("utf-8")
		return line.strip().split(self.__delim)[1]
