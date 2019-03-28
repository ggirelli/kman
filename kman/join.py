
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
'''

from enum import Enum
from heapq import merge
from kman.batch import Batch
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
		if self.mode == self.MODES.UNIQUE:
			self.__join_function = self.__join_unique
		elif self.mode == self.MODES.SEQ_COUNT:
			self.__join_function = self.__join_sequence_count
		elif self.mode == self.MODES.VEC_COUNT:
			self.__join_function = self.__join_vector_count

	@property
	def mode(self):
		return self.__mode
	@property
	def join_function(self):
		return self.__join_function

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

		with open(outpath, "w+") as OH:
			first_record = next(crawler)
			current_seq = first_record[1]
			current_headers = [first_record[0]]
			for record in tqdm(crawler, total = nrecords):
				if current_seq == record[1]:
					current_headers.append(record[0])
				else:
					self.join_function(OH, (current_seq, current_headers))

	def __join_unique(self, OH, batch):
		print(batch)
		import sys; sys.exit()
		pass

	def __join_sequence_count(self, OH, batch):
		print(batch)
		import sys; sys.exit()
		pass

	def __join_vector_count(self, OH, batch):
		print(batch)
		import sys; sys.exit()
		pass
