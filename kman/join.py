
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
		VCOUNT = 3
	DEFAULT_MODE = MODES.UNIQUE

	__mode = DEFAULT_MODE

	def __init__(self, mode = None):
		super(KJoiner, self).__init__()
		if type(mode) != type(None):
			assert mode in self.MODES
			self.__mode = mode

	@property
	def mode():
		return self.__mode

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
			for record in tqdm(crawler, total = nrecords):
				time.sleep(0.01)
				pass
