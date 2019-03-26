
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
'''

from enum import Enum

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

	@property
	def mode():
		return self.__mode

	def join(self, batches):
		pass
