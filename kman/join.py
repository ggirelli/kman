
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
'''

from enum import Enum

class KJoiner(object):
	"""docstring for KJoiner"""

	class MODE(Enum):
		UNIQUE = 1
		SEQ_COUNT = 2
		VCOUNT = 3

	def __init__(self):
		super(KJoiner, self).__init__()
