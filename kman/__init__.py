"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from kman import const, io
from kman import batch, batcher
from kman import abundance, seq, join
from importlib.metadata import version

try:
    __version__ = version(__name__)
except Exception as e:
    raise e

__all__ = [
    "__version__",
    "abundance",
    "batch",
    "batcher",
    "const",
    "io",
    "join",
    "seq",
]
