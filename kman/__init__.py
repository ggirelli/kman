"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from kman import asserts, const, io
from kman import batch, batcher
from kman import abundance, seq, join

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    pass

__all__ = [
    "__version__",
    "abundance",
    "asserts",
    "batch",
    "batcher",
    "const",
    "io",
    "join",
    "seq",
]
