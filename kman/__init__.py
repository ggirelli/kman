"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from importlib.metadata import version

from kman import abundance, batch, batcher, const, io, join, seq

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
