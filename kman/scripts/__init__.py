"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from kman.scripts import arguments
from kman.scripts import kmer, kmer_batch, kmer_count, kmer_uniq

import logging
from rich.logging import RichHandler  # type: ignore

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

__all__ = ["arguments", "kmer", "kmer_batch", "kmer_count", "kmer_uniq"]
