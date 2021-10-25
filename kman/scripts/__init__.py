"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging

from rich.logging import RichHandler  # type: ignore

from kman.scripts import arguments, kmer, kmer_batch, kmer_count, kmer_uniq

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

__all__ = ["arguments", "kmer", "kmer_batch", "kmer_count", "kmer_uniq"]
