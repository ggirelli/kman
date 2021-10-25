"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import resource
import tempfile
from typing import Optional

import click  # type: ignore

from kman.batcher import BatcherThreading, FastaBatcher, load_batches
from kman.const import CONTEXT_SETTINGS
from kman.io import input_file_exists, set_tempdir
from kman.join import KJoiner, KJoinerThreading
from kman.scripts import arguments as args


@click.command(
    name="count",
    context_settings=CONTEXT_SETTINGS,
    help="""
Count occurrences of all k-mers from INPUT.

\b
Counting modes:
       SEQ_COUNT a tabulation-separated table with sequence and count
       VEC_COUNT a vector for each record in the input fasta, containing the
                 occurrence count of each k-mer (start position-bound)
VEC_COUNT_MASKED a vector for each record in the input fasta, containing the
                 occurrence count of each k-mer (start position-bound). Only
                 occurrences in other records are counted.

When using a vector (VEC_*) mode, the output will be stored in a folder with the
name of the specified OUTPUT path, after removing the extension.
The INPUT file can be gzipped.
""",
)
@args.input_path()
@args.output_path(file_okay=True)
@args.k()
@args.reverse()
@args.scan_mode()
@args.batch_size()
@args.batch_mode()
@args.previous_batches()
@args.count_mode()
@args.memory_mode()
@args.threads()
@args.tmp()
@args.re_sort()
def run(
    input_path: str,
    output_path: str,
    k: int,
    reverse: bool = False,
    scan_mode: str = FastaBatcher.MODE.KMERS.name,
    batch_size: int = 1000000,
    batch_mode: str = BatcherThreading.FEED_MODE.APPEND.name,
    previous_batches: Optional[str] = None,
    count_mode: str = KJoiner.MODE.SEQ_COUNT.name,
    memory_mode: str = KJoiner.MEMORY.NORMAL.name,
    threads: int = 1,
    tmp: str = tempfile.gettempdir(),
    re_sort: bool = False,
) -> None:
    input_file_exists(input_path)
    set_tempdir(tmp)

    if previous_batches is not None:
        batches = load_batches(previous_batches, threads, re_sort)
    else:
        batches = (
            FastaBatcher(
                scan_mode=FastaBatcher.MODE[scan_mode],
                reverse=reverse,
                size=batch_size,
                threads=threads,
            )
            .do(input_path, k, BatcherThreading.FEED_MODE[batch_mode])
            .collection
        )

    prep_joiner(
        KJoinerThreading(KJoiner.MODE[count_mode], KJoiner.MEMORY[memory_mode]),
        len(batches),
        threads,
    ).join(batches, output_path)

    logging.info("That's all! :smiley:")


def prep_joiner(
    joiner: KJoinerThreading, n_batches: int, threads: int = 1
) -> KJoinerThreading:
    """Instantiate k-way joiner for counting.

    Args:
        joiner (KJoinerThreading): pre-instantiated joiner.
        n_batches (int): number of input batches.
        threads (int, optional): for parallelization. Defaults to 1.

    Returns:
        KJoinerThreading
    """
    joiner.threads = threads
    joiner.batch_size = max(2, int(n_batches / threads))
    rlim_min, rlim_max = resource.getrlimit(resource.RLIMIT_NOFILE)
    joiner.batch_size = min(joiner.batch_size, rlim_max)
    resource.setrlimit(
        resource.RLIMIT_NOFILE, (max(rlim_min, joiner.batch_size), rlim_max)
    )
    return joiner
