"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import pathlib
import resource
import tempfile
from typing import Callable, Optional

import click  # type: ignore

from kmermaid.batcher import BatcherThreading, FastaBatcher, load_batches
from kmermaid.const import CONTEXT_SETTINGS
from kmermaid.io import input_file_exists, set_tempdir
from kmermaid.join import KJoiner, KJoinerThreading
from kmermaid.scripts import arguments as args


def add_click_hooks(f: Callable[..., None]):
    """Click hook decorator.

    :param f: run function
    :type f: Callable[..., None]
    :return: decorated f
    :rtype: Callable[..., None]
    """

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
    def wrapper(*args, **kwargs):
        """Decorator wrapper."""
        f(*args, **kwargs)

    return wrapper


@add_click_hooks
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
    """Run kmer_count script.

    :param input_path: input FASTA file
    :type input_path: str
    :param output_path: output file
    :type output_path: str
    :param k: k-mer length
    :type k: int
    :param reverse: perform reverse complement operation, defaults to False
    :type reverse: bool
    :param scan_mode: scanning mode, defaults to FastaBatcher.MODE.KMERS.name
    :type scan_mode: str
    :param batch_size: records per batch, defaults to 1000000
    :type batch_size: int
    :param batch_mode: batching mode, defaults to BatcherThreading.FEED_MODE.APPEND.name
    :type batch_mode: str
    :param previous_batches: path to folder with previous batches, defaults to None
    :type previous_batches: Optional[str]
    :param count_mode: counting mode, defaults to KJoiner.MODE.SEQ_COUNT.name
    :type count_mode: str
    :param memory_mode: memory handling mode, defaults to KJoiner.MEMORY.NORMAL.name
    :type memory_mode: str
    :param threads: for parallelization, defaults to 1
    :type threads: int
    :param tmp: path to temporary folder, defaults to tempfile.gettempdir()
    :type tmp: str
    :param re_sort: re-sort previous batches, defaults to False
    :type re_sort: bool
    """
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
            .do(pathlib.Path(input_path), k, BatcherThreading.FEED_MODE[batch_mode])
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
    """Instantiate k-joiner for counting.

    :param joiner: pre-instantiated joiner.
    :type joiner: KJoinerThreading
    :param n_batches: number of batches in input.
    :type n_batches: int
    :param threads: for parallelization, defaults to 1
    :type threads: int
    :return: updated joiner instance.
    :rtype: KJoinerThreading
    """
    joiner.threads = threads
    joiner.batch_size = max(2, n_batches // threads)
    rlim_min, rlim_max = resource.getrlimit(resource.RLIMIT_NOFILE)
    joiner.batch_size = min(joiner.batch_size, rlim_max)
    resource.setrlimit(
        resource.RLIMIT_NOFILE, (max(rlim_min, joiner.batch_size), rlim_max)
    )
    return joiner
