"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import pathlib
import tempfile
from typing import Callable, Optional

import click  # type: ignore

from kmermaid.batcher import BatcherThreading, FastaBatcher, load_batches
from kmermaid.const import CONTEXT_SETTINGS
from kmermaid.io import input_file_exists, set_tempdir
from kmermaid.join import KJoinerThreading
from kmermaid.scripts import arguments as args


def add_click_hooks(f: Callable[..., None]):
    """Click hook decorator.

    :param f: run function
    :type f: Callable[..., None]
    :return: decorated f
    :rtype: Callable[..., None]
    """

    @click.command(
        name="uniq",
        context_settings=CONTEXT_SETTINGS,
        help="Extract all k-mers that appear only once in the INPUT fasta file.",
    )
    @args.input_path()
    @args.output_path(file_okay=True)
    @args.k()
    @args.reverse()
    @args.scan_mode()
    @args.batch_size()
    @args.batch_mode()
    @args.previous_batches()
    @args.threads()
    @args.tmp()
    @args.re_sort()
    def wrapper(*args, **kwargs):
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
    threads: int = 1,
    tmp: str = tempfile.gettempdir(),
    re_sort: bool = False,
) -> None:
    """Run kmer_uniq script.

    :param input_path: input FASTA file
    :type input_path: str
    :param output_path: path to output file
    :type output_path: str
    :param k: k-mer length
    :type k: int
    :param reverse: perform reverse-complement operation, defaults to False
    :type reverse: bool
    :param scan_mode: scanning mode, defaults to FastaBatcher.MODE.KMERS.name
    :type scan_mode: str
    :param batch_size: records per batch, defaults to 1000000
    :type batch_size: int
    :param batch_mode: batching mode, defaults to BatcherThreading.FEED_MODE.APPEND.name
    :type batch_mode: str
    :param previous_batches: path to folder with previous batches, defaults to None
    :type previous_batches: Optional[str]
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

    get_joiner(len(batches), threads).join(batches, output_path)

    logging.info("That's all! :smiley:")


def get_joiner(n_batches: int, threads: int = 1) -> KJoinerThreading:
    """Instantiate k-joiner for uniquing.

    :param n_batches: number of input batches.
    :type n_batches: int
    :param threads: for parallelization, defaults to 1
    :type threads: int
    :return: updated joiner instance.
    :rtype: KJoinerThreading
    """
    joiner = KJoinerThreading()
    joiner.threads = threads
    joiner.batch_size = max(2, int(n_batches / threads))
    return joiner
