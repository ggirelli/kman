"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import tempfile
from typing import Optional

import click  # type: ignore

from kman.batcher import BatcherThreading, FastaBatcher, load_batches
from kman.const import CONTEXT_SETTINGS
from kman.io import input_file_exists, set_tempdir
from kman.join import KJoinerThreading
from kman.scripts import arguments as args


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
