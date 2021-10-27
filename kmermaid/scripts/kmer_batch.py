"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import os
import pathlib
import tempfile
from typing import Callable

import click  # type: ignore

from kmermaid.batcher import BatcherThreading, FastaBatcher
from kmermaid.const import CONTEXT_SETTINGS
from kmermaid.io import copy_batches, input_file_exists, set_tempdir
from kmermaid.scripts import arguments as args


def add_click_hooks(f: Callable[..., None]):
    """Click hook decorator.

    :param f: run function
    :type f: Callable[..., None]
    :return: decorated f
    :rtype: Callable[..., None]
    """

    @click.command(
        name="batch",
        context_settings=CONTEXT_SETTINGS,
        help="""
    Generate batches of k-mers from an INPUT fasta file.

    Batches are written to an OUTPUT folder, which must be empty or non-existent.
    The INPUT file can be gzipped.
    """,
    )
    @args.input_path()
    @args.output_path(dir_okay=True)
    @args.k()
    @args.reverse()
    @args.scan_mode()
    @args.batch_size()
    @args.batch_mode()
    @args.threads()
    @args.tmp()
    @args.compress()
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
    threads: int = 1,
    tmp: str = tempfile.gettempdir(),
    compress: bool = False,
) -> None:
    """Run kmer_batch script.

    :param input_path: path to input FASTA
    :type input_path: str
    :param output_path: path to output folder
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
    :param threads: for parallelization, defaults to 1
    :type threads: int
    :param tmp: path to temporary folder, defaults to tempfile.gettempdir()
    :type tmp: str
    :param compress: compress output, defaults to False
    :type compress: bool
    """
    prepare_run(input_path, output_path, tmp)

    try:
        copy_batches(
            FastaBatcher(
                scan_mode=FastaBatcher.MODE[scan_mode],
                reverse=reverse,
                size=batch_size,
                threads=threads,
            )
            .do(pathlib.Path(input_path), k, BatcherThreading.FEED_MODE[batch_mode])
            .collection,
            output_path,
            compress,
        )
    except IOError as e:
        logging.error(f"Unable to write to output directory '{output_path}'.\n{e}")

    logging.info("That's all! :smiley:")


def prepare_run(input_path: str, output_path: str, tmp: str) -> None:
    """Prepare output folders and check input before running.

    :param input_path: path to input FASTA.
    :type input_path: str
    :param output_path: path to output folder.
    :type output_path: str
    :param tmp: path to temporary folder.
    :type tmp: str
    :raises AssertionError: if output folder exists or is not empty.
    """
    input_file_exists(input_path)
    if os.path.isdir(output_path) and len(os.listdir(output_path)) != 0:
        raise AssertionError("output folder must be empty or non-existent.")
    set_tempdir(tmp)
    os.makedirs(output_path, exist_ok=True)
