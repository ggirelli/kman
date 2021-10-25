"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from kman.const import CONTEXT_SETTINGS
from kman.batcher import BatcherThreading, FastaBatcher
from kman.io import copy_batches, input_file_exists, set_tempdir
from kman.scripts import arguments as args
import logging
import os
import tempfile


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
    prepare_run(input_path, output_path, tmp)

    try:
        copy_batches(
            FastaBatcher(
                scan_mode=FastaBatcher.MODE[scan_mode],
                reverse=reverse,
                size=batch_size,
                threads=threads,
            )
            .do(input_path, k, BatcherThreading.FEED_MODE[batch_mode])
            .collection,
            output_path,
            compress,
        )
    except IOError as e:
        logging.error(f"Unable to write to output directory '{output_path}'.\n{e}")

    logging.info("That's all! :smiley:")


def prepare_run(input_path: str, output_path: str, tmp: str) -> None:
    """Prepare output folders and checks input before running.

    Args:
        input_path (str): path to input FASTA
        output_path (str): path to output folder
        tmp (str): path to temporary folder

    Raises:
        AssertionError: if output folder exists or is not empty
    """
    input_file_exists(input_path)
    if os.path.isdir(output_path) and len(os.listdir(output_path)) == 0:
        raise AssertionError("output folder must be empty or non-existent.")
    set_tempdir(tmp)
    os.makedirs(output_path, exist_ok=True)
