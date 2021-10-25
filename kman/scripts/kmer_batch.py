"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
import gzip
from kman.const import CONTEXT_SETTINGS
from kman.batcher import BatcherThreading, FastaBatcher
from kman.io import set_tempdir
import logging
import os
import shutil
import tempfile
from tqdm import tqdm  # type: ignore


@click.command(
    name="batch",
    context_settings=CONTEXT_SETTINGS,
    help="""
Generate batches of k-mers from an INPUT fasta file.

Batches are written to an OUTPUT folder, which must be empty or non-existent.
The INPUT file can be gzipped.
""",
)
@click.argument(
    "input_path",
    metavar="INPUT",
    type=click.Path(exists=True, file_okay=True, readable=True),
)
@click.argument(
    "output_path",
    metavar="OUTPUT",
    type=click.Path(exists=False, dir_okay=True, writable=True),
)
@click.argument("k", type=click.INT)
@click.option(
    "--reverse",
    "-r",
    is_flag=True,
    help="Include also reverse-complemented sequences",
)
@click.option(
    "--scan-mode",
    "-s",
    type=click.Choice([m.name for m in list(FastaBatcher.MODE)], case_sensitive=True),
    default=FastaBatcher.MODE.KMERS.name,
    help=f"""HELP PAGE MISSING!!! Default: {FastaBatcher.MODE.KMERS.name}""",
)
@click.option(
    "--batch-size",
    "-b",
    type=click.INT,
    default=1000000,
    help="Number of k-mers per batch. Default: 1000000",
)
@click.option(
    "--batch-mode",
    "-m",
    type=click.Choice(
        [m.name for m in list(BatcherThreading.FEED_MODE)], case_sensitive=True
    ),
    default=BatcherThreading.FEED_MODE.APPEND.name,
    help=f"""HELP PAGE MISSING!!! Default: {BatcherThreading.FEED_MODE.APPEND.name}""",
)
@click.option(
    "--threads",
    "-t",
    type=click.INT,
    default=1,
    help="Number of threads for parallelization.",
)
@click.option(
    "--tmp",
    "-T",
    type=click.Path(exists=True),
    default=tempfile.gettempdir(),
    help=f"""Temporary folder path. Default: "{tempfile.gettempdir()}""",
)
@click.option(
    "--compress",
    "-C",
    is_flag=True,
    help="Compress output files.",
)
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
    if not os.path.isfile(input_path):
        raise AssertionError(f"input file not found: {input_path}")
    if os.path.isdir(output_path) and len(os.listdir(output_path)) == 0:
        raise AssertionError("output folder must be empty or non-existent.")
    set_tempdir(tmp)

    batcher = FastaBatcher(size=batch_size, threads=threads)
    batcher.mode = FastaBatcher.FEED_MODE[scan_mode]
    batcher.doReverseComplement = reverse
    batcher.do(input_path, k, BatcherThreading.FEED_MODE[batch_mode])
    os.makedirs(output_path, exist_ok=True)

    try:
        run_batching(batcher, output_path, compress)
    except IOError as e:
        logging.error(f"Unable to write to output directory '{output_path}'.\n{e}")

    logging.info("That's all! :smiley:")


def run_batching(
    batcher: FastaBatcher, output_path: str, compress: bool = False
) -> None:
    batch_list = tqdm(
        (batch for batch in batcher.collection if os.path.isfile(batch.tmp)),
        total=len(batcher.collection),
    )
    for current_batch in batch_list:
        if compress:
            with gzip.open(
                os.path.join(output_path, f"{os.path.basename(current_batch.tmp)}.gz"),
                "wb",
            ) as OH:
                with open(current_batch.tmp, "rb") as IH:
                    for line in IH:
                        OH.write(line)
        else:
            shutil.copy(current_batch.tmp, output_path)
