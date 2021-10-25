"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from kman.const import CONTEXT_SETTINGS
from kman.batcher import BatcherThreading, FastaBatcher
from kman.io import set_tempdir
from kman.join import KJoiner, KJoinerThreading
import logging
import os
import resource
import tempfile
from typing import Optional


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
@click.argument(
    "input_path",
    metavar="INPUT",
    type=click.Path(exists=True, file_okay=True, readable=True),
)
@click.argument(
    "output_path",
    metavar="OUTPUT",
    type=click.Path(exists=False, file_okay=True, writable=True),
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
    "--previous-batches",
    "-B",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="Path to folder with previously generated batches.",
)
@click.option(
    "--count-mode",
    "-m",
    type=click.Choice(
        [m.name for m in list(KJoiner.MODE) if "COUNT" in m.name], case_sensitive=True
    ),
    default=KJoiner.MODE.SEQ_COUNT.name,
    help=f"""Default: "{KJoiner.MODE.SEQ_COUNT.name}""",
)
@click.option(
    "--memory-mode",
    "-M",
    type=click.Choice([m.name for m in list(KJoiner.MEMORY)], case_sensitive=True),
    default=KJoiner.MEMORY.NORMAL.name,
    help=f"""HELP PAGE MISSING!!! Default: "{KJoiner.MEMORY.NORMAL.name}""",
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
    "--re-sort",
    "-R",
    is_flag=True,
    help="Force batch re-sorting, when loaded with -B.",
)
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
    if not os.path.isfile(input_path):
        raise AssertionError(f"input file not found: {input_path}")
    set_tempdir(tmp)

    if previous_batches is not None:
        if not os.path.isdir(previous_batches) or len(os.listdir(previous_batches)) > 0:
            raise AssertionError(
                f"folder with previous batches empty or not found: {previous_batches}"
            )
        logging.info(f"Loading previous batches from '{previous_batches}'...")
        batches = BatcherThreading.from_files(previous_batches, threads, reSort=re_sort)
    else:
        batcher = FastaBatcher(size=batch_size, threads=threads)
        batcher.mode = FastaBatcher.FEED_MODE[scan_mode]
        batcher.doReverseComplement = reverse
        batcher.do(input_path, k, BatcherThreading.FEED_MODE[batch_mode])
        batches = batcher.collection

    joiner = KJoinerThreading(KJoiner.MODE[count_mode], KJoiner.MEMORY[memory_mode])
    joiner.threads = threads

    joiner.batch_size = max(2, int(len(batches) / threads))
    rlimits = resource.getrlimit(resource.RLIMIT_NOFILE)
    joiner.batch_size = min(joiner.batch_size, rlimits[1])
    resource.setrlimit(
        resource.RLIMIT_NOFILE, (max(rlimits[0], joiner.batch_size), rlimits[1])
    )

    joiner.join(batches, output_path)
    logging.info("That's all! :smiley:")
