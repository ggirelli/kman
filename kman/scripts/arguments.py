"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from kman.batcher import BatcherThreading, FastaBatcher
from kman.join import KJoiner
import tempfile


def input_path() -> click.Argument:
    return click.argument(
        "input_path",
        metavar="INPUT",
        type=click.Path(exists=True, file_okay=True, readable=True),
    )


def output_path(file_okay=False, dir_okay=False) -> click.Argument:
    return click.argument(
        "output_path",
        metavar="OUTPUT",
        type=click.Path(
            exists=False, file_okay=file_okay, dir_okay=dir_okay, writable=True
        ),
    )


def k() -> click.Argument:
    return click.argument("k", type=click.INT)


def reverse() -> click.Option:
    return click.option(
        "--reverse",
        "-r",
        is_flag=True,
        help="Include also reverse-complemented sequences",
    )


def scan_mode() -> click.Option:
    return click.option(
        "--scan-mode",
        "-s",
        type=click.Choice(
            [m.name for m in list(FastaBatcher.MODE)], case_sensitive=True
        ),
        default=FastaBatcher.MODE.KMERS.name,
        help=f"""HELP PAGE MISSING!!! Default: {FastaBatcher.MODE.KMERS.name}""",
    )


def batch_size() -> click.Option:
    return click.option(
        "--batch-size",
        "-b",
        type=click.INT,
        default=1000000,
        help="Number of k-mers per batch. Default: 1000000",
    )


def batch_mode() -> click.Option:
    return click.option(
        "--batch-mode",
        "-m",
        type=click.Choice(
            [m.name for m in list(BatcherThreading.FEED_MODE)], case_sensitive=True
        ),
        default=BatcherThreading.FEED_MODE.APPEND.name,
        help=f"""HELP PAGE MISSING!!!
        Default: {BatcherThreading.FEED_MODE.APPEND.name}""",
    )


def previous_batches() -> click.Option:
    return click.option(
        "--previous-batches",
        "-B",
        type=click.Path(exists=True, dir_okay=True, readable=True),
        help="Path to folder with previously generated batches.",
    )


def count_mode() -> click.Option:
    return click.option(
        "--count-mode",
        "-m",
        type=click.Choice(
            [m.name for m in list(KJoiner.MODE) if "COUNT" in m.name],
            case_sensitive=True,
        ),
        default=KJoiner.MODE.SEQ_COUNT.name,
        help=f"""Default: "{KJoiner.MODE.SEQ_COUNT.name}""",
    )


def memory_mode() -> click.Option:
    return click.option(
        "--memory-mode",
        "-M",
        type=click.Choice([m.name for m in list(KJoiner.MEMORY)], case_sensitive=True),
        default=KJoiner.MEMORY.NORMAL.name,
        help=f"""HELP PAGE MISSING!!! Default: "{KJoiner.MEMORY.NORMAL.name}""",
    )


def threads() -> click.Option:
    return click.option(
        "--threads",
        "-t",
        type=click.INT,
        default=1,
        help="Number of threads for parallelization.",
    )


def tmp() -> click.Option:
    return click.option(
        "--tmp",
        "-T",
        type=click.Path(exists=True),
        default=tempfile.gettempdir(),
        help=f"""Temporary folder path. Default: "{tempfile.gettempdir()}""",
    )


def compress() -> click.Option:
    return click.option(
        "--compress",
        "-C",
        is_flag=True,
        help="Compress output files.",
    )


def re_sort() -> click.Option:
    return click.option(
        "--re-sort",
        "-R",
        is_flag=True,
        help="Force batch re-sorting, when loaded with -B.",
    )
