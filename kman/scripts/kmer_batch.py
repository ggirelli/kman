"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
import gzip
from kman.asserts import enable_rich_assert
from kman.batcher import BatcherThreading, FastaBatcher
from kman.scripts import arguments as ap
import logging
import os
import shutil
import tempfile
from tqdm import tqdm  # type: ignore


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "batch",
        description="""
Generate batches of k-mers from a fasta file.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Create k-mer batches from a fasta file.",
    )

    parser.add_argument(
        "input",
        type=str,
        help="""
    Path to input fasta file. Can be gzipped (ending in ".gz")""",
    )
    parser.add_argument(
        "output",
        type=str,
        help="""
        Path to output folder, which must be empty or non-existent.""",
    )
    parser.add_argument(
        "k",
        type=int,
        help="""Oligonucleotide (substring) length in nucleotides.""",
    )

    parser.add_argument(
        "-R",
        dest="do_reverse",
        action="store_const",
        const=True,
        default=False,
        help="""Reverse complement sequences.""",
    )

    advanced = parser.add_argument_group("advanced arguments")
    advanced.add_argument(
        "-s",
        type=str,
        help=f'''Choose scanning mode. See description for more details.
        Default: "{FastaBatcher.MODE.KMERS.name}"''',
        default=FastaBatcher.MODE.KMERS.name,
        choices=[m.name for m in list(FastaBatcher.MODE)],
    )
    advanced.add_argument(
        "-b", type=int, default=1e6, help="""Number of kmers per batch. Default: 1e6"""
    )
    advanced.add_argument(
        "-m",
        type=str,
        help=f'''Choose batching mode. See description for more details.
        Default: "{BatcherThreading.FEED_MODE.APPEND.name}"''',
        default=BatcherThreading.FEED_MODE.APPEND.name,
        choices=[m.name for m in list(BatcherThreading.FEED_MODE)],
    )
    advanced.add_argument(
        "-t", type=int, default=1, help="""Number of threads for parallelization."""
    )
    advanced.add_argument(
        "-T",
        type=str,
        default=tempfile.gettempdir(),
        help=f'''Temporary folder path. Default: "{tempfile.gettempdir()}"''',
    )
    advanced.add_argument(
        "-C",
        dest="do_compress",
        action="store_const",
        const=True,
        default=False,
        help="""Compress output files.""",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    args.o = args.output
    args.m = BatcherThreading.FEED_MODE[args.m]
    args.s = FastaBatcher.MODE[args.s]

    assert os.path.isfile(
        args.input
    ), f"path to an existing fasta file expected, file not found: '{args.input}'"

    if os.path.isdir(args.o):
        assert 0 == len(
            os.listdir(args.o)
        ), "output folder must be empty or non-existent."

    return args


def run_batching(args: argparse.Namespace, batcher: FastaBatcher) -> None:
    batchList = tqdm([b for b in batcher.collection if os.path.isfile(b.tmp)])
    if args.do_compress:
        for b in batchList:
            gzname = os.path.join(args.o, f"{os.path.basename(b.tmp)}.gz")
            OH = gzip.open(gzname, "wb")
            with open(b.tmp, "rb") as IH:
                for line in IH:
                    OH.write(line)
            OH.close()
    else:
        for b in batchList:
            shutil.copy(b.tmp, args.o)


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    if not os.path.isdir(args.T):
        os.makedirs(args.T, exist_ok=True)
    tempfile.tempdir = args.T

    batcher = FastaBatcher(size=args.b, threads=args.t)
    batcher.mode = args.s
    batcher.doReverseComplement = args.do_reverse
    batcher.do(args.input, args.k, args.m)
    os.makedirs(args.o, exist_ok=True)

    try:
        run_batching(args, batcher)
    except IOError as e:
        logging.error(f"Unable to write to output directory '{args.o}'.\n{e}")

    logging.info("That's all! :smiley:")
