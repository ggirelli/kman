"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from kman.asserts import enable_rich_assert
from kman.batcher import BatcherThreading, FastaBatcher
from kman.join import KJoinerThreading
from kman.scripts import arguments as ap
import logging
import os
import tempfile


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "uniq",
        description="""
Extract all k-mer (i.e., k-characters substrings) that appear in the input
sequence only once.
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
        Path to output fasta file.""",
    )
    parser.add_argument(
        "k",
        type=int,
        required=True,
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
        "-B",
        type=str,
        help="""Path to folder with previously
        generated batches. Useful to skip the batching step.""",
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
        "--resort",
        dest="do_resort",
        action="store_const",
        const=True,
        default=False,
        help="""Force batch re-sorting, when loaded
        with -B.""",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert os.path.isfile(
        args.input
    ), f"path to an existing fasta file expected, file not found: '{args.input}'"

    if not os.path.isdir(args.T):
        os.makedirs(args.T, exist_ok=True)
    tempfile.tempdir = args.T

    args.m = BatcherThreading.FEED_MODE[args.m]
    args.s = FastaBatcher.MODE[args.s]

    if args.B is not None:
        assert os.path.isdir(args.B), f"batching output folder not found: '{args.B}'"
        assert 0 < len(
            os.listdir(args.B)
        ), f"the batching output folder is empty: '{args.B}'"

    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    if args.B is not None:
        logging.info("Loading previously generated batches from '%s'..." % args.B)
        batches = BatcherThreading.from_files(args.B, args.t, reSort=args.do_resort)
    else:
        batcher = FastaBatcher(size=args.b, threads=args.t)
        batcher.mode = args.s
        batcher.do(args.input, args.k, args.m)
        batcher.doReverseComplement = args.do_reverse
        batches = batcher.collection

    joiner = KJoinerThreading()
    joiner.threads = args.t
    joiner.batch_size = max(2, int(len(batches) / args.t))
    joiner.join(batches, args.output)
    logging.info("That's all! :smiley:")
