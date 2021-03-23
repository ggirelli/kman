"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from kman.asserts import enable_rich_assert
from kman.batcher import BatcherThreading, FastaBatcher
from kman.join import KJoiner, KJoinerThreading
from kman.scripts import arguments as ap
import logging
import os
import resource
import tempfile


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "count",
        description="""
Count occurrences of all k-mer (i.e., k-characters substrings) that appear in
the input sequence. Depending on the mode, a different output is generated:

         SEQ_COUNT a tabulation-separated table with sequence and count

         VEC_COUNT a vector for each record in the input fasta, containing the
                   occurrence count of each k-mer (start position-bound)

  VEC_COUNT_MASKED a vector for each record in the input fasta, containing the
                   occurrence count of each k-mer (start position-bound). Only
                   occurrences in other records are counted.

When using a VEC_* mode, the output will be stored in a folder with the name of
the specified output, after removing the extension.
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
        Path to output tsv file.""",
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
        "-m",
        type=str,
        help=f'''Choose batching mode. See description for more details.
        Default: "{BatcherThreading.FEED_MODE.APPEND.name}"''',
        default=BatcherThreading.FEED_MODE.APPEND.name,
        choices=[m.name for m in list(BatcherThreading.FEED_MODE)],
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
        "-c",
        type=str,
        help=f'''Choose counting mode. See description for more details.
        Default: "{KJoiner.MODE.SEQ_COUNT.name}"''',
        default=KJoiner.MODE.SEQ_COUNT.name,
        choices=[m.name for m in list(KJoiner.MODE) if "COUNT" in m.name],
    )
    advanced.add_argument(
        "-M",
        type=str,
        help=f'''Choose memory mode. See description for more details.
        Default: "{KJoiner.MEMORY.NORMAL.name}"''',
        default=KJoiner.MEMORY.NORMAL.name,
        choices=[m.name for m in list(KJoiner.MEMORY)],
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

    args.s = FastaBatcher.MODE[args.s]
    args.m = BatcherThreading.FEED_MODE[args.m]
    args.c = KJoiner.MODE[args.c]
    args.M = KJoiner.MEMORY[args.M]

    if args.B is not None:
        assert os.path.isdir(args.B), f"batching output folder not found: '{args.B}'"
        assert 0 < len(
            os.listdir(args.B)
        ), f"the batching output folder is empty: '{args.B}'"

    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    if args.B is not None:
        logging.info(f"Loading previously generated batches from '{args.B}'...")
        batches = BatcherThreading.from_files(args.B, args.t, reSort=args.do_resort)
    else:
        batcher = FastaBatcher(size=args.b, threads=args.t)
        batcher.mode = args.s
        batcher.doReverseComplement = args.do_reverse
        batcher.do(args.input, args.k, args.m)
        batches = batcher.collection

    joiner = KJoinerThreading(args.c, args.M)
    joiner.threads = args.t

    joiner.batch_size = max(2, int(len(batches) / args.t))
    rlimits = resource.getrlimit(resource.RLIMIT_NOFILE)
    joiner.batch_size = min(joiner.batch_size, rlimits[1])
    resource.setrlimit(
        resource.RLIMIT_NOFILE, (max(rlimits[0], joiner.batch_size), rlimits[1])
    )

    joiner.join(batches, args.output)
    logging.info("That's all! :smiley:")
