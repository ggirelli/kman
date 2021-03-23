"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from kman.const import __version__
from kman.scripts import arguments as ap
from kman import scripts
import sys


def default_parser(*args) -> None:
    print("kmer -h for usage details.")
    sys.exit()


def main():
    parser = argparse.ArgumentParser(
        description=f"""
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/kman
Code:       http://github.com/ggirelli/kman

K-mer management tools.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.set_defaults(parse=default_parser)
    parser = ap.add_version_option(parser)

    subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )

    scripts.kmer_batch.init_parser(subparsers)
    scripts.kmer_count.init_parser(subparsers)
    scripts.kmer_uniq.init_parser(subparsers)

    args = parser.parse_args()
    args = args.parse(args)
    args.run(args)
