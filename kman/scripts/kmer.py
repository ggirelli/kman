"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from kman import __version__
from kman.const import CONTEXT_SETTINGS
from kman.scripts import kmer_batch, kmer_count, kmer_uniq
import sys


@click.group(
    name="radiant",
    context_settings=CONTEXT_SETTINGS,
    help=f"""\b
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/kman
Code:       http://github.com/ggirelli/kman

K-mer management tools.
""",
)
@click.version_option(__version__)
def main():
    """This is just an entry point."""
    pass


main.add_command(kmer_batch.run)
main.add_command(kmer_count.run)
main.add_command(kmer_uniq.run)
