"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore

from kmermaid import __version__
from kmermaid.const import CONTEXT_SETTINGS
from kmermaid.scripts import kmer_batch, kmer_count, kmer_uniq


@click.group(
    name="radiant",
    context_settings=CONTEXT_SETTINGS,
    help=f"""\b
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/kmermaid
Code:       http://github.com/ggirelli/kmermaid

K-mer management tools.
""",
)
@click.version_option(__version__)
def main():
    """This is just an entry point."""


main.add_command(kmer_batch.run)
main.add_command(kmer_count.run)
main.add_command(kmer_uniq.run)
