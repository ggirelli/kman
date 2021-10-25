"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
"""

import gzip
from kman.batch import Batch
import io
import os
import shutil
import tempfile
from tqdm import tqdm
from typing import List, Tuple


class SmartFastaParser(object):
    """Fasta parser with minimally open buffer.

    Extends the SimpleFastaParser function to avoid 'too many open files' issue.
    The input buffer is opened only when needed, and the byte-based position in
    the Fasta file is kept in memory. This continuous opening/closing adds a bit
    of an overhead but is useful when many Fasta files must be parsed at the
    same time.

    Variables:
            __compressed {bool} -- whether the Fasta is compressed.
            __pos {number} -- byte-based location in the Fasta file.
    """

    __compressed = False
    __pos = 0

    def __init__(self, FH):
        super(SmartFastaParser, self).__init__()
        if str == type(FH):
            if FH.endswith(".gz"):
                self.__compressed = True
                self.__FH = gzip.open(FH, "rt")
            else:
                self.__FH = open(FH, "r+")
        elif io.TextIOWrapper == type(FH):
            self.__FH = FH
            if self.__FH.name.endswith(".gz"):
                self.__compressed = True
        else:
            assert False, "type error."

    def __reopen(self):
        """Re-open the buffer and seek the last recorded position."""
        if self.__FH.closed:
            if self.__compressed:
                self.__FH = gzip.open(self.__FH.name, "rt")
            else:
                self.__FH = open(self.__FH.name, "r+")
        self.__FH.seek(self.__pos)

    def __skip_blank_and_comments(self) -> Tuple[str, bool]:
        # Skip any text before the first record (e.g. blank lines, comments)
        while True:
            line = self.__FH.readline()
            self.__pos = self.__FH.tell()
            if line == "":
                return ("", False)  # Premature end of file, or just empty?
            if line[0] == ">":
                break
        return (line, True)

    def __parse_sequence(self) -> List[str]:
        lines = []
        line = self.__FH.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            else:
                self.__pos = self.__FH.tell()
            lines.append(line.rstrip())
            line = self.__FH.readline()
        return lines

    def parse(self):
        """Iterate over Fasta records as string tuples.

        For each record a tuple of two strings is returned, the FASTA title
        line (without the leading '>' character), and the sequence (with any
        whitespace removed). The title line is not divided up into an
        identifier (the first word) and comment or description.

        Additionally, keep the Fasta handler open only when strictly necessary.
        """
        self.__reopen()
        line, parsable = self.__skip_blank_and_comments()
        assert parsable, "premature end of file or empty file"

        while True:
            self.__reopen()
            if line is None:
                line = self.__FH.readline()

            if not line:
                return

            if line[0] != ">":
                raise ValueError(
                    "Records in Fasta files should start with '>' character"
                )
            title = line[1:].rstrip()
            seq_lines = self.__parse_sequence()

            self.__FH.close()
            yield title, "".join(seq_lines).replace(" ", "").replace("\r", "")

            line = None

        assert False, "Should not reach this line"

    @staticmethod
    def parse_file(path):
        return SmartFastaParser(path).parse()


def set_tempdir(path: str, create: bool = True) -> None:
    if not os.path.isdir(path) and create:
        os.makedirs(path, exist_ok=True)
    else:
        raise AssertionError(f"folder not found: {path}")
    tempfile.tempdir = path


def input_file_exists(path: str) -> None:
    if not os.path.isfile(path):
        raise AssertionError(f"input file not found: {path}")


def copy_batches(
    batches: List[Batch], output_path: str, compress: bool = False
) -> None:
    """Copy generated batches to output folder.

    Args:
        batches (List[Batch]): generated batches.
        output_path (str): path to output folder.
        compress (bool, optional): gzip batches. Defaults to False.
    """
    batch_list = tqdm(
        (batch for batch in batches if os.path.isfile(batch.tmp)),
        total=len(batches),
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
