"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
"""

import gzip
import pathlib
from abc import abstractmethod
from typing import IO, Any, Iterator, List, Tuple, Union

FASTA_SIMPLE_RECORD = Tuple[str, str]


class ParserBase:
    """Interface for file parser classes.

    :param _FH: buffer handle for input file
    :type _FH: IO
    """

    _FH: IO
    _OUTPUT_TYPE = Union[Any]
    OUTPUT_TYPE = Any

    @abstractmethod
    def __init__(self, path: pathlib.Path):
        """Initialize parser.

        :param path: path to file
        :type path: pathlib.Path
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @abstractmethod
    def parse(self) -> Iterator[_OUTPUT_TYPE]:
        """Parse a file.

        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def parse_file(path: pathlib.Path) -> Iterator[_OUTPUT_TYPE]:
        """Parse a file (static method).

        :param path: path to file
        :type path: pathlib.Path
        :yield: record iterator
        :rtype: Iterator[Any]
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError


class FastaParserBase(ParserBase):
    """Interface for FASTA file parser classes.

    :param _FH: buffer handle for input FASTA file
    :type _FH: IO
    """

    _OUTPUT_TYPE = Union[FASTA_SIMPLE_RECORD]
    OUTPUT_TYPE = _OUTPUT_TYPE
    _FH: IO

    @abstractmethod
    def __init__(self, path: pathlib.Path):
        """Initialize parser.

        :param path: path to FASTA file
        :type path: pathlib.Path
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @abstractmethod
    def parse(self) -> Iterator[_OUTPUT_TYPE]:
        """Parse a FASTA file.

        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def parse_file(path: pathlib.Path) -> Iterator[_OUTPUT_TYPE]:
        """Parse a FASTA file (static method).

        :param path: path to FASTA file
        :type path: pathlib.Path
        :yield: fasta record iterator
        :rtype: Iterator[FASTA_SIMPLE_RECORD]
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError


class SmartFastaParser(FastaParserBase):
    """Fasta parser with minimally open buffer.

    Extends the SimpleFastaParser function to avoid 'too many open files' issue. The
    input buffer is opened only when needed, and the byte-based position in the FASTA
    file is kept in memory. This continuous opening/closing adds a bit of an overhead
    but is useful when many FASTA files must be parsed at the same time.

    :param __compressed: whether the Fasta is compressed
    :type __compressed: bool
    :param __last_position: byte-based location in the Fasta file
    :type __last_position: int
    """

    __compressed: bool
    __last_position: int = 0

    def __init__(self, path: pathlib.Path):
        """Initialize parser.

        :param path: path to FASTA file
        :type path: pathlib.Path
        :raises AssertionError: if file does not exist
        """
        if not path.exists():
            raise AssertionError(f"file not found: {path}")
        if ".gz" in path.suffixes:
            self.__compressed = True
            self._FH = gzip.open(path, "rt")
        else:
            self.__compressed = False
            self._FH = path.open("r+")
        self._FH.close()

    def is_open(self) -> bool:
        """Check if the handle is still open.

        :return: if the handle is still open
        :rtype: bool
        """
        return not self._FH.closed

    def __readline(self) -> str:
        """Read next line. Empty if EOF is reached.

        :return: next line
        :rtype: str
        """
        if self.__compressed:
            return self._FH.readline().decode()
        return self._FH.readline()

    def __seek_last_position(self) -> bool:
        """Move buffer to last recorded position.

        :return: if the seek was successful
        :rtype: bool
        """
        if not self.is_open():
            return False
        self._FH.seek(self.__last_position)
        return True

    def __reopen(self) -> None:
        """Re-open the buffer and seek the last recorded position."""
        if not self.is_open():
            if self.__compressed:
                self._FH = gzip.open(self._FH.name, "rt")
            else:
                self._FH = open(self._FH.name, "r+")
            self.__seek_last_position()

    def __skip_blank_and_comments(self) -> Tuple[str, bool]:
        """Skip any text before the first record (e.g. blank lines, comments)

        :return: next content line, and False if EOF or empty line was reached
        :rtype: Tuple[str, bool]
        """
        while True:
            line = self.__readline()
            if self._FH.tell() == self.__last_position:
                return ("", False)  # EOF
            self.__last_position = self._FH.tell()
            if line.startswith(">"):
                return (line, True)

    def __parse_sequence(self) -> str:
        """Read all sequence lines for current record.

        :return: joined sequence lines
        :rtype: str
        """
        lines: List[str] = []
        while True:
            line = self.__readline()
            if not line:
                break
            if line.startswith(">"):
                break
            self.__last_position = self._FH.tell()
            lines.append(line.rstrip())
        return "".join(lines)

    def parse(self) -> Iterator[FastaParserBase._OUTPUT_TYPE]:
        """Iterate over FASTA records as string tuples.

        For each record, a tuple of two strings is returned: the FASTA header (without
        the leading '>' character), and the sequence (with any whitespace removed). The
        header line is not divided up into identifier (the first word) and comments or
        description. Additionally, keep the FASTA handler open only when strictly
        necessary (closed between record readings).

        :raises ValueError: when fasta does not pass format check
        :yield: (header, sequence)
        :rtype: Tuple[str, str]
        :raises AssertionError: if file is not parsable
        :raises ValueError: if a header line does not start with '>'
        :raises AssertionError: if the EOF is wrongly parsed
        """
        self.__reopen()

        line, parsable = self.__skip_blank_and_comments()
        if not parsable:
            raise AssertionError("premature end of file or no header found")
        header_found: bool = True

        while True:
            self.__reopen()

            line = self.__readline() if header_found else line
            if not line:
                return
            if not line.startswith(">"):
                raise ValueError("header lines must start with '>'")

            header = line[1:].rstrip()
            sequence = self.__parse_sequence().replace(" ", "").replace("\r", "")

            self._FH.close()
            header_found = False
            yield header, sequence

    @staticmethod
    def parse_file(path: pathlib.Path) -> Iterator[FastaParserBase._OUTPUT_TYPE]:
        """Parse a FASTA file in a smarter way.

        :param path: path to FASTA file
        :type path: str
        :yield: FASTA record iterator
        :rtype: Iterator[FASTA_SIMPLE_RECORD]
        """
        yield from SmartFastaParser(path).parse()
