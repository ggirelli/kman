"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: batch systems
"""

import gzip
import os
import tempfile
import time
from typing import IO, Any, Iterator, List, Type

from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore

from kmermaid.parsers import SmartFastaParser
from kmermaid.seq import KMer


class Batch(object):
    """Batch container.

    Records in the Batch are accessible through the record_gen and sorted
    methods, which source either from memory or from written files. A Batch
    cannot be resized. After full size is reached, a new Batch should be created

    Variables:
            _fread {str} -- name of type method to read a record
            _fwrite {str} -- name of type method to convert a record to string
            _keyAttr {str} -- name of type attribute to extract that provides the
                              key for sorting
            _written {bool} -- if the batch was written to file
            _i {number} -- current record location
            _tmp_dir {tempfile.TemporaryDirectory}
            _tmp {tempfile.TemporaryFile}
            isFasta {bool} -- whether the output should be in fasta format
            suffix {str} -- extension for the output temporary file
    """

    _fread = "from_file"
    _fwrite = "as_fasta"
    _keyAttr = "seq"

    _written = False
    _i: int = 0
    _tmp_dir = None
    _tmp = None
    isFasta = True
    suffix = ".fa"

    def __init__(self, t: Type, tmpDir: str, size: int = 1):
        """Initialize Batch

        :param t: record type
        :type t: Type
        :param tmpDir: path to temporary folder
        :type tmpDir: str
        :param size: records per batch, defaults to 1
        :type size: int
        :raises AssertionError: if size is lower than 1
        """
        super().__init__()
        if size < 1:
            raise AssertionError
        self.__size = int(size)
        self._remaining = self.__size
        self.__type = t
        self._tmp_dir = tmpDir
        self.__records = [None] * self.__size

    @property
    def is_written(self):
        return self._written

    @property
    def current_size(self) -> int:
        return self._i

    @property
    def size(self):
        return self.__size

    @property
    def remaining(self):
        return self._remaining

    @property
    def collection(self):
        if isinstance(self.__records, type(None)):
            return None
        return self.__records.copy()

    @property
    def type(self):
        return self.__type

    @property
    def tmp(self):
        if self._tmp is None:
            with tempfile.NamedTemporaryFile(
                mode="w+",
                dir=self._tmp_dir,
                prefix=str(hash(time.time())),
                suffix=self.suffix,
            ) as TH:
                self._tmp = TH.name
        return self._tmp

    @property
    def info(self) -> str:
        """Prints Batch information in a readable format.

        :return: batch information
        :rtype: str
        """
        info = "%s\ntype: %s\nsize: %d" % (self.tmp, self.type, self.size)
        info += "\ni: %d\nremaining: %d" % (self.current_size, self.remaining)
        info += "\nwritten: %r\n" % self.is_written
        return info

    @property
    def keyAttr(self):
        return self._keyAttr

    @keyAttr.setter
    def keyAttr(self, k):
        if not isinstance(k, str):
            raise AssertionError
        if not hasattr(self.type, k):
            raise AssertionError
        self._keyAttr = k

    @property
    def fread(self):
        return self._fread

    @fread.setter
    def fread(self, f):
        if not isinstance(f, str):
            raise AssertionError
        if not hasattr(self.type, f):
            raise AssertionError
        self._fread = f

    @property
    def fwrite(self):
        return self._fwrite

    @fwrite.setter
    def fwrite(self, f):
        if not isinstance(f, str):
            raise AssertionError
        if not hasattr(self.type, f):
            raise AssertionError
        self._fwrite = f

    def sorted(self, smart: bool = False) -> Any:
        """Generator of sorted record.

        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :return: generator
        :rtype: Any
        """
        if self.isFasta:
            return sorted(
                self.record_gen(smart), key=lambda x: getattr(x, self.keyAttr)
            )
        return sorted(self.record_gen(smart))

    def _record_gen_from_handle(self, TH: IO, smart=False):
        if self.isFasta:
            fasta_parser = SmartFastaParser.parse_file if smart else SimpleFastaParser
            for record in fasta_parser(TH):
                yield getattr(self.__type, self.fread)(record)
        else:
            for line in TH:
                yield getattr(self.__type, self.fread)(line)

    def _record_gen_from_file(self, smart: bool = False) -> Any:
        """Generator of records, read from file.

        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :yield: record
        :rtype: Any
        """
        if self.tmp.endswith(".gz"):
            TH = gzip.open(self.tmp, "rt")
        else:
            TH = open(self.tmp, "r+")
        yield from self._record_gen_from_handle(TH, smart)
        if not TH.closed:
            TH.close()

    def record_gen(self, smart: bool = False) -> Any:
        """Generator of records.

        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :yield: record
        :rtype: Any
        """
        if self.is_written:
            yield from self._record_gen_from_file(smart)
        else:
            for record in self.__records:
                if type(None) != type(record):
                    yield record

    def check_record(self, record: Any) -> None:
        """Check that record type matches the Batch.

        :param record
        :type record: Any
        :raises AssertionError: if record type is not compatible
        """
        if type(record) != self.type:
            raise AssertionError(f"record must be {self.type}, not {type(record)}.")

    def add(self, record: Any):
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record
        :type record: Any
        :raises AssertionError: if batch is full or locally stored
        """
        if self.is_full():
            raise AssertionError("this batch is full.")
        if self.is_written:
            raise AssertionError("this batch has been stored locally.")
        self.check_record(record)
        self.__records[self._i] = record
        self._i += 1
        self._remaining -= 1

    def add_all(self, recordGen):
        """Add all records from a generator to the current Batch.

        :param recordGen: record generator
        :type recordGen: [type]
        """
        for record in recordGen:
            self.add(record)

    def to_write(self, doSort: bool = False) -> List[Any]:
        """Generator of writeable records.

        Generator function that converts batch records into writeable ones.

        :param doSort: sort when writing, defaults to False
        :type doSort: bool
        :return: output status
        :rtype: List[Any]
        """
        if doSort:
            return [
                getattr(r, self.fwrite)()
                for r in self.sorted()
                if not type(None) == type(r)
            ]
        return [
            getattr(r, self.fwrite)()
            for r in self.record_gen()
            if type(None) != type(r)
        ]

    def write(self, doSort: bool = False, force: bool = False) -> None:
        """Write batch to disk.

        :param doSort: sort while writing, defaults to False
        :type doSort: bool
        :param force: overwrite, defaults to False
        :type force: bool
        """
        if not self.is_written or force:
            output = [
                x + "\n" if not x.endswith("\n") else x for x in self.to_write(doSort)
            ]
            with open(self.tmp, "w+") as TH:
                TH.write("".join(output))
            self.__records = [None]
            self._written = True

    @staticmethod
    def from_file(
        path: str,
        t: Type = KMer,
        isFasta: bool = True,
        smart: bool = False,
        reSort: bool = False,
    ) -> "Batch":
        """Generate a batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch.
        :type path: str
        :param t: record type, defaults to KMer
        :type t: Type
        :param isFasta: if the input is FASTA, defaults to True
        :type isFasta: bool
        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :param reSort: re-sort while reading, defaults to False
        :type reSort: bool
        :return: read batch
        :rtype: Batch
        """
        if isFasta:
            FH = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r+")
            if smart:
                size = max(2, sum(1 for record in SmartFastaParser(FH).parse()))
            else:
                size = max(2, sum(1 for record in SimpleFastaParser(FH)))
        else:
            FH = open(path, "r+")
            size = max(2, sum(1 for line in FH))

        batch = Batch(t, os.path.dirname(path), size)
        batch._tmp = FH.name
        batch._i = size
        batch._remaining = 0
        batch._written = True
        batch.isFasta = isFasta

        if reSort:
            batch.write(doSort=True, force=True)

        FH.close()
        return batch

    @staticmethod
    def from_batcher(
        batch_type: Type, size: int = 1, tmp: str = tempfile.gettempdir()
    ) -> "Batch":
        """Initialize a batch using a Batcher for default values.

        :param batch_type: type of items in the batch
        :type batch_type: Type
        :param size: records per batch, defaults to 1
        :type size: int
        :param tmp: path to temporary folder
        :type tmp: str
        :return: new empty batch
        :rtype: Batch
        :raises AssertionError: if size is 0 or negative.
        """
        if size < 1:
            raise AssertionError(f"size cannot be 0 or negative: {size}")
        return Batch(batch_type, tmp, size)

    def reset(self):
        """Reset the batch.

        Empties current records collection and any written file.
        """
        if self.is_written:
            os.remove(self.tmp)
        self._written = False
        self._i = 0
        self._remaining = self.size
        self.__records = [None] * self.__size

    def is_full(self) -> bool:
        """Whether the Batch collection is full.

        :return: filling status
        :rtype: bool
        """
        return self.remaining == 0

    def unwrite(self):
        """Unwrites the batch from storage.

        Reads records from stored file to memory, if the batch is not full.
        """
        if not self.is_full() and self.is_written:
            self.__records = [None] * self.size
            self.__records[: self.current_size] = list(self.record_gen())
            self._written = False
            os.remove(self.tmp)


class BatchAppendable(Batch):
    """Batch container.

    Records in the Batch are accessible through the record_gen and sorted
    methods, which source either from memory or from written files. A Batch
    cannot be resized. After full size is reached, a new Batch should be created
    """

    def __init__(self, t: Type, tmpDir: str, size: int = 1):
        """Initialize an appenable Batch.

        :param t: record type
        :type t: Type
        :param tmpDir: temporary folder
        :type tmpDir: str
        :param size: records per batch, defaults to 1
        :type size: int
        """
        super().__init__(t, tmpDir, size)
        self._written = True

    @property
    def info(self) -> str:
        """Prints Batch information in a readable format.

        :return: batch details
        :rtype: str
        """
        info = "%s\ntype: %s\nsize: %d" % (self.tmp, self.type, self.size)
        info += "\ni: %d\nremaining: %d" % (self.current_size, self.remaining)
        info += "\nappending: True\n"
        return info

    def record_gen(self, smart: bool = False):
        """Generator of records.

        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :yield: record
        :rtype: Iterator[Any]
        """
        if self.current_size != 0:
            yield from self._record_gen_from_file(smart)

    def add(self, record: Any) -> None:
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record: record
        :type record: Any
        :raises AssertionError: if batch is full
        """
        if self.is_full():
            raise AssertionError("this batch is full.")
        super().check_record(record)
        with open(self.tmp, "a+") as OH:
            output = getattr(record, self.fwrite)()
            if not output.endswith("\n"):
                output += "\n"
            OH.write(output)
        self._i += 1
        self._remaining -= 1

    def add_all(self, recordGen: Iterator[Any]) -> None:
        """Add recirds from a generator to the current Batch

        :param recordGen: record generator
        :type recordGen: Iterator[Any]
        """
        for record in recordGen:
            self.add(record)

    def write(self, doSort: bool = False, force: bool = False) -> None:
        """Write batch to disk, only if re-sorting is needed.

        :param doSort: sort while writing, defaults to False
        :type doSort: bool
        :param force: overwrite
        :type force: bool
        """
        if (not self.is_written or force) and doSort:
            with open(self.tmp, "w+") as TH:
                TH.write(
                    "".join(
                        [
                            x + "\n" if not x.endswith("\n") else x
                            for x in self.to_write(doSort)
                        ]
                    )
                )

    @staticmethod
    def from_file(
        path: str,
        t: Type = KMer,
        isFasta: bool = True,
        smart: bool = False,
        reSort: bool = False,
    ) -> "Batch":
        """Generate a Batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch
        :type path: str
        :param t: record type, defaults to KMer
        :type t: Type
        :param isFasta: is input FASYA, defaults to True
        :type isFasta: bool
        :param smart: use smarter IO, defaults to False
        :type smart: bool
        :param reSort: re-sort while reading, defaults to False
        :type reSort: bool
        :return: Batch
        :rtype: Batch
        """
        if isFasta:
            FH = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r+")
            if smart:
                size = max(2, sum(1 for record in SmartFastaParser(FH).parse()))
            else:
                size = max(2, sum(1 for record in SimpleFastaParser(FH)))
        else:
            FH = open(path, "r+")
            size = max(2, sum(1 for line in FH))

        batch = BatchAppendable(t, os.path.dirname(path), size)
        batch._tmp = FH.name
        batch._i = size
        batch._remaining = 0
        batch._written = True

        FH.close()
        return batch

    @staticmethod
    def from_batcher(
        batch_type: Type, size: int = 1, tmp: str = tempfile.gettempdir()
    ) -> "Batch":
        """Initialize a batch using a Batcher for default values.

        :param batch_type: type of items in the batch
        :type batch_type: Type
        :param size: records per batch, defaults to 1
        :type size: int
        :param tmp: path to temporary folder
        :type tmp: str
        :return: new empty batch
        :rtype: Batch
        :raises AssertionError: if size is 0 or negative.
        """
        if size < 1:
            raise AssertionError(f"size cannot be 0 or negative: {size}")
        return BatchAppendable(batch_type, tmp, size)

    def reset(self):
        """Reset the batch.

        Empties current records collection and any written file.
        """
        os.remove(self.tmp)
        self._i = 0
        self._remaining = self.size

    def unwrite(self, *args, **kwargs):
        pass
