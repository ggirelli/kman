"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: batch systems
"""

import copy
import os
import pathlib
import tempfile
import time
from typing import Any, Callable, Iterator, List, Optional, Type, Union

from kmermaid.parsers import ParserBase


class Batch:
    """Container for batched records.

    Records in the Batch are accessible through the record_gen and sorted
    methods, which source either from memory or from written files. A Batch
    cannot be resized. After full size is reached, a new Batch should be created

    :param _record_type: class type of records to be batched
    :type _record_type: Type
    :param _collection: list of batched records
    :type _collection: List[_record_type]
    :param _current_size: records currently in the batch instance, defaults to 0
    :type _current_size: int
    :param _size: maximum records allowed in the batch instance
    :type _size: int
    :param _temp: path to temporary folder
    :type _temp: pathlib.Path
    :param _written: if the batch was written to disk, defaults to False
    :type _written: bool
    :param parser: for file parsing
    :type parser: ParserBase
    :param ftorec: function to convert parser output to record
    :type ftorec: Callable[[_parser_output_type], _record_type]
    :param ftostr: function to make record writable
    :type ftostr: Callable[[_record_type], str]
    """

    _record_type: Union[Any]
    _parser_output_type: Union[Any]
    _collection: List[Optional["_record_type"]] = []

    _current_size: int = 0
    _size: int

    _temp: pathlib.Path = pathlib.Path(
        tempfile.NamedTemporaryFile(
            mode="w+",
            prefix=str(hash(time.time())),
            suffix=".fa",
        ).name
    )
    _written: bool = False

    _parser: ParserBase
    _ftorec: Callable[["_parser_output_type"], "_record_type"]
    _ftostr: Callable[["_record_type"], str]

    def __init__(
        self,
        record_type: Type,
        parser: ParserBase,
        ftorec: Callable[["_parser_output_type"], "_record_type"],
        ftostr: Callable[["_record_type"], str],
        size: int = 1,
    ):
        """Initialize Batch

        :param record_type: record type
        :type record_type: Type
        :param parser: for file parsing
        :type parser: ParserBase
        :param ftorec: function to convert parser output to record
        :type ftorec: Callable[[_parser_output_type], _record_type]
        :param ftostr: function to make record writable
        :type ftostr: Callable[[_record_type], str]
        :param size: records per batch, defaults to 1
        :type size: int
        :raises AssertionError: if size is lower than 1
        """
        super().__init__()
        if size < 1:
            raise AssertionError
        self._record_type = record_type
        self._parser = parser
        self._parser_output_type = parser.OUTPUT_TYPE
        self._ftorec = ftorec
        self._ftostr = ftostr
        self._init_collection(size)

    def _init_collection(self, size: int) -> None:
        """Initialize batch record collection.

        :param size: records per batch
        :type size: int
        """
        self._collection = [None] * self._size
        self._size = size

    @property
    def record_type(self):
        return self._record_type

    @property
    def collection(self) -> List["_record_type"]:
        return self._collection.copy()

    @property
    def current_size(self) -> int:
        return self._current_size

    @property
    def size(self) -> int:
        return self._size

    @property
    def remaining(self) -> int:
        return self._size - self._current_size

    @property
    def temp(self) -> pathlib.Path:
        return self._temp

    @property
    def parser(self) -> ParserBase:
        return copy.copy(self._parser)

    def is_written(self) -> bool:
        """Whether the batch was written to disk.

        :return: written to disk.
        :rtype: bool
        """
        return self._written

    def is_full(self) -> bool:
        """Whether the batch collection is full.

        :return: if filled
        :rtype: bool
        """
        return self.remaining == 0

    def info(self) -> str:
        """Print batch information in a readable format.

        :return: batch information
        :rtype: str
        """
        return f"""{self._temp}
        Record type: {self._record_type}
        Size: {self.current_size}/{self.size} ({self.remaining} remaining)
        Written: {self.is_written()}"""

    def record_gen(
        self,
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> Iterator["_record_type"]:
        """Generator of records.

        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :yield: record
        :rtype: Any
        :raises KeyError: if sort_by is not a valid key
        """
        if not self.is_written():
            records = (record for record in self._collection if record is not None)
        else:
            records = (
                self._ftorec(record) for record in self._parser.parse_file(self.temp)
            )
        if do_sort:
            if sort_by is not None:
                if sort_by not in dir(self._record_type):
                    raise KeyError(
                        f"attribute '{sort_by}' not found in {self._record_type}"
                    )
                yield from sorted(records, key=lambda x: getattr(x, str(sort_by)))
            else:
                yield from sorted(records)
        else:
            yield from records

    def check_record(self, record: "_record_type") -> None:
        """Check that record type matches the Batch.

        :param record
        :type record: _record_type
        :raises TypeError: if record type is not compatible
        """
        if not isinstance(record, self._record_type):
            raise TypeError(
                f"record must be '{self._record_type}', not '{type(record)}'."
            )

    def add(self, record: "_record_type"):
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record: rercord to be added
        :type record: _record_type
        :raises AssertionError: if batch is full or locally stored
        """
        if self.is_full():
            raise AssertionError("this batch is full.")
        if self.is_written:
            raise AssertionError("this batch has been stored locally.")
        self.check_record(record)
        self._collection[self._current_size] = record
        self._current_size += 1

    def add_all(self, record_generator: Iterator["_record_type"]) -> None:
        """Add all records from a generator to the current Batch.

        :param record_generator: record generator
        :type record_generator: Iterator[_record_type]
        """
        for record in record_generator:
            self.add(record)

    def to_str(
        self,
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> Iterator[str]:
        """Generator of writable records.

        Generator function that converts batch records into writable ones.

        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :yield: record as writable string
        :rtype: Iterator[str]
        """
        yield from (
            self._ftostr(record) for record in self.record_gen(do_sort, sort_by)
        )

    def write(
        self, do_sort: bool = False, sort_by: Optional[str] = None, force: bool = False
    ) -> None:
        """Write batch to disk.

        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :param force: overwrite, defaults to False
        :type force: bool
        """
        if not self.is_written() or force:
            with open(self._temp, "w+") as TH:
                TH.write(
                    "".join(
                        [
                            f"{x}\n" if not x.endswith("\n") else x
                            for x in self.to_str(do_sort, sort_by)
                        ]
                    )
                )
            self._collection = [None]
            self._written = True

    @staticmethod
    def from_file(
        path: pathlib.Path,
        record_type: Type,
        parser: ParserBase,
        ftorec: Callable[["_parser_output_type"], "_record_type"],
        ftostr: Callable[["_record_type"], str],
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> "Batch":
        """Generate a batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch.
        :type path: pathlib.Path
        :param record_type: record type
        :type record_type: Type
        :param parser: for file parsing
        :type parser: ParserBase
        :param ftorec: function to convert parser output to record
        :type ftorec: Callable[[_parser_output_type], _record_type]
        :param ftostr: function to make record writable
        :type ftostr: Callable[[_record_type], str]
        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :return: new batch
        :rtype: Batch
        """
        size = max(2, len(list(parser.parse_file(path))))
        batch = Batch(record_type, parser, ftorec, ftostr, size)
        batch._temp = path
        batch._current_size = size
        batch._written = True

        if do_sort:
            batch.write(do_sort, sort_by, force=True)

        return batch

    def reset(self):
        """Reset the batch.

        Empties current records collection and any written file.
        """
        if self.is_written():
            os.remove(self.tmp)
        self._written = False
        self._current_size = 0
        self._collection = [None] * self._size

    def unwrite(self):
        """Unwrites the batch from storage.

        Reads records from stored file to memory, if the batch is not full.
        """
        if not self.is_full() and self.is_written():
            self._collection = [None] * self.size
            self._collection[: self.current_size] = list(self.record_gen())
            self._written = False
            os.remove(self.tmp)


class BatchAppendable(Batch):
    """Batch container.

    Records in the Batch are accessible through the record_gen and sorted
    methods, which source either from memory or from written files. A Batch
    cannot be resized. After full size is reached, a new Batch should be created
    """

    _record_type: Union[Any]
    _parser_output_type: Union[Any]
    _parser: ParserBase

    _collection: List[Optional["_record_type"]] = [None]
    _written: bool = True

    def _init_collection(self, size: int) -> None:
        """Initialize batch record collection.

        :param size: records per batch
        :type size: int
        """
        self._collection = [None]
        self._size = size

    @property
    def info(self) -> str:
        """Prints Batch information in a readable format.

        :return: batch details
        :rtype: str
        """
        return f"{super().info()}\nAppending: True\n"

    def record_gen(
        self,
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> Iterator["_record_type"]:
        """Generator of records.

        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :yield: record
        :rtype: _record_type
        :raises KeyError: if sort_by is not a valid key
        """
        records = (
            self._ftorec(record) for record in self._parser.parse_file(self.temp)
        )
        if do_sort:
            if sort_by is None:
                if sort_by not in dir(self._record_type):
                    raise KeyError(
                        f"attribute '{sort_by}' not found in {self._record_type}"
                    )
                yield from sorted(records, key=lambda x: getattr(x, str(sort_by)))
            else:
                yield from sorted(records)
        else:
            yield from records

    def add(self, record: "_record_type") -> None:
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record: to add
        :type record: _record_type
        :raises AssertionError: if batch is full
        """
        if self.is_full():
            raise AssertionError("this batch is full")
        self.check_record(record)
        with open(self._temp, "a+") as OH:
            output = self._ftostr(record)
            output = output if output.endswith("\n") else f"{output}\n"
            OH.write(output)
        self._current_size += 1

    def write(
        self, do_sort: bool = False, sort_by: Optional[str] = None, force: bool = False
    ) -> None:
        """Write batch to disk, only if re-sorting is needed.


        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :param force: overwrite
        :type force: bool
        """
        if force and do_sort:
            with open(self.temp, "w+") as TH:
                TH.write(
                    "".join(
                        [
                            f"{x}\n" if not x.endswith("\n") else x
                            for x in self.to_str(do_sort, sort_by)
                        ]
                    )
                )

    @staticmethod
    def from_file(
        path: pathlib.Path,
        record_type: Type,
        parser: ParserBase,
        ftorec: Callable[["_parser_output_type"], "_record_type"],
        ftostr: Callable[["_record_type"], str],
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> "Batch":
        """Generate a batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch.
        :type path: pathlib.Path
        :param record_type: record type
        :type record_type: Type
        :param parser: for file parsing
        :type parser: ParserBase
        :param ftorec: function to convert parser output to record
        :type ftorec: Callable[[_parser_output_type], _record_type]
        :param ftostr: function to make record writable
        :type ftostr: Callable[[_record_type], str]
        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :return: new batch
        :rtype: Batch
        """
        size = max(2, len(list(parser.parse_file(path))))
        batch = BatchAppendable(record_type, parser, ftorec, ftostr, size)
        batch._temp = path
        batch._current_size = size

        if do_sort:
            batch.write(do_sort, sort_by, force=True)

        return batch

    def reset(self):
        """Reset the batch.

        Empties current records collection and any written file.
        """
        os.remove(self._temp)
        self._current_size = 0

    def unwrite(self):
        """Do nothing for Appendable batches."""
