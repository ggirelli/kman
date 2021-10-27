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
from abc import abstractmethod
from typing import Any, Iterator, List, Optional, Type

from kmermaid.parsers import ParserBase


class Sortable:
    """Class that can be sorted."""

    @abstractmethod
    def __lt__(self, other: Any) -> bool:
        """Check if current instance is lower than another instance.

        :param other: other instance
        :type other: Any
        :return: if current instance is lower
        :rtype: bool
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError


class Batchable(Sortable):
    """Class that can be collected in a Batch."""

    @staticmethod
    @abstractmethod
    def from_raw(raw: Any, /, **kwargs) -> "Batchable":
        """Raw form to Batchable object conversion.

        :param raw: to be converted to Batchable object
        :type raw: Any
        :param **kwargs: additional keyword arguments
        :return: a Batchable object
        :rtype: Batchable
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @abstractmethod
    def to_str(self) -> str:
        """Batchable object to String conversion.

        :return: String form of self
        :rtype: str
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError

    @abstractmethod
    def supports_parser(self, parser: Type[ParserBase]) -> bool:
        """Whether a parser is supported.

        :param parser: parser to check
        :type parser: Type[ParserBase]
        :return: if parser is supported
        :rtype: bool
        :raises NotImplementedError: abstract method
        """
        raise NotImplementedError


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
    """

    _record_type: Type[Batchable]
    _collection: List[Batchable] = []

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

    def __init__(
        self,
        record_type: Type[Batchable],
        parser: ParserBase,
        size: int = 1,
    ):
        """Initialize Batch

        :param record_type: record type
        :type record_type: Type[Batchable]
        :param parser: for file parsing
        :type parser: ParserBase
        :param size: records per batch, defaults to 1
        :type size: int
        :raises AssertionError: if size is lower than 1
        """
        super().__init__()
        if size < 1:
            raise AssertionError("a batch must allocate more than one record")
        self._record_type = record_type
        self._parser = parser
        self._size = size

    @property
    def record_type(self):
        return self._record_type

    @property
    def collection(self) -> List[Batchable]:
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
    ) -> Iterator[Batchable]:
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
                self._record_type.from_raw(record)
                for record in self._parser.parse_file(self.temp)
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

    def check_record(self, record: Batchable) -> None:
        """Check that record type matches the Batch.

        :param record
        :type record: Batchable
        :raises TypeError: if record type is not compatible
        """
        if not isinstance(record, self._record_type):
            raise TypeError(
                f"record must be '{self._record_type}', not '{type(record)}'."
            )

    def add(self, record: Batchable):
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record: rercord to be added
        :type record: Batchable
        :raises AssertionError: if batch is full or locally stored
        """
        if self.is_full():
            raise AssertionError("this batch is full.")
        if self.is_written():
            raise AssertionError("this batch has been stored locally.")
        self.check_record(record)
        self._collection[self._current_size] = record
        self._current_size += 1

    def add_all(self, record_generator: Iterator[Batchable]) -> None:
        """Add all records from a generator to the current Batch.

        :param record_generator: record generator
        :type record_generator: Iterator[Batchable]
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
        yield from (record.to_str() for record in self.record_gen(do_sort, sort_by))

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
            self._collection = []
            self._written = True

    @staticmethod
    def from_file(
        path: pathlib.Path,
        record_type: Type[Batchable],
        parser: ParserBase,
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> "Batch":
        """Generate a batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch.
        :type path: pathlib.Path
        :param record_type: record type
        :type record_type: Type[Batchable]
        :param parser: for file parsing
        :type parser: ParserBase
        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :return: new batch
        :rtype: Batch
        """
        size = max(2, len(list(parser.parse_file(path))))
        batch = Batch(record_type, parser, size)
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

    _record_type: Type[Batchable]
    _parser: ParserBase

    _collection: List[Batchable] = []
    _written: bool = True

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
    ) -> Iterator[Batchable]:
        """Generator of records.

        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :yield: record
        :rtype: Batchable
        :raises KeyError: if sort_by is not a valid key
        """
        records = (
            self._record_type.from_raw(record)
            for record in self._parser.parse_file(self.temp)
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

    def add(self, record: Batchable) -> None:
        """Add a record to the current batch.

        Does not work if the batch is full. Also, the record type must match the
        batch type. Automatically selects whether to append in memory or to the
        temporary file (if self.is_appending).

        :param record: to add
        :type record: Batchable
        :raises AssertionError: if batch is full
        """
        if self.is_full():
            raise AssertionError("this batch is full")
        self.check_record(record)
        with open(self._temp, "a+") as OH:
            output = record.to_str()
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
        record_type: Type[Batchable],
        parser: ParserBase,
        do_sort: bool = False,
        sort_by: Optional[str] = None,
    ) -> "Batch":
        """Generate a batch from a file.

        Used to link an existing file to a Batch instance.

        :param path: path to previously generated batch.
        :type path: pathlib.Path
        :param record_type: record type
        :type record_type: Type[Batchable]
        :param parser: for file parsing
        :type parser: ParserBase
        :param do_sort: sort records
        :type do_sort: bool
        :param sort_by: key for sorting
        :type sort_by: Optional[str]
        :return: new batch
        :rtype: Batch
        """
        size = max(2, len(list(parser.parse_file(path))))
        batch = BatchAppendable(record_type, parser, size)
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
