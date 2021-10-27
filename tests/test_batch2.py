"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""
from __future__ import annotations

import inspect
import os
import pathlib
from typing import Iterator, Type

from kmermaid.batch2 import Batch, Batchable, BatchAppendable
from kmermaid.parsers import ParserBase


class BatchableTest(Batchable):
    _s: str

    def __init__(self, s: str):
        self._s = s

    @property
    def key(self) -> str:
        return self._s

    @staticmethod
    def from_raw(raw: str, /, **kwargs) -> BatchableTest:
        return BatchableTest(raw)

    def to_str(self) -> str:
        return self._s

    def supports_parser(self, parser: Type[ParserBase]) -> bool:
        return inspect.signature(parser.parse).return_annotation == Iterator[str]

    def __lt__(self, other: BatchableTest) -> bool:
        return self.to_str() < other.to_str()


class ParserTest(ParserBase):
    def __init__(self, path: pathlib.Path):
        self._FH = path.open("r+")

    def parse(self) -> Iterator[str]:
        yield from self._FH

    @staticmethod
    def parse_file(path: pathlib.Path) -> Iterator[str]:
        yield from ParserTest(path).parse()


def test_Batch():
    b = Batch(BatchableTest, ParserTest, 5)

    assert b.size == 5
    assert b.remaining == 5
    assert b.current_size == 0
    assert b.is_written() is False
    assert BatchableTest == b.record_type
    assert not list(b.record_gen())
    assert not list(b.record_gen(do_sort=True))
    assert isinstance(b.check_record(BatchableTest("test")), type(None))

    b.add(BatchableTest("First record"))
    assert b.current_size == 1
    assert b.size - b.current_size == b.remaining

    b.add_all([BatchableTest("Second record"), BatchableTest("Third record")])
    assert b.current_size == 3
    assert b.size - b.current_size == b.remaining

    b.add_all([BatchableTest("4th record")])

    assert ["First record", "Second record", "Third record", "4th record"] == list(
        b.to_str()
    )

    b.write()
    assert os.path.isfile(b.temp)
    assert b.is_written()
    assert b.collection == [None] * b.size
    assert b.current_size == 4
    assert len(list(b.record_gen())) == 4

    b2 = b.from_file(pathlib.Path(b.temp), BatchableTest, ParserTest, False)
    assert b2.current_size == 4
    assert b.temp == b2.temp
    assert b2.is_written()

    assert list(b.to_str()) == list(b2.to_str())

    b.unwrite()
    assert b.current_size == 4
    assert not b.is_written()

    b.add(BatchableTest("5th record"))
    assert b.is_full()

    recList = [
        "4th record\n",
        "5th record",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    assert recList == list(b.to_str(do_sort=True))

    b.write()
    b.reset()
    assert b.current_size == 0
    assert b.size == 5
    assert b.remaining == 5
    assert not b.is_written()
    assert not os.path.isfile(b.temp)


def test_BatchAppendable():
    b = BatchAppendable(BatchableTest, ParserTest, 5)

    print((b.temp, b.temp.exists()))

    assert b.size == 5
    assert b.remaining == 5
    assert b.current_size == 0
    assert b.is_written() is True
    assert BatchableTest == b.record_type
    assert not list(b.record_gen())
    assert not list(b.record_gen(do_sort=True))

    b.check_record(BatchableTest("test"))
    b.add(BatchableTest("First record"))
    assert b.current_size == 1
    assert b.size - b.current_size == b.remaining

    b.add_all(
        (x for x in [BatchableTest("Second record"), BatchableTest("Third record")])
    )
    assert b.current_size == 3
    assert b.size - b.current_size == b.remaining

    b.add_all((x for x in [BatchableTest("4th record")]))

    assert [
        "First record\n",
        "Second record\n",
        "Third record\n",
        "4th record\n",
    ] == list(b.to_str())

    b2 = b.from_file(pathlib.Path(b.temp), BatchableTest, ParserTest)
    assert b2.current_size == 4
    assert b.temp == b2.temp
    if not b2.is_written():
        raise AssertionError

    assert list(b.to_str()) == list(b2.to_str())

    b.unwrite()
    assert b.current_size == 4
    assert b.is_written()

    b.add(BatchableTest("5th record"))
    assert b.is_full()

    recList = [
        "4th record\n",
        "5th record\n",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    assert recList == list(b.to_str(do_sort=True)), (
        recList,
        list(b.record_gen(do_sort=True)),
    )

    b.write()
    b.reset()
    assert b.current_size == 0
    assert b.size == 5
    assert b.remaining == 5
    assert b.is_written
    assert not os.path.isfile(b.temp)
