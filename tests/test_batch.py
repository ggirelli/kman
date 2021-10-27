"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import os
import pathlib

from kmermaid.batch import Batch, BatchAppendable


def test_Batch():
    b = Batch(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    assert b.size == 5
    assert b.remaining == 5
    assert b.current_size == 0
    assert b.is_written is False
    assert str == b.type
    assert not list(b.record_gen())
    assert not list(b.sorted())
    assert isinstance(b.check_record("test"), type(None))

    try:
        b.check_record(1)
    except AssertionError:
        pass
    else:
        raise AssertionError("record type must be tested")

    try:
        b.add(1)
    except AssertionError:
        pass
    else:
        raise AssertionError("record type must be tested when adding it")

    b.add("First record")
    assert b.current_size == 1
    assert b.size - b.current_size == b.remaining

    b.add_all(["Second record", "Third record"])
    assert b.current_size == 3
    assert b.size - b.current_size == b.remaining

    b.add_all(["4th record"])

    assert list(b.record_gen()) == list(b.to_write())

    b.write()
    assert os.path.isfile(b.tmp)
    assert b.is_written
    assert b.collection == [None]
    assert b.current_size == 4
    assert len(list(b.record_gen())) == 4

    b2 = b.from_file(pathlib.Path(b.tmp), str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    assert b2.current_size == 4
    assert b.tmp == b2.tmp
    assert b2.is_written

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    assert rec1 == rec2

    b.unwrite()
    assert b.current_size == 4
    assert not b.is_written

    b.add("5th record")
    assert b.is_full()

    recList = [
        "4th record\n",
        "5th record",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    assert recList == list(b.sorted())

    b.write()
    b.reset()
    assert b.current_size == 0
    assert b.size == 5
    assert b.remaining == 5
    assert not b.is_written
    assert not os.path.isfile(b.tmp)


def test_BatchAppendable():
    b = BatchAppendable(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    assert b.size == 5
    assert b.remaining == 5
    assert b.current_size == 0
    assert b.is_written is True
    assert str == b.type
    assert not list(b.record_gen())
    assert not list(b.sorted())

    continue_test_BatchAppendable_2(b)


def continue_test_BatchAppendable_2(b: Batch):
    b.check_record("test")
    try:
        b.check_record(1)
    except AssertionError:
        pass
    else:
        raise AssertionError("record type must be tested")

    try:
        b.add(1)
    except AssertionError:
        pass
    else:
        raise AssertionError("record type must be tested when adding it")

    continue_test_BatchAppendable_3(b)


def continue_test_BatchAppendable_3(b: Batch):
    b.add("First record")
    assert b.current_size == 1
    assert b.size - b.current_size == b.remaining

    b.add_all(["Second record", "Third record"])
    assert b.current_size == 3
    assert b.size - b.current_size == b.remaining

    b.add_all(["4th record"])

    assert list(b.record_gen()) == list(b.to_write())

    b2 = b.from_file(pathlib.Path(b.tmp), str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    assert b2.current_size == 4
    assert b.tmp == b2.tmp
    if not b2.is_written:
        raise AssertionError

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    assert rec1 == rec2

    b.unwrite()
    assert b.current_size == 4
    assert b.is_written

    b.add("5th record")
    assert b.is_full()

    recList = [
        "4th record\n",
        "5th record\n",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    assert recList == list(b.sorted()), (recList, list(b.sorted()))

    b.write()
    b.reset()
    assert b.current_size == 0
    assert b.size == 5
    assert b.remaining == 5
    assert b.is_written
    assert not os.path.isfile(b.tmp)
