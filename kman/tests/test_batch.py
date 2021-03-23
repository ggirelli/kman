"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from kman.batch import Batch, BatchAppendable
import os


def test_Batch():
    b = Batch(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    assert 5 == b.size
    assert 5 == b.remaining
    assert 0 == b.current_size
    assert b.is_written is False
    assert str == b.type
    assert 0 == len(list(b.record_gen()))
    assert 0 == len(list(b.sorted()))

    assert isinstance(b.check_record("test"), type(None))
    try:
        b.check_record(1)
    except AssertionError:
        pass
    else:
        assert False, "record type must be tested"

    try:
        b.add(1)
    except AssertionError:
        pass
    else:
        assert False, "record type must be tested when adding it"

    b.add("First record")
    assert 1 == b.current_size
    assert b.size - b.current_size == b.remaining

    b.add_all(["Second record", "Third record"])
    assert 3 == b.current_size
    assert b.size - b.current_size == b.remaining

    b.add_all(["4th record"])

    assert list(b.record_gen()) == list(b.to_write())

    b.write()
    assert os.path.isfile(b.tmp)
    assert b.is_written
    assert isinstance(b.collection, type(None))
    assert 4 == b.current_size
    assert 4 == len(list(b.record_gen()))

    b2 = b.from_file(b.tmp, str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    assert 4 == b2.current_size
    assert b.tmp == b2.tmp
    assert b2.is_written

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    assert rec1 == rec2

    b.unwrite()
    assert 4 == b.current_size
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
    assert 0 == b.current_size
    assert 5 == b.size
    assert 5 == b.remaining
    assert not b.is_written
    assert not os.path.isfile(b.tmp)


def test_BatchAppendable():
    b = BatchAppendable(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    assert 5 == b.size
    assert 5 == b.remaining
    assert 0 == b.current_size
    assert b.is_written is True
    assert str == b.type
    assert 0 == len(list(b.record_gen()))
    assert 0 == len(list(b.sorted()))

    assert isinstance(b.check_record("test"), type(None))
    try:
        b.check_record(1)
    except AssertionError:
        pass
    else:
        assert False, "record type must be tested"

    try:
        b.add(1)
    except AssertionError:
        pass
    else:
        assert False, "record type must be tested when adding it"

    b.add("First record")
    assert 1 == b.current_size
    assert b.size - b.current_size == b.remaining

    b.add_all(["Second record", "Third record"])
    assert 3 == b.current_size
    assert b.size - b.current_size == b.remaining

    b.add_all(["4th record"])

    assert list(b.record_gen()) == list(b.to_write())

    b2 = b.from_file(b.tmp, str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    assert 4 == b2.current_size
    assert b.tmp == b2.tmp
    assert b2.is_written

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    assert rec1 == rec2

    b.unwrite()
    assert 4 == b.current_size
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
    assert recList == list(b.sorted())

    b.write()
    b.reset()
    assert 0 == b.current_size
    assert 5 == b.size
    assert 5 == b.remaining
    assert b.is_written
    assert not os.path.isfile(b.tmp)
