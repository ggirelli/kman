"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import os

from kmermaid.batch import Batch, BatchAppendable


def test_Batch():
    b = Batch(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    if b.size != 5:
        raise AssertionError
    if b.remaining != 5:
        raise AssertionError
    if b.current_size != 0:
        raise AssertionError
    if b.is_written is not False:
        raise AssertionError
    if str != b.type:
        raise AssertionError
    if len(list(b.record_gen())) != 0:
        raise AssertionError
    if len(list(b.sorted())) != 0:
        raise AssertionError

    if not isinstance(b.check_record("test"), type(None)):
        raise AssertionError
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
    if b.current_size != 1:
        raise AssertionError
    if b.size - b.current_size != b.remaining:
        raise AssertionError

    b.add_all(["Second record", "Third record"])
    if b.current_size != 3:
        raise AssertionError
    if b.size - b.current_size != b.remaining:
        raise AssertionError

    b.add_all(["4th record"])

    if list(b.record_gen()) != list(b.to_write()):
        raise AssertionError

    b.write()
    if not os.path.isfile(b.tmp):
        raise AssertionError
    if not b.is_written:
        raise AssertionError
    if b.collection != [None]:
        raise AssertionError
    if b.current_size != 4:
        raise AssertionError
    if len(list(b.record_gen())) != 4:
        raise AssertionError

    b2 = b.from_file(b.tmp, str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    if b2.current_size != 4:
        raise AssertionError
    if b.tmp != b2.tmp:
        raise AssertionError
    if not b2.is_written:
        raise AssertionError

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    if rec1 != rec2:
        raise AssertionError

    b.unwrite()
    if b.current_size != 4:
        raise AssertionError
    if b.is_written:
        raise AssertionError

    b.add("5th record")
    if not b.is_full():
        raise AssertionError

    recList = [
        "4th record\n",
        "5th record",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    if recList != list(b.sorted()):
        raise AssertionError

    b.write()
    b.reset()
    if b.current_size != 0:
        raise AssertionError
    if b.size != 5:
        raise AssertionError
    if b.remaining != 5:
        raise AssertionError
    if b.is_written:
        raise AssertionError
    if os.path.isfile(b.tmp):
        raise AssertionError


def test_BatchAppendable():
    b = BatchAppendable(str, ".", 5)
    b.isFasta = False
    b.fwrite = "__str__"
    b.fread = "__str__"
    b.keyAttr = "__str__"

    if b.size != 5:
        raise AssertionError
    if b.remaining != 5:
        raise AssertionError
    if b.current_size != 0:
        raise AssertionError
    if b.is_written is not True:
        raise AssertionError
    if str != b.type:
        raise AssertionError
    if len(list(b.record_gen())) != 0:
        raise AssertionError
    if len(list(b.sorted())) != 0:
        raise AssertionError

    if not isinstance(b.check_record("test"), type(None)):
        raise AssertionError
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
    if b.current_size != 1:
        raise AssertionError
    if b.size - b.current_size != b.remaining:
        raise AssertionError

    b.add_all(["Second record", "Third record"])
    if b.current_size != 3:
        raise AssertionError
    if b.size - b.current_size != b.remaining:
        raise AssertionError

    b.add_all(["4th record"])

    if list(b.record_gen()) != list(b.to_write()):
        raise AssertionError

    b2 = b.from_file(b.tmp, str, False)
    b2.isFasta = False
    b2.fwrite = "__str__"
    b2.fread = "__str__"
    b2.keyAttr = "__str__"
    if b2.current_size != 4:
        raise AssertionError
    if b.tmp != b2.tmp:
        raise AssertionError
    if not b2.is_written:
        raise AssertionError

    rec1 = list(b.record_gen())
    rec2 = list(b2.record_gen())
    if rec1 != rec2:
        raise AssertionError

    b.unwrite()
    if b.current_size != 4:
        raise AssertionError
    if not b.is_written:
        raise AssertionError

    b.add("5th record")
    if not b.is_full():
        raise AssertionError

    recList = [
        "4th record\n",
        "5th record\n",
        "First record\n",
        "Second record\n",
        "Third record\n",
    ]
    if recList != list(b.sorted()):
        raise AssertionError(recList, list(b.sorted()))

    b.write()
    b.reset()
    if b.current_size != 0:
        raise AssertionError
    if b.size != 5:
        raise AssertionError
    if b.remaining != 5:
        raise AssertionError
    if not b.is_written:
        raise AssertionError
    if os.path.isfile(b.tmp):
        raise AssertionError
