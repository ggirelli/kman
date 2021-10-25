"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import oligo_melting as om  # type: ignore

from kmermaid.seq import KMer, Sequence, SequenceCoords, SequenceCount


def test_SequenceCoords_start():
    try:
        SequenceCoords("chr1", -1, 1, SequenceCoords.STRAND.MINUS)
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("start must be tested to be >=0")


def test_SequenceCoords_end():
    try:
        SequenceCoords("chr1", 0, -1, SequenceCoords.STRAND.MINUS)
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("end must be tested to be >=0")


def test_SequenceCoords_strand():
    try:
        SequenceCoords("chr1", 0, 1, "minus")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("strand must be tested to be from SequenceCoords.STRAND")


def test_SequenceCoords_2():
    ref = "chr1"
    start = 0
    end = 1e3
    strand = SequenceCoords.STRAND.PLUS
    sc = SequenceCoords(ref, start, end, strand)
    if ref != sc.ref:
        raise AssertionError
    if start != sc.start:
        raise AssertionError
    if end != sc.end:
        raise AssertionError
    if strand != sc.strand:
        raise AssertionError

    if strand != SequenceCoords.rev(SequenceCoords.STRAND.MINUS):
        raise AssertionError
    if SequenceCoords.STRAND.MINUS != SequenceCoords.rev(strand):
        raise AssertionError

    strRep = "%s:%d-%d:%s" % (ref, start, end, strand.label)
    if strRep != str(sc):
        raise AssertionError
    if sc != SequenceCoords.from_str(strRep):
        raise AssertionError


def test_KMer_start():
    try:
        KMer("chr1", -1, 10, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("start must be tested to be >=0")


def test_KMer_end():
    try:
        KMer("chr1", 10, -1, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("end must be tested to be >=0")


def test_KMer_strand():
    try:
        KMer("chr1", 0, 11, "ACGATCGATCG", strand="minus")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("strand must be tested to be from SequenceCoords.STRAND")


def test_KMer_NATYPE():
    try:
        KMer("chr1", 0, 11, "ACGATCGATCG", t="DNA")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")


def test_KMer_length():
    try:
        KMer("chr1", 0, 1, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("length must be tested for match")


def test_KMer():
    seq = "ACGATCGATCG"
    k = KMer("chr1", 0, len(seq), seq)
    sc = SequenceCoords("chr1", 0, 11, SequenceCoords.STRAND.PLUS)
    if sc != k.coords:
        raise AssertionError
    if str(sc) != k.header:
        raise AssertionError
    if seq != k.seq:
        raise AssertionError

    if k != KMer.from_fasta((k.header, k.seq)):
        raise AssertionError
    if ">%s\n%s\n" % (k.header, k.seq) != k.as_fasta():
        raise AssertionError
    if "%s\t%s" % (k.header, k.seq) != str(k):
        raise AssertionError

    if not k.is_ab_checked():
        raise AssertionError


def test_Sequence():
    try:
        Sequence("ACGATCGATCG", "DNA")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")

    ref = "stest"
    s = Sequence("ACGAT", om.NATYPES.DNA, ref)
    k4mer = [KMer(ref, 0, 4, "ACGA"), KMer(ref, 1, 5, "CGAT")]

    if k4mer != list(s.kmers(4)):
        raise AssertionError(k4mer, list(s.kmers(4)))
    if k4mer != list(s.kmerator(s.text, 4, s.natype, ref)):
        raise AssertionError

    s = Sequence("ACGATCGATCG", om.NATYPES.DNA)
    if s != Sequence("ACGATCGATCG", om.NATYPES.DNA):
        raise AssertionError
    if s == Sequence("ACGATCGATCG", om.NATYPES.RNA):
        raise AssertionError

    batches = [("ACGAT", 0), ("ATCGA", 3), ("GATCG", 6)]
    if batches != list(s.batches(3, 5)):
        raise AssertionError
    if batches != list(s.batcher(s.text, 3, 5)):
        raise AssertionError

    s = Sequence("ACGATCGATCG", om.NATYPES.DNA, "ref")
    klist = [
        [KMer("ref", 0, 0 + 4, "ACGA"), KMer("ref", 1, 1 + 4, "CGAT")],
        [KMer("ref", 2, 2 + 4, "GATC"), KMer("ref", 3, 3 + 4, "ATCG")],
        [KMer("ref", 4, 4 + 4, "TCGA"), KMer("ref", 5, 5 + 4, "CGAT")],
        [KMer("ref", 6, 6 + 4, "GATC"), KMer("ref", 7, 7 + 4, "ATCG")],
    ]
    if klist != [list(g) for g in s.kmers_batched(4, 5)]:
        raise AssertionError
    if klist != [
        list(g) for g in s.kmerator_batched(s.text, 4, s.natype, 5, s.name)
    ]:
        raise AssertionError

    mStrand = SequenceCoords.STRAND.MINUS
    klist = [
        [
            KMer("ref", 0, 0 + 4, "ACGA"),
            KMer("ref", 0, 0 + 4, "TCGT", strand=mStrand),
            KMer("ref", 1, 1 + 4, "CGAT"),
            KMer("ref", 1, 1 + 4, "ATCG", strand=mStrand),
        ],
        [
            KMer("ref", 2, 2 + 4, "GATC"),
            KMer("ref", 2, 2 + 4, "GATC", strand=mStrand),
            KMer("ref", 3, 3 + 4, "ATCG"),
            KMer("ref", 3, 3 + 4, "CGAT", strand=mStrand),
        ],
        [
            KMer("ref", 4, 4 + 4, "TCGA"),
            KMer("ref", 4, 4 + 4, "TCGA", strand=mStrand),
            KMer("ref", 5, 5 + 4, "CGAT"),
            KMer("ref", 5, 5 + 4, "ATCG", strand=mStrand),
        ],
        [
            KMer("ref", 6, 6 + 4, "GATC"),
            KMer("ref", 6, 6 + 4, "GATC", strand=mStrand),
            KMer("ref", 7, 7 + 4, "ATCG"),
            KMer("ref", 7, 7 + 4, "CGAT", strand=mStrand),
        ],
    ]
    if klist != [
        list(g) for g in s.kmerator_batched(s.text, 4, s.natype, 5, s.name, True)
    ]:
        raise AssertionError


def test_SequenceCount():
    try:
        SequenceCount("ACGATCGATCG", [1, 2, 3], "DNA")
    except AssertionError:
        pass
    else:
        if not False:
            raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")

    sc1 = str(SequenceCoords("chr1", 0, 1e3, SequenceCoords.STRAND.PLUS))
    sc2 = str(SequenceCoords("chr1", 1e3, 2e3, SequenceCoords.STRAND.PLUS))
    sco = SequenceCount("ACGATCGATCG", [sc1, sc2], om.NATYPES.DNA)
    if sco.header != [sc1, sc2]:
        raise AssertionError
    if sco.seq != sco.text:
        raise AssertionError

    strRepr = "%s\t%s" % (sco.seq, " ".join(sco.header))
    if str(sco) != strRepr:
        raise AssertionError
    if sco != sco.from_text(strRepr):
        raise AssertionError
    if str(sco) + "\n" != sco.as_text():
        raise AssertionError
