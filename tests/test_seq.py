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
        raise AssertionError("start must be tested to be >=0")


def test_SequenceCoords_end():
    try:
        SequenceCoords("chr1", 0, -1, SequenceCoords.STRAND.MINUS)
    except AssertionError:
        pass
    else:
        raise AssertionError("end must be tested to be >=0")


def test_SequenceCoords_strand():
    try:
        SequenceCoords("chr1", 0, 1, "minus")
    except AssertionError:
        pass
    else:
        raise AssertionError("strand must be tested to be from SequenceCoords.STRAND")


def test_SequenceCoords_2():
    ref = "chr1"
    start = 0
    end = 1e3
    strand = SequenceCoords.STRAND.PLUS
    sc = SequenceCoords(ref, start, end, strand)
    assert ref == sc.ref
    assert start == sc.start
    assert end == sc.end
    assert strand == sc.strand

    assert strand == SequenceCoords.rev(SequenceCoords.STRAND.MINUS)
    assert SequenceCoords.STRAND.MINUS == SequenceCoords.rev(strand)

    strRep = "%s:%d-%d:%s" % (ref, start, end, strand.label)
    assert strRep == str(sc)
    assert sc == SequenceCoords.from_str(strRep)


def test_KMer_start():
    try:
        KMer("chr1", -1, 10, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        raise AssertionError("start must be tested to be >=0")


def test_KMer_end():
    try:
        KMer("chr1", 10, -1, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        raise AssertionError("end must be tested to be >=0")


def test_KMer_strand():
    try:
        KMer("chr1", 0, 11, "ACGATCGATCG", strand="minus")
    except AssertionError:
        pass
    else:
        raise AssertionError("strand must be tested to be from SequenceCoords.STRAND")


def test_KMer_NATYPE():
    try:
        KMer("chr1", 0, 11, "ACGATCGATCG", t="DNA")
    except AssertionError:
        pass
    else:
        raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")


def test_KMer_length():
    try:
        KMer("chr1", 0, 1, "ACGATCGATCG")
    except AssertionError:
        pass
    else:
        raise AssertionError("length must be tested for match")


def test_KMer():
    seq = "ACGATCGATCG"
    k = KMer("chr1", 0, len(seq), seq)
    sc = SequenceCoords("chr1", 0, 11, SequenceCoords.STRAND.PLUS)
    assert sc == k.coords
    assert str(sc) == k.header
    assert seq == k.seq

    assert k == KMer.from_fasta((k.header, k.seq))
    assert ">%s\n%s\n" % (k.header, k.seq) == k.as_fasta()
    assert "%s\t%s" % (k.header, k.seq) == str(k)

    assert k.is_ab_checked()


def test_Sequence():
    try:
        Sequence("ACGATCGATCG", "DNA")
    except AssertionError:
        pass
    else:
        raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")

    ref = "stest"
    s = Sequence("ACGAT", om.NATYPES.DNA, ref)
    k4mer = [KMer(ref, 0, 4, "ACGA"), KMer(ref, 1, 5, "CGAT")]

    assert k4mer == list(s.kmers(4)), (k4mer, list(s.kmers(4)))
    assert k4mer == list(s.kmerator(s.text, 4, s.natype, ref))

    s = Sequence("ACGATCGATCG", om.NATYPES.DNA)
    assert s == Sequence("ACGATCGATCG", om.NATYPES.DNA)
    assert s != Sequence("ACGATCGATCG", om.NATYPES.RNA)

    batches = [("ACGAT", 0), ("ATCGA", 3), ("GATCG", 6)]
    assert batches == list(s.batches(3, 5))
    assert batches == list(s.batcher(s.text, 3, 5))

    s = Sequence("ACGATCGATCG", om.NATYPES.DNA, "ref")
    klist = [
        [KMer("ref", 0, 0 + 4, "ACGA"), KMer("ref", 1, 1 + 4, "CGAT")],
        [KMer("ref", 2, 2 + 4, "GATC"), KMer("ref", 3, 3 + 4, "ATCG")],
        [KMer("ref", 4, 4 + 4, "TCGA"), KMer("ref", 5, 5 + 4, "CGAT")],
        [KMer("ref", 6, 6 + 4, "GATC"), KMer("ref", 7, 7 + 4, "ATCG")],
    ]
    assert klist == [list(g) for g in s.kmers_batched(4, 5)]
    assert klist == [
        list(g) for g in s.kmerator_batched(s.text, 4, s.natype, 5, s.name)
    ]

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
    assert klist == [
        list(g) for g in s.kmerator_batched(s.text, 4, s.natype, 5, s.name, True)
    ]


def test_SequenceCount():
    try:
        SequenceCount("ACGATCGATCG", [1, 2, 3], "DNA")
    except AssertionError:
        pass
    else:
        raise AssertionError("nucl. acid type must be tested to be from om.NATYPES")

    sc1 = str(SequenceCoords("chr1", 0, 1e3, SequenceCoords.STRAND.PLUS))
    sc2 = str(SequenceCoords("chr1", 1e3, 2e3, SequenceCoords.STRAND.PLUS))
    sco = SequenceCount("ACGATCGATCG", [sc1, sc2], om.NATYPES.DNA)
    assert sco.header == [sc1, sc2]
    assert sco.seq == sco.text

    strRepr = "%s\t%s" % (sco.seq, " ".join(sco.header))
    assert str(sco) == strRepr
    assert sco == sco.from_text(strRepr)
    assert str(sco) + "\n" == sco.as_text()
