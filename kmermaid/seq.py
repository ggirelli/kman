"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: sequence management systems
"""

import logging
import re
from enum import Enum, unique
from typing import Iterator, List, Tuple

import oligo_melting as om  # type: ignore


class SequenceCoords(object):
    """Reference genome window coordinates system.

    Variables:
            regexp {sre.SRE_PATTERN} -- regular expression to parse input strings
    """

    @unique
    class STRAND(Enum):
        """Reference genome strand.

        To be used only for strandedness. Coordinates always refer to the PLUS
        strand.

        Extends:
                Enum

        Variables:
                PLUS {number} -- positive strand
                MINUS {number} -- negative strand
        """

        PLUS = 0
        MINUS = 1

        @property
        def label(self):
            return "+-"[self.value]

    regexp = re.compile(
        "".join(
            ["^(?P<ref>.+):", "(?P<start>[0-9]+)-(?P<end>[0-9]+):(?P<strand>[\\+-])$"]
        )
    )

    def __init__(self, ref, start, end, strand=STRAND.PLUS):
        super(SequenceCoords, self).__init__()
        if start < 0:
            raise AssertionError
        if end < 0:
            raise AssertionError
        if not isinstance(strand, self.STRAND):
            raise AssertionError
        self._ref = ref
        self._start = start
        self._end = end
        self._strand = strand

    @property
    def ref(self):
        return self._ref

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    def __eq__(self, other):
        return all(
            [
                isinstance(other, SequenceCoords),
                self.ref == other.ref,
                self.start == other.start,
                self.end == other.end,
                self.strand == other.strand,
            ]
        )

    @staticmethod
    def rev(strand: "SequenceCoords.STRAND") -> "SequenceCoords.STRAND":
        """Provide reverse strand.

        :param strand: current strand
        :type strand: SequenceCoords.STRAND
        :return: reverse strand
        :rtype: SequenceCoords.STRAND
        """
        if SequenceCoords.STRAND.PLUS == strand:
            return SequenceCoords.STRAND.MINUS
        return SequenceCoords.STRAND.PLUS

    def __repr__(self):
        return "%s:%d-%d:%s" % (self.ref, self.start, self.end, self.strand.label)

    @staticmethod
    def from_str(s: str) -> "SequenceCoords":
        """Build a SequenceCoords object from a string.

        :param s: input string
        :type s: str
        :return: class instance
        :rtype: SequenceCoords
        :raises AssertionError: if string is not compatible with SequenceCoords.regexp
        """
        regexp_match = SequenceCoords.regexp.search(s)
        if regexp_match is None:
            raise AssertionError(f"incompatible string: {s}")
        ref, start, end, strand = regexp_match.group("ref", "start", "end", "strand")
        return SequenceCoords(
            ref,
            int(start),
            int(end),
            list(SequenceCoords.STRAND)[
                [x.label for x in list(SequenceCoords.STRAND)].index(strand)
            ],
        )


class Sequence(om.Sequence):
    """Nucleic acid sequence.

    Extends the homonym class from oligo-melting adding k-mer generators and
    batcher methods.

    Extends:
            om.Sequence

    Variables:
            doReverseComplement {bool} -- whether to generate the reverse complement
                                          when crawling through the sequence.
    """

    doReverseComplement = False

    def __init__(self, seq, t, name=None):
        if not isinstance(t, om.NATYPES):
            raise AssertionError("sequence type must be from om.NATYPES")
        super().__init__(seq, t, name)

    def kmers(self, k: int) -> Iterator["Sequence"]:
        """Extract k-mers from Sequence.

        :param k: k-mer length.
        :type k: int
        :return: Iterator.
        :rtype: Iterator[Sequence]
        """
        return self.kmerator(
            self.text, k, self.natype, self.name, rc=self.doReverseComplement
        )

    def batches(self, k: int, batchSize: int) -> Iterator[Tuple[str, int]]:
        """Split the sequence in batches ready to be fed to the kmerator.

        This is useful when parallelizing, due to the fact that generators
        cannot be pickled and can cause problems when using a process-based
        back-end. Instead, pass the sequences to the processes and run kmerator
        there.

        :param k: k-mer length
        :type k: int
        :param batchSize: records per batch
        :type batchSize: int
        :return: iterator of (header, sequence)
        :rtype: Iterator[Tuple[str, int]]
        """
        return self.batcher(self.text, k, batchSize)

    def kmers_batched(
        self, k: int, batchSize: int = 1
    ) -> Iterator[Iterator["Sequence"]]:
        """Extract batches of k-mers from Sequence.

        :param k: k-mer length
        :type k: int
        :param batchSize: records per batch, defaults to 1
        :type batchSize: int
        :yield: k-mer batch generator
        :rtype: Iterator[Sequence]
        :raises AssertionError: if batchSize is lower than 1
        """
        if batchSize < 1:
            raise AssertionError
        if batchSize == 1:
            yield self.kmers(k)
        else:
            yield from self.kmerator_batched(
                self.text,
                k,
                self.natype,
                batchSize,
                self.name,
                rc=self.doReverseComplement,
            )

    @staticmethod
    def __kmer_yielding(i, seq, prefix, k, t, offset, strand, rc) -> Iterator["KMer"]:
        yield KMer(
            prefix,
            i + offset,
            i + offset + k,
            seq[i : i + k],
            t,
            strand=strand,
        )

    @staticmethod
    def __kmer_yielding_with_rc(
        i, seq, prefix, k, t, offset, strand, rc
    ) -> Iterator["KMer"]:
        yield KMer(
            prefix,
            i + offset,
            i + offset + k,
            seq[i : i + k],
            t,
            strand=strand,
        )
        yield KMer(
            prefix,
            i + offset,
            i + offset + k,
            Sequence.mkrc(seq[i : i + k], t),
            t,
            strand=SequenceCoords.rev(strand),
        )

    @staticmethod
    def yield_kmers(
        seq: str,
        prefix: str,
        k: int,
        t: om.NATYPES,
        offset: int,
        strand: SequenceCoords.STRAND,
        rc: bool,
    ) -> Iterator["KMer"]:
        """Extract k-mers from sequence.

        :param seq: sequence
        :type seq: str
        :param prefix: reference record name
        :type prefix: str
        :param k: k-mer length
        :type k: int
        :param t: nucleic acid type
        :type t: om.NATYPES
        :param offset: current location for shifting
        :type offset: int
        :param strand: strandedness
        :type strand: SequenceCoords.STRAND
        :param rc: perform reverse-complement operation
        :type rc: bool
        :yield: k-mers
        :rtype: Iterator[Kmer]
        """
        seq = seq.upper()
        kmer_yielder = (
            Sequence.__kmer_yielding_with_rc if rc else Sequence.__kmer_yielding
        )
        for i in range(len(seq) - k + 1):
            if not Sequence.check_ab(seq[i : i + k], om.AB_NA[t]):
                logging.warning(
                    " ".join(
                        [
                            "skipped sequence with unexpected character:",
                            seq[i : i + k],
                        ]
                    )
                )
                continue
            yield from kmer_yielder(i, seq, prefix, k, t, offset, strand, rc)

    @staticmethod
    def kmerator(
        seq: str,
        k: int,
        t: om.NATYPES,
        prefix: str = "ref",
        offset: "int" = 0,
        strand: SequenceCoords.STRAND = SequenceCoords.STRAND.PLUS,
        rc=False,
    ) -> Iterator["KMer"]:
        """Extract k-mers from sequence.

        :param seq: sequence
        :type seq: str
        :param k: k-mer length
        :type k: int
        :param t: nucleic acid type
        :type t: om.NATYPES
        :param prefix: reference record name, defaults to "ref"
        :type prefix: str
        :param offset: current location for shifting, defaults to 0
        :type offset: int
        :param strand: strandedness, defaults to SequenceCoords.STRAND.PLUS
        :type strand: SequenceCoords.STRAND
        :param rc: do perform reverse-complement operation, defaults to False
        :type rc: bool
        :return: k-mer iterator
        :rtype: Iterator[KMer]
        """
        return iter(Sequence.yield_kmers(seq, prefix, k, t, offset, strand, rc))

    @staticmethod
    def batcher(seq: str, k: int, batchSize: int) -> Iterator[Tuple[str, int]]:
        """Split the sequence in batches ready to be fed to the kmerator.

        This is useful when parallelizing, due to the fact that generators
        cannot be pickled and can cause problems when using a process-based
        back-end. Instead, pass the sequences to the processes and run kmerator
        there.

        :param seq: sequence
        :type seq: str
        :param k: k-mer length
        :type k: int
        :param batchSize: records per batch
        :type batchSize: int
        :yield: k-mer sequence and start position
        :rtype: Iterator[Tuple[str, int]]
        """
        start = 0
        while start < len(seq) - k + 1:
            end = min(len(seq), start + batchSize)
            yield (seq[start:end], start)
            start += batchSize - k + 1

    @staticmethod
    def kmerator_batched(
        seq: str, k: int, t: om.NATYPES, batchSize: int = 1, prefix="ref", rc=False
    ) -> Iterator[Iterator["Sequence"]]:
        """Extract batches of k-mers from sequence.

        :param seq: sequence
        :type seq: str
        :param k: k-mer length
        :type k: int
        :param t: nucleic acid type
        :type t: om.NATYPES
        :param batchSize: records per batch, defaults to 1
        :type batchSize: int
        :param prefix: reference record name, defaults to "ref"
        :type prefix: str
        :param rc: do perform reverse-complement, defaults to False
        :type rc: bool
        :yield: k-mer iterator
        :rtype: Iterator[Kmer]
        :raises AssertionError: if batchSize is < 1
        """
        if batchSize < 1:
            raise AssertionError
        if batchSize == 1:
            yield Sequence.kmerator(seq, k, t, prefix, rc=rc)
        for (seq2beKmered, i) in Sequence.batcher(seq, k, batchSize):
            yield Sequence.kmerator(seq2beKmered, k, t, prefix, offset=i, rc=rc)


class KMer(Sequence):
    """K-mer object.

    Extends the Sequence class providing coordinate system and alphabet check.

    Extends:
            Sequence
    """

    def __init__(
        self,
        chrom,
        start,
        end,
        seq,
        t=om.NATYPES.DNA,
        strand=SequenceCoords.STRAND.PLUS,
    ):
        if len(seq) != end - start:
            raise AssertionError
        super().__init__(seq, t)
        self._coords = SequenceCoords(chrom, start, end, strand)

    @property
    def coords(self):
        return self._coords

    @property
    def header(self):
        return str(self.coords)

    @property
    def seq(self):
        return self.text

    def __eq__(self, other):
        if not self.coords == other.coords:
            return False
        return super().__eq__(other)

    @staticmethod
    def from_fasta(record: Tuple[str, str], t: om.NATYPES = om.NATYPES.DNA) -> "KMer":
        """Reads a KMer from a Fasta record.

        :param record: (header, sequence)
        :type record: Tuple[str, str]
        :param t: nucleic acid type, defaults to om.NATYPES.DNA
        :type t: om.NATYPES
        :return: k-mer instance
        :rtype: KMer
        """
        coords = SequenceCoords.from_str(record[0])
        return KMer(
            coords.ref, coords.start, coords.end, record[1], t, strand=coords.strand
        )

    @staticmethod
    def from_file(*args, **kwargs):
        return KMer.from_fasta(*args, **kwargs)

    def as_fasta(self) -> str:
        """Fasta-like representation.

        :return: fasta format
        :rtype: str
        """
        return ">%s\n%s\n" % (self.header, self.seq)

    def __repr__(self):
        return "%s\t%s" % (self.header, self.seq)

    def is_ab_checked(self) -> bool:
        """Check if AB is fully respected.

        Checks if the KMer sequence respects the AB. I.e., it does not contain
        any foreign characters.

        :return: if alphabet is respected
        :rtype: bool
        """
        return all(c in self.ab[0] for c in set(self.text))


class SequenceCount(Sequence):
    """Sequence counting system.

    Retains a sequence and the headers of all the regions where it is present.

    Extends:
            Sequence

    Variables:
            __headers {list} -- list of headers.
    """

    __headers: List[str] = []

    def __init__(self, seq, headers, t=om.NATYPES.DNA):
        super().__init__(seq, t)
        if not all(isinstance(h, str) for h in headers):
            raise AssertionError
        self.__headers = headers

    @property
    def header(self):
        return self.__headers.copy()

    @property
    def seq(self):
        return self.text

    @staticmethod
    def from_text(line: str, t: om.NATYPES = om.NATYPES.DNA) -> "SequenceCount":
        """Reads a KMer from a FASTA record.

        :param line: sequence
        :type line: str
        :param t: nucleic acid type, defaults to om.NATYPES.DNA
        :type t: om.NATYPES, optional
        :return: class instance
        :rtype: SequenceCount
        """
        seq, headers = line.strip().split("\t")
        return SequenceCount(seq, headers.split(" "), t)

    @staticmethod
    def from_file(*args, **kwargs):
        return SequenceCount.from_text(*args, **kwargs)

    def __repr__(self):
        return "%s\t%s" % (self.seq, " ".join(self.header))

    def as_text(self):
        return str(self) + "\n"
