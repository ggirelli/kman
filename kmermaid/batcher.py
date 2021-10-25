"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batching
"""

import gzip
import itertools
import logging
import multiprocessing as mp
import os
import tempfile
from enum import Enum
from typing import IO, Any, List, Optional, Tuple, Type, Union

import oligo_melting as om  # type: ignore
from joblib import Parallel, delayed  # type: ignore
from tqdm import tqdm  # type: ignore

from kmermaid.batch import Batch
from kmermaid.parsers import SmartFastaParser
from kmermaid.seq import KMer, Sequence

TMP_DIR = tempfile.TemporaryDirectory


class BatcherBase(object):
    """Basic batching system.

    Builds a collection of equally sized batches. In each Batch, a record is an
    instance of a class of the given size. Batch record type is consistent
    across the collection.

    Variables:
            DEFAULT_BATCH_SIZE {int} -- default batch size
            DEFAULT_NATYPE {om.NATYPES} -- default nucleic acid type
            _tmp {TMP_DIR}
            _batches {list} -- list of Batch instances
            __size {int} -- batch size
            _type {type} -- batched record type
            __natype {om.NATYPES} -- nucleic acid type
    """

    DEFAULT_BATCH_SIZE = int(1e6)
    DEFAULT_NATYPE = om.NATYPES.DNA

    _tmpH: TMP_DIR
    _tmp: str
    _batches: List[Batch]
    __size = DEFAULT_BATCH_SIZE
    _type: Type[Sequence] = KMer
    __natype = DEFAULT_NATYPE

    def __init__(
        self,
        size: int,
        natype: Optional[om.NATYPES] = None,
        tmp: Optional[Union[TMP_DIR, str]] = None,
    ):

        """Initializes BatcherBase

        :param size: batching size
        :type size: int
        :param natype: nucleic acid type, defaults to None
        :type natype: Optional[om.NATYPES], optional
        :param tmp: temporary directory, defaults to None
        :type tmp: Optional[Union[TMP_DIR, str]], optional
        """
        super().__init__()
        self.size = size
        self.natype = natype
        if isinstance(tmp, TMP_DIR):
            self._tmpH = tmp
            self._tmp = tmp.name
        elif isinstance(tmp, str):
            self._tmp = tmp
        else:
            self._tmp = tempfile.gettempdir()
        self._batches = [Batch.from_batcher(self.type, self.size, self.tmp)]

    @property
    def size(self) -> int:
        return self.__size

    @size.setter
    def size(self, size: Union[int, None]) -> None:
        if size is not None:
            if size < 1:
                raise AssertionError
            self.__size = size

    @property
    def type(self) -> om.NATYPES:
        return self._type

    @property
    def natype(self):
        return self.__natype

    @natype.setter
    def natype(self, natype: Optional[om.NATYPES]) -> None:
        if natype is not None:
            if natype not in om.NATYPES:
                raise AssertionError
            self.__natype = natype

    @property
    def collection(self):
        return self._batches

    @property
    def tmp(self):
        if self._tmp is None:
            self._tmpH = TMP_DIR(prefix="kmermaidBatch")
            self._tmp = self._tmpH.name
        return self._tmp

    def new_batch(self):
        """Add a new empty batch to the current collection."""
        if self.collection[-1].is_full():
            self.collection[-1].write()
            self._batches.append(Batch.from_batcher(self.type, self.size, self.tmp))

    def add_record(self, record: Any):
        """Add a record to the current collection.

        :param record: to be added.
        :type record: Any
        """
        self.new_batch()  # Add new batch if needed
        self.collection[-1].add(record)

    def write_all(
        self, f: str = "as_fasta", doSort: bool = False, verbose: bool = False
    ):
        """Write all batches to disk.

        :param f: name of string representation method, defaults to "as_fasta"
        :type f: str
        :param doSort: sort while writing, defaults to False
        :type doSort: bool
        :param verbose: be verbose, defaults to False
        :type verbose: bool
        """
        biList = range(len(self.collection))
        description = "Writing"
        if doSort:
            description += "&Sorting"
        if verbose:
            biList = tqdm(biList, desc=description)
        for bi in biList:
            if self.collection[bi].current_size != 0:
                self.collection[bi].write(f, doSort)


class BatcherThreading(BatcherBase):
    """Parallelized batching system.

    Extends BatcherBase for parallelization.

    Extends:
            BatcherBase

    Variables:
            __threads {number} -- number of threads for parallelization
    """

    class FEED_MODE(Enum):
        """Feeding modes.

        Used with feed_collection() method.

        Extends:
                Enum

        Variables:
                REPLACE {number} -- replace current collection.
                FLOW {number} -- flow records into current collection.
                APPEND {number} -- append to current collection.
        """

        REPLACE = 1
        FLOW = 2
        APPEND = 3

    __threads = 1

    def __init__(
        self,
        size: int,
        threads: int = 1,
        natype: Optional[om.NATYPES] = None,
        tmp: Optional[str] = None,
    ):
        """Initialize BatcherThreading

        :param size: records per batch
        :type size: int
        :param threads: for parallelization, defaults to 1
        :type threads: int
        :param natype: nucleic acid type, defaults to None
        :type natype: Optional[om.NATYPE], optional
        :param tmp: temporary folder, defaults to None
        :type tmp: Optional[str], optional
        """
        super().__init__(size, natype, tmp)
        self.threads = threads

    @property
    def threads(self):
        return self.__threads

    @threads.setter
    def threads(self, t):
        self.__threads = max(1, min(t, mp.cpu_count()))

    def __flow_batches(self, collection) -> None:
        for _ in tqdm(range(len(collection)), desc="Flowing"):
            batch = collection.pop()
            for record in batch.record_gen():
                self.add_record(record)
            batch.reset()

    def feed_collection(
        self,
        new_collection: List[Batch],
        mode: "BatcherThreading.FEED_MODE" = FEED_MODE.FLOW,
    ):
        """Feed a new collection to the batcher

        :param new_collection: new collection
        :type new_collection: List[Batch]
        :param mode: feed mode, defaults to FEED_MODE.FLOW
        :type mode: BatcherThreading
        :raises AssertionError: if any batch type is not compatible
        """
        if any(b.type != self.type for b in new_collection):
            raise AssertionError
        if mode == self.FEED_MODE.REPLACE:
            self._batches = new_collection
        elif mode == self.FEED_MODE.FLOW:
            self.__flow_batches(new_collection)
        elif mode == self.FEED_MODE.APPEND:
            self._batches.extend(new_collection)

    @staticmethod
    def from_files(
        dirPath: str,
        threads: int,
        t: Type = KMer,
        isFasta: bool = True,
        reSort: bool = False,
    ) -> List[Batch]:
        """Load files from a directory into a list of Batch objects.

        :param dirPath: path to folder with batches
        :type dirPath: str
        :param threads: for parallelization
        :type threads: int
        :param t: record type, defaults to KMer
        :type t: Type
        :param isFasta: is the input fasta, defaults to True
        :type isFasta: bool
        :param reSort: sort while reading, defaults to False
        :type reSort: bool
        :return: list of read batches
        :rtype: List[Batch]
        :raises AssertionError: if input is not a folder
        """
        if not os.path.isdir(dirPath):
            raise AssertionError
        threads = max(1, min(threads, mp.cpu_count()))
        if threads == 1:
            return [
                Batch.from_file(os.path.join(dirPath, fname), t, isFasta)
                for fname in tqdm(os.listdir(dirPath))
            ]
        return Parallel(n_jobs=threads, verbose=11)(
            delayed(Batch.from_file)(
                os.path.join(dirPath, fname), t, isFasta, reSort=reSort
            )
            for fname in os.listdir(dirPath)
        )


class FastaBatcher(BatcherThreading):
    """FASTA file k-mer batching.

    Divides k-mer from the records of a FASTA file into batches.

    Extends:
            BatcherThreading

    Variables:
            _doReverseComplement {bool} -- whether to batch also the reverse
                                           complement of the sequences
    """

    class MODE(Enum):
        KMERS = 1
        RECORDS = 2

    _doReverseComplement: bool
    _mode: MODE

    def __init__(
        self,
        scan_mode: MODE = MODE.KMERS,
        reverse: bool = False,
        threads: int = 1,
        size: int = BatcherThreading.DEFAULT_BATCH_SIZE,
        natype: om.NATYPES = BatcherThreading.DEFAULT_NATYPE,
        tmp: str = tempfile.gettempdir(),
    ):
        """Initialize FastaBatcher.

        :param scan_mode: scanning mode, defaults to MODE.KMERS
        :type scan_mode: MODE
        :param reverse: reverse-complement sequences, defaults to False
        :type reverse: bool
        :param threads: for parallelization, defaults to 1
        :type threads: int
        :param size: records per batch, defaults to BatcherThreading.DEFAULT_BATCH_SIZE
        :type size: int
        :param natype: nucleic acid type, defaults to BatcherThreading.DEFAULT_NATYPE
        :type natype: om.NATYPES
        :param tmp: temporary folder, defaults to tempfile.gettempdir()
        :type tmp: str
        """
        super().__init__(threads, size, natype, tmp)
        self.mode = scan_mode
        self.doReverseComplement = reverse

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, m):
        if m not in self.MODE:
            raise AssertionError
        self._mode = m

    @property
    def doReverseComplement(self):
        return self._doReverseComplement

    @doReverseComplement.setter
    def doReverseComplement(self, rc):
        if type(True) != type(rc):
            raise AssertionError
        self._doReverseComplement = rc

    def __do_over_kmers(
        self,
        FH: IO,
        k: int,
        feedMode: BatcherThreading.FEED_MODE = BatcherThreading.FEED_MODE.APPEND,
    ):
        """Parallelize over kmers.

        :param FH: Fasta file handle
        :type FH: IO
        :param k: k-mer length
        :type k: int
        :param feedMode: feeding mode, defaults to BatcherThreading.FEED_MODE.APPEND
        :type feedMode: BatcherThreading.FEED_MODE, optional
        """
        batcher = FastaRecordBatcher.from_parent(self)
        for record in SmartFastaParser(FH).parse():
            batcher.do(record, k)
            for batch in batcher.collection:
                batch.unwrite()
        self.feed_collection(batcher.collection, feedMode)
        self.write_all(doSort=True, verbose=True)

    def __do_over_records(
        self,
        FH: IO,
        k: int,
        feedMode: BatcherThreading.FEED_MODE = BatcherThreading.FEED_MODE.APPEND,
    ):
        """Parallelize over FASTA records.

        Run a non-parallelized RecordBatcher for each fasta record, in parallel.

        :param FH: fasta file handle.
        :type FH: IO
        :param k: k-mer length
        :type k: int
        :param feedMode: feeding mode, defaults to BatcherThreading.FEED_MODE.APPEND
        :type feedMode: BatcherThreading.FEED_MODE, optional
        """

        def do_record(
            size: int, natype: om.NATYPES, tmp: str, record: Tuple[str, str], k: int
        ) -> List[Batch]:
            """Batch a single record.

            Function to be passed to parallel(delayed(...))

            :param size: records per batcher
            :type size: int
            :param natype: nucleic acid type
            :type natype: om.NATYPES
            :param tmp: temporary folder
            :type tmp: str
            :param record: fasta record
            :type record: Tuple[str, str]
            :param k: k-mer length
            :type k: int
            :return: List[Batch]
            """
            batcher = FastaRecordBatcher(1, size, natype, tmp)
            batcher.do(record, k, False)
            return batcher.collection

        batchCollections = Parallel(n_jobs=self.threads, verbose=11)(
            delayed(do_record)(self.size, self.natype, self.tmp, record, k)
            for record in SmartFastaParser(FH).parse()
        )

        self.feed_collection(list(itertools.chain(*batchCollections)), feedMode)

        def do_sort_write(b):
            b.write(True)

        Parallel(n_jobs=self.threads, verbose=11)(
            delayed(do_sort_write)(b) for b in self.collection if b.current_size != 0
        )

    def do(
        self,
        fasta: str,
        k: int,
        feedMode: BatcherThreading.FEED_MODE = BatcherThreading.FEED_MODE.APPEND,
    ):
        """Start batching a fasta file.

        Batch a fasta file up to the specified number of k-mers.

        :param fasta: path to fasta file
        :type fasta: str
        :param k: k-mer length
        :type k: int
        :param feedMode: feeding mode, defaults to BatcherThreading.FEED_MODE.APPEND
        :type feedMode: BatcherThreading.FEED_MODE, optional
        :raises AssertionError: if input is not a file
        :raises AssertionError: if k is lower than or equal to 1
        """
        if not os.path.isfile(fasta):
            raise AssertionError(f"input file not found: {fasta}")
        if k <= 1:
            raise AssertionError(f"k must be >= 1, got {k} instead.")

        FH = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, "r+")
        if self._mode == self.MODE.KMERS:
            self.__do_over_kmers(FH, k, feedMode)
        elif self._mode == self.MODE.RECORDS:
            self.__do_over_records(FH, k, feedMode)

        FH.close()


class FastaRecordBatcher(BatcherThreading):
    """FASTA record batchin system.

    Divides k-mer from a single FASTA record into batches.

    Extends:
            BatcherThreading

    Variables:
            _doReverseComplement {bool} -- whether to batch also the reverse
                                           complement of the sequences
    """

    _doReverseComplement = False

    def __init__(
        self,
        size: int,
        threads: int = 1,
        natype: om.NATYPES = om.NATYPES.DNA,
        tmp: str = tempfile.gettempdir(),
    ):
        """Initialize FastaRecordBatcher.

        :param threads: for parallelization, defaults to 1
        :type threads: int
        :param size: records per batch, defaults to None
        :type size: int
        :param natype: nucleic acid type, defaults to None
        :type natype: om.NATYPES
        :param tmp: temporary folder, defaults to None
        :type tmp: str
        """
        super().__init__(threads, size, natype, tmp)

    @property
    def doReverseComplement(self):
        return self._doReverseComplement

    @doReverseComplement.setter
    def doReverseComplement(self, rc):
        if type(True) != type(rc):
            raise AssertionError
        self._doReverseComplement = rc

    def do(self, record: Tuple[str, str], k: int, verbose: bool = True):
        """Start batching a fasta record.

        Requires a fasta record with header and sequence.

        :param record: (header, sequence)
        :type record: Tuple[str, str]
        :param k: k-mer length
        :type k: int
        :param verbose: be verbose, defaults to True
        :type verbose: bool
        """
        record_name = record[0].split(" ")[0]
        if verbose:
            logging.info(f"Batching record '{record_name}'...")
        if self.threads == 1:
            kmerGen = Sequence.kmerator(
                record[1], k, self.natype, record_name, rc=self.doReverseComplement
            )
            kmerGen = tqdm(kmerGen) if verbose else kmerGen
            for kmer in kmerGen:
                if kmer.is_ab_checked():
                    self.add_record(kmer)
        else:
            batches = Parallel(n_jobs=self.threads, verbose=11)(
                delayed(FastaRecordBatcher.build_batch)(seq, record_name, k, self, i)
                for (seq, i) in Sequence.batcher(record[1], k, self.size)
            )
            self.feed_collection(batches, self.FEED_MODE.APPEND)
        self.write_all()

    @staticmethod
    def build_batch(
        seq: str, name: str, k: int, batcher: "FastaRecordBatcher", i: int = 0
    ) -> Batch:
        """Build a Batch.

        :param seq: sequence
        :type seq: str
        :param name: header
        :type name: str
        :param k: k-mer length
        :type k: int
        :param batcher: parent
        :type batcher: FastaRecordBatcher
        :param i: position offset, defaults to 0
        :type i: int
        :return: built Batch
        :rtype: Batch
        """
        batch = Batch.from_batcher(batcher.type, batcher.size, batcher.tmp)
        recordGen = Sequence.kmerator(
            seq, k, batcher.natype, name, i, rc=batcher.doReverseComplement
        )
        batch.add_all((k for k in recordGen if k.is_ab_checked()))
        batch.write(doSort=True)
        return batch

    @staticmethod
    def from_parent(parent: "FastaBatcher") -> "FastaRecordBatcher":
        """Initialize FastaRecordBatcher.

        A parent batcher class can be specified, whose attributes are inherited.

        :param parent: parent batcher to inherit attributes from
        :type parent: FastaBatcher
        :return: FastaRecordBatcher
        :rtype: FastaRecordBatcher
        """
        batcher = FastaRecordBatcher(
            parent.threads, parent.size, parent.natype, parent.tmp
        )
        batcher._doReverseComplement = parent.doReverseComplement
        return batcher


def load_batches(
    previous_batches: str, threads: int = 1, re_sort: bool = False
) -> List[Batch]:
    """Load previously generated batches.

    :param previous_batches: path to folder containing previous batches.
    :type previous_batches: str
    :param threads: for parallelization, defaults to 1
    :type threads: int
    :param re_sort: re-sort while reading, defaults to False
    :type re_sort: bool
    :raises AssertionError: input folder must exist and be non-empty
    :return: read batches
    :rtype: List[Batch]
    """
    if not os.path.isdir(previous_batches) or len(os.listdir(previous_batches)) > 0:
        raise AssertionError(
            f"folder with previous batches empty or not found: {previous_batches}"
        )
    logging.info(f"Loading previous batches from '{previous_batches}'...")
    return BatcherThreading.from_files(previous_batches, threads, reSort=re_sort)
