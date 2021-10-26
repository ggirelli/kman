"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
"""

import logging
import multiprocessing as mp
import tempfile
from enum import Enum
from heapq import merge
from itertools import chain
from typing import (
    IO,
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    Type,
)

import numpy as np  # type: ignore
import oligo_melting as om  # type: ignore
from joblib import Parallel, delayed  # type: ignore
from tqdm import tqdm  # type: ignore

from kmermaid.abundance import AbundanceVector, AbundanceVectorLocal
from kmermaid.batch import Batch, BatchAppendable
from kmermaid.batcher import BatcherThreading
from kmermaid.seq import SequenceCoords, SequenceCount


class Crawler:
    """Crawling system.

    Allows for crawling along Batch records, or to group the records based on
    their sequence and crawl directly through the groups.

    Variables:
            doSort {bool} -- whether sorting should be performed while crawling
            verbose {bool} -- whether to be verbose
    """

    doSort = False
    doSmart = False
    verbose = True
    desc = ""

    @staticmethod
    def count_records(batches: List[Batch]) -> int:
        """Count records across batches.

        :param batches: list of batches
        :type batches: List[Batch]
        :return: number of records
        :rtype: int
        """
        return sum(b.current_size for b in batches)

    def do_records(self, batches: List[Batch]) -> Iterator[Tuple[str, str]]:
        """Crawl through the batches.

        Produces a generator function that yields (r.header, r.seq) for each
        record r across batches. Batches are crawled through their r.record_gen
        if self.doSort==False, otherwise using r.sorted.

        :param batches: list of batches
        :type batches: List[Type[Batch]]
        :yield: iterator across all batches
        :rtype: Iterator[Tuple[str, str]]
        :raises AssertionError: when the batches are not recognized
        """
        if any(type(b) not in [Batch, BatchAppendable] for b in batches):
            raise AssertionError()

        if self.doSort:
            generators = [
                ((str(r.header), str(r.seq)) for r in b.sorted(self.doSmart))
                for b in batches
                if type(None) != type(b)
            ]

        else:
            generators = [
                ((str(r.header), str(r.seq)) for r in b.record_gen(self.doSmart))
                for b in batches
                if type(None) is not type(b)
            ]

        yield from merge(*generators, key=lambda x: x[1])

    def do_batch(self, batches: List[Batch]) -> Iterator[Tuple[List[str], str]]:
        """Group records from batches based on sequence.

        Crawls into groups of records from input batches.

        :param batches: list of batches
        :type batches: List[Type[Batch]]
        :yield: (header, sequence)
        :rtype: Tuple[List[str], str]
        """
        crawler = self.do_records(batches)

        try:
            first_record = next(crawler)
        except StopIteration:
            logging.error("nothing to crawl")
            return

        current_seq = first_record[1]
        current_headers = [first_record[0]]

        crawler = (
            tqdm(crawler, initial=1, desc=self.desc, total=self.count_records(batches))
            if self.verbose
            else crawler
        )

        for record in crawler:
            if current_seq == record[1]:
                current_headers.append(record[0])
            else:
                yield (current_headers, current_seq)
                current_seq = record[1]
                current_headers = [record[0]]

        yield (current_headers, current_seq)


class KJoiner:
    """K-way joining system.

    Perform K-way joining of Batches of elements compatible with the Crawler
    class. See the different MODEs of joining for more details.

    Variables:
            DEFAULT_MODE {MODE} -- default join mode
            __mode {MODE} -- join mode
            __join_function {function} -- join method
    """

    class MODE(Enum):
        """Modes of joining.

        Extends:
                Enum

        Variables:
                UNIQUE {number} -- preserve only unique sequences
                SEQ_COUNT {number} -- count sequences
                VEC_COUNT {number} -- produce abundance vector
                VEC_COUNT_MASKED {number} -- produce abundance vector w/o counting
                                             abundance in same reference record
        """

        UNIQUE = 1
        SEQ_COUNT = 2
        VEC_COUNT = 3
        VEC_COUNT_MASKED = 4

    class MEMORY(Enum):
        """Modes of memory management.

        Extends:
                Enum

        Variables:
                NORMAL {number} -- store AbundanceVector in the RAM
                LOCAL {number} -- store AbundanceVector locally in HDF5 files
        """

        NORMAL = 1
        LOCAL = 2

    DEFAULT_MODE = MODE.UNIQUE
    DEFAULT_MEMORY = MEMORY.NORMAL

    __mode = DEFAULT_MODE
    __memory = DEFAULT_MEMORY
    __join_function = None

    def __init__(self, mode: MODE = None, memory: MEMORY = None):
        """Initialize KJoiner

        :param mode: joining mode, defaults to None
        :type mode: MODE
        :param memory: storage mode, defaults to None
        :type memory: MEMORY
        :raises AssertionError: if joining mode is unrecognized
        :raises AssertionError: if storage mode is unrecognized
        """
        super().__init__()
        if mode is not None:
            if mode not in self.MODE:
                raise AssertionError
            self.__mode = mode
        if memory is not None:
            if memory not in self.MEMORY:
                raise AssertionError
            self.__memory = memory
        self.__set_join_function()

    @property
    def mode(self):
        return self.__mode

    @mode.setter
    def mode(self, mode):
        if mode not in self.MODE:
            raise AssertionError
        self.__mode = mode
        self.__set_join_function()

    @property
    def memory(self):
        return self.__memory

    @memory.setter
    def memory(self, memory):
        if memory not in self.MEMORY:
            raise AssertionError
        self.__memory = memory
        self.__set_join_function()

    @property
    def join_function(self):
        return self.__join_function

    def __set_join_function(self):
        """Select the appropriate join function, based on the current mode."""
        if self.mode == self.MODE.UNIQUE:
            self.__join_function = self.join_unique
        elif self.mode == self.MODE.SEQ_COUNT:
            self.__join_function = self.join_sequence_count
        elif self.mode == self.MODE.VEC_COUNT:
            self.__join_function = self.join_vector_count
        elif self.mode == self.MODE.VEC_COUNT_MASKED:
            self.__join_function = self.join_vector_count_masked

    @staticmethod
    def join_unique(
        headers: List[str], seq: str, OH: IO, **kwargs
    ) -> Optional[Tuple[str, str]]:
        """Perform uniq joining, retain only unique records.

        :param headers: headers of records with the current sequence
        :type headers: List[str]
        :param seq: sequence
        :type seq: str
        :param OH: output file handle
        :type OH: IO
        :param **kwargs: capture additional arguments
        :return: (header, sequence)
        :rtype: Optional[Tuple[str, str]]
        """
        if len(headers) != 1:
            return None
        batch = (headers[0], seq)
        OH.write(">%s\n%s\n" % batch)
        return batch

    @staticmethod
    def join_sequence_count(
        headers: List[str], seq: str, OH: IO, **kwargs
    ) -> Tuple[str, int]:
        """Perform sequence counting through joining.

        Counts sequence occurrences.

        :param headers: headers of records with current sequence
        :type headers: List[str]
        :param seq: sequence
        :type seq: str
        :param OH: output file handle
        :type OH: IO
        :param **kwargs: capture additional arguments
        :return: (sequence, count)
        :rtype: Tuple[str, int]
        """
        batch = (seq, len(headers))
        OH.write("%s\t%d\n" % batch)
        return batch

    @staticmethod
    def join_vector_count(
        headers: List[str], seq: str, OH: IO, vector: AbundanceVector, **kwargs
    ) -> None:
        """Generate abundance vectors through joining.

        :param headers: headers of records with current sequence
        :type headers: List[str]
        :param seq: sequence
        :type seq: str
        :param OH: output file handle
        :type OH: IO
        :param **kwargs: capture additional arguments
        :param vector: vector to populate
        :type vector: AbundanceVector
        """
        hcount = len(headers)
        for header in headers:
            coords = SequenceCoords.from_str(header)
            vector.add_count(
                coords.ref, coords.strand.label, int(coords.start), hcount, len(seq)
            )

    @staticmethod
    def join_vector_count_masked(
        headers: List[str], seq: str, OH: IO, vector: AbundanceVector, **kwargs
    ) -> None:
        """Generate abundance vectors through joining.

        :param headers: headers of records with current sequence
        :type headers: List[str]
        :param seq: sequence
        :type seq: str
        :param OH: output file handle
        :type OH: IO
        :param **kwargs: capture additional arguments
        :param vector: vector to populate
        :type vector: AbundanceVector
        """
        coords = [SequenceCoords.from_str(h) for h in headers]
        if len(coords) != 1:
            refList, refCounts = np.unique([h.ref for h in coords], return_counts=True)

            if len(refList) != 1:
                for h in coords:
                    hcount = refCounts[refList != h.ref].sum()
                    vector.add_count(
                        h.ref, h.strand.label, int(h.start), hcount, len(seq)
                    )

    def _pre_join(self, outpath: str) -> Dict[str, Any]:
        """Prepares for joining.

        Perform appropriate action (mode-based):
        - open buffer to output file
        - create AbundanceVector instance

        :param outpath: path to output
        :type outpath: str
        :return: keyword arguments for join functions
        :rtype: Dict[str, Any]
        :raises AssertionError: ifmemory mode is neither NORMAL nor LOCAL
        """
        if not self.mode.name.startswith("VEC_"):
            return dict(OH=open(outpath, "w+"))

        kwargs: Dict[str, Any] = {"OH": outpath}
        if self.memory == self.MEMORY.NORMAL:
            kwargs["vector"] = AbundanceVector()
        elif self.memory == self.MEMORY.LOCAL:
            kwargs["vector"] = AbundanceVectorLocal()
        elif self.memory not in [self.MEMORY.NORMAL, self.MEMORY.LOCAL]:
            raise AssertionError
        return kwargs

    def _post_join(self, **kwargs) -> None:
        """Wraps up after joining.

        Perform appropriate actiond (mode-based):
        - close buffer to output file
        - write AbundanceVector to file

        :param **kwargs: capture additional arguments
        """
        if not self.mode.name.startswith("VEC_"):
            kwargs["OH"].close()
        else:
            kwargs["vector"].write_to(kwargs["OH"])

    def join(self, batches: List[Batch], outpath: str) -> None:
        """Perform k-joining of batches.

        :param batches: list of batches
        :type batches: List[Batch]
        :param outpath: path to output
        :type outpath: str
        """
        kwargs = self._pre_join(outpath)

        crawler = Crawler()
        print("Joining...")
        for batch in crawler.do_batch(batches):
            self.join_function(*batch, **kwargs)

        self._post_join(**kwargs)


class KJoinerThreading(KJoiner):
    """Parallelized K-joining system.

    Extends KJoiner class for parallelization.

    Extends:
            KJoiner

    Variables:
            _tmp {tempfile.TemporaryDirectory}
            _threads {number} -- number of threads for parallelization
            __batch_size {number} -- number of batches to join per thread
            __doSort {bool} -- whether to sort while joining
    """

    _tmp = None
    _threads = 1
    __batch_size = 10
    __doSort = False

    @property
    def doSort(self):
        return self.__doSort

    @doSort.setter
    def doSort(self, doSort):
        if type(True) is not type(doSort):
            raise AssertionError
        self.__doSort = doSort

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, t):
        self._threads = max(1, min(t, mp.cpu_count()))

    @property
    def batch_size(self):
        return self.__batch_size

    @batch_size.setter
    def batch_size(self, batch_size):
        if type(0) != type(batch_size):
            raise AssertionError
        if batch_size < 2:
            raise AssertionError
        self.__batch_size = batch_size

    @property
    def tmp(self):
        if self._tmp is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="kmermaidJoin")
        return self._tmp

    def __parallel_join(self, recordBatches: List[Batch], outpath: str) -> None:
        """Join sequenceCount batches in parallel.

        :param recordBatches: list of batches
        :type recordBatches: List[Batch]
        :param outpath: path to output
        :type outpath: str
        """
        kwargs = self._pre_join(outpath)

        batcher = SeqCountBatcher.from_parent(self, self.batch_size)
        batcher.doSort = self.doSort
        print("Intermediate batching...")
        batcher.do(recordBatches)
        print("Joining...")
        batcher.join(self.join_function, **kwargs)

        self._post_join(**kwargs)

    def join(self, batches: List[Batch], outpath: str) -> None:
        """Perform k-joining of batches.

        :param batches: list of batches
        :type batches: List[Batch]
        :param outpath: path to output
        :type outpath: str
        """
        if self.threads == 1:
            super().join(batches, outpath)
        else:
            self.__parallel_join(batches, outpath)


class SeqCountBatcher(BatcherThreading):
    """SequenceCount batchin system.

    Batching system to convert record batches into Batches of SequenceCount
    records, which are optimized for uniquing/counting.

    Extends:
            BatcherThreading

    Variables:
            _type {type} -- Batch record type
            __doSort {bool} -- whether to perform sorting while batching
    """

    _type = SequenceCount
    __doSort = False

    def __init__(
        self,
        n_batches: int = 10,
        threads: int = 1,
        size: int = 1000000,
        natype: om.NATYPES = om.NATYPES.DNA,
        tmp: str = tempfile.gettempdir(),
    ):
        """Initialize SeqCountBatcher instance.

        :param n_batches: number of batches at a time, defaults to 10
        :type n_batches: int
        :param threads: for parallelization, defaults to 1
        :type threads: int
        :param size: records per batch, defaults to 1000000
        :type size: int
        :param natype: nucleic acid type, defaults to om.NATYPES.DNA
        :type natype: om.NATYPES
        :param tmp: temporary folder, defaults to tempfile.gettempdir()
        :type tmp: str
        """
        super().__init__(threads, size, natype, tmp)
        self.n_batches = n_batches

    @property
    def doSort(self):
        return self.__doSort

    @doSort.setter
    def doSort(self, doSort):
        if type(True) is not type(doSort):
            raise AssertionError
        self.__doSort = doSort

    def do(self, recordBatch: List[Batch]) -> None:
        """Start batching the records.

        Batch seq.Sequence sub-class batch.Batch records into seq.SequenceCounts
        batch.Batch instances.

        :param recordBatch: list of batches
        :type recordBatch: List[Batch]
        """
        batchList = [
            recordBatch[i : min(len(recordBatch), i + self.n_batches)]
            for i in range(0, len(recordBatch), self.n_batches)
        ]
        batches = Parallel(n_jobs=self.threads, verbose=11)(
            delayed(SeqCountBatcher.build_batch)(
                batchedRecords, self.type, self.tmp, self.doSort
            )
            for batchedRecords in batchList
        )
        self.feed_collection(batches, self.FEED_MODE.REPLACE)

    @staticmethod
    def build_batch(
        recordBatchList: List[Batch],
        recordType: Type,
        tmpDir: str,
        doSort: bool = False,
    ) -> Batch:
        """Build a batch.

        :param recordBatchList: list of batches.
        :type recordBatchList: List[Batch]
        :param recordType: batched record type.
        :type recordType: Type
        :param tmpDir: path to temporary folder.
        :type tmpDir: str
        :param doSort: re-sort batches, defaults to False
        :type doSort: bool
        :return: new batch
        :rtype: Batch
        """
        crawling = Crawler()
        crawling.doSort = doSort
        crawling.doSmart = True
        crawling.verbose = False

        batch_size = max(crawling.count_records(recordBatchList), 1)
        batch = BatchAppendable(recordType, tmpDir, batch_size)
        batch.isFasta = False
        batch.suffix = ".txt"
        batch.fwrite = "as_text"
        batch.add_all(
            (
                SequenceCount(seq, headers)
                for (headers, seq) in crawling.do_batch(recordBatchList)
            )
        )
        batch.write()

        return batch

    def join(self, fjoin: Callable, **kwargs) -> None:
        """Join SequenceCount batches.

        :param fjoin: join function
        :type fjoin: Callable
        :param **kwargs: additional arguments for join function
        """
        crawler = Crawler()
        crawler.doSmart = True
        crawler.desc = "Final joining..."
        for headers, seq in crawler.do_batch(self.collection):
            headers = list(chain(*headers))
            fjoin(headers, seq, **kwargs)

    @staticmethod
    def from_parent(parent: KJoinerThreading, n_batches: int) -> "SeqCountBatcher":
        """Instantiate class from parent.

        :param parent: parent joiner
        :type parent: KJoinerThreading
        :param n_batches: number of batches
        :type n_batches: int
        :raises AssertionError: if parent type is incompatible
        :return: new instance
        :rtype: SeqCountBatcher
        """
        if type(parent) is not KJoinerThreading:
            raise AssertionError
        return SeqCountBatcher(n_batches, parent.threads, tmp=parent.tmp)
