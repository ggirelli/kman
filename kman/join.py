"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for batch joining
"""

import multiprocessing as mp
import tempfile
from enum import Enum
from heapq import merge
from itertools import chain

import numpy as np  # type: ignore
from joblib import Parallel, delayed  # type: ignore
from tqdm import tqdm  # type: ignore

from kman.abundance import AbundanceVector, AbundanceVectorLocal
from kman.batch import Batch, BatchAppendable
from kman.batcher import BatcherThreading
from kman.seq import SequenceCoords, SequenceCount


class Crawler(object):
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

    def __init__(self):
        super().__init__()

    def count_records(self, batches):
        """Count records across batches.

        Arguments:
                batches {list} -- list of Batches

        Returns:
                int -- number of records
        """
        return sum(b.current_size for b in batches)

    def do_records(self, batches):
        """Crawl through the batches.

        Produces a generator function that yields (r.header, r.seq) for each
        record r across batches. Batches are crawled through their r.record_gen
        if self.doSort==False, otherwise using r.sorted.

        Arguments:
                batches {list} -- list of Batches

        Returns:
                generator -- record generator
        """
        assert all(type(b) in [Batch, BatchAppendable] for b in batches)

        if self.doSort:
            generators = [
                ((r.header, r.seq) for r in b.sorted(self.doSmart))
                for b in batches
                if type(None) != type(b)
            ]

        else:
            generators = [
                ((r.header, r.seq) for r in b.record_gen(self.doSmart))
                for b in batches
                if not type(None) == type(b)
            ]

        return merge(*generators, key=lambda x: x[1])

    def do_batch(self, batches):
        """Group records from batches based on sequence.

        Crawls into groups of records from input batches.

        Arguments:
                batches {list} -- list of Batches

        Yields:
                tuple -- (headers, sequence)
        """
        crawler = self.do_records(batches)

        first_record = next(crawler)
        current_seq = first_record[1]
        current_headers = [first_record[0]]

        if self.verbose:
            crawler = tqdm(
                crawler, initial=1, desc=self.desc, total=self.count_records(batches)
            )

        for record in crawler:
            if current_seq == record[1]:
                current_headers.append(record[0])
            else:
                yield (current_headers, current_seq)
                current_seq = record[1]
                current_headers = [record[0]]

        yield (current_headers, current_seq)


class KJoiner(object):
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

    def __init__(self, mode=None, memory=None):
        """Initialize KJoiner.

        Keyword Arguments:
                mode {KJoiner.MODE} -- (default: {None})
        """
        super().__init__()
        if mode is not None:
            assert mode in self.MODE
            self.__mode = mode
        if memory is not None:
            assert memory in self.MEMORY
            self.__memory = memory
        self.__set_join_function()

    @property
    def mode(self):
        return self.__mode

    @mode.setter
    def mode(self, mode):
        assert mode in self.MODE
        self.__mode = mode
        self.__set_join_function()

    @property
    def memory(self):
        return self.__memory

    @memory.setter
    def memory(self, memory):
        assert memory in self.MEMORY
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
    def join_unique(headers, seq, OH, **kwargs):
        """Perform unique joining.

        Retains only unique records.

        Arguments:
                OH {io.TextIOWrapper} -- buffer to output file
                headers {list} -- list of headers with seq
                seq {str} -- sequence
        """
        if len(headers) == 1:
            batch = (headers[0], seq)
            OH.write(">%s\n%s\n" % batch)
            return batch

    @staticmethod
    def join_sequence_count(headers, seq, OH, **kwargs):
        """Perform sequence counting through joining.

        Counts sequence occurrences.

        Arguments:
                OH {io.TextIOWrapper} -- buffer to output file
                headers {list} -- list of headers with seq
                seq {str} -- sequence
        """
        batch = (seq, len(headers))
        OH.write("%s\t%d\n" % batch)
        return batch

    @staticmethod
    def join_vector_count(headers, seq, OH, vector, **kwargs):
        """Generate abundance vectors through joining.

        Arguments:
                OH {io.TextIOWrapper} -- buffer to output file
                headers {list} -- list of headers with seq
                seq {str} -- sequence
                vector {AbundanceVector}
        """
        hcount = len(headers)
        for header in headers:
            coords = SequenceCoords.from_str(header)
            vector.add_count(
                coords.ref, coords.strand.label, int(coords.start), hcount, len(seq)
            )

    @staticmethod
    def join_vector_count_masked(headers, seq, OH, vector, **kwargs):
        """Generate abundance vectors through joining.

        Arguments:
                OH {io.TextIOWrapper} -- buffer to output file
                headers {list} -- list of headers with seq
                seq {str} -- sequence
                vector {AbundanceVector}
        """
        headers = [SequenceCoords.from_str(h) for h in headers]
        if len(headers) != 1:
            refList, refCounts = np.unique([h.ref for h in headers], return_counts=True)

            if len(refList) != 1:
                for h in headers:
                    hcount = refCounts[refList != h.ref].sum()
                    vector.add_count(
                        h.ref, h.strand.label, int(h.start), hcount, len(seq)
                    )

    def _pre_join(self, outpath):
        """Prepares for joining.

        Perform appropriate actiond (mode-based):
        - open buffer to output file
        - create AbundanceVector instance

        Arguments:
                outpath {str} -- path to output

        Returns:
                dict -- keyword arguments for join function
        """
        kwargs = {"OH": outpath}
        if not self.mode.name.startswith("VEC_"):
            kwargs["OH"] = open(outpath, "w+")
        elif self.memory == self.MEMORY.NORMAL:
            kwargs["vector"] = AbundanceVector()
        elif self.memory == self.MEMORY.LOCAL:
            kwargs["vector"] = AbundanceVectorLocal()
        else:
            assert self.memory in [self.MEMORY.NORMAL, self.MEMORY.LOCAL]
        return kwargs

    def _post_join(self, **kwargs):
        """Wraps up after joining.

        Perform appropriate actiond (mode-based):
        - close buffer to output file
        - write AbundanceVector to file

        Arguments:
                **kwargs {dict} -- join function keyword arguments
        """
        if not self.mode.name.startswith("VEC_"):
            kwargs["OH"].close()
        else:
            kwargs["vector"].write_to(kwargs["OH"])

    def join(self, batches, outpath, doSort=False):
        """Join batches.

        Perform k-joining of batches.

        Arguments:
                batches {list} -- list of Batch instances
                outpath {str} -- path to output file

        Keyword Arguments:
                doSort {bool} -- whether batches need to be sorted
                                                 (default: {False})
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

    def __init__(self, mode=None, memory=None):
        super().__init__(mode, memory)

    @property
    def doSort(self):
        return self.__doSort

    @doSort.setter
    def doSort(self, doSort):
        assert type(True) == type(doSort)
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
        assert type(0) == type(batch_size)
        assert batch_size >= 2
        self.__batch_size = batch_size

    @property
    def tmp(self):
        if self._tmp is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="kmanJoin")
        return self._tmp

    def __parallel_join(self, recordBatches, outpath):
        """Joins sequenceCount batches in paralle.

        Arguments:
                recordBatches {list} -- list of Batches
                outpath {str} -- path to output
        """
        kwargs = self._pre_join(outpath)

        batcher = SeqCountBatcher.from_parent(self, self.batch_size)
        batcher.doSort = self.doSort
        print("Intermediate batching...")
        batcher.do(recordBatches)
        print("Joining...")
        batcher.join(self.join_function, **kwargs)

        self._post_join(**kwargs)

    def join(self, batches, outpath):
        """Join batches.

        Perform k-joining of batches.

        Arguments:
                batches {list} -- list of Batch instances
                outpath {str} -- path to output file
        """
        if self.threads == 1:
            super().join(batches, outpath, self.doSort)
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

    def __init__(self, n_batches=10, threads=1, size=None, natype=None, tmp=None):
        """Initialize SeqCountBatcher instance.

        Keyword Arguments:
                n_batches {number} -- number of batches at a time (default: {10})
                threads {int} -- number of threads for parallelization
                                 (default: {1})
                size {int} -- batch size
                natype {om.NATYPES} -- nucleic acid type
                tmp {tempfile.TemporaryDirectory}
        """
        super().__init__(threads, size, natype, tmp)
        self.n_batches = n_batches

    @property
    def doSort(self):
        return self.__doSort

    @doSort.setter
    def doSort(self, doSort):
        assert type(True) == type(doSort)
        self.__doSort = doSort

    def do(self, recordBatch):
        """Start batching the records.

        Batch seq.Sequence sub-class batch.Batch records into seq.SequenceCounts
        batch.Batch instances.

        Arguments:
                recordBatch {list} -- list of Batches
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
    def build_batch(recordBatchList, recordType, tmpDir, doSort=False):
        """Builds a Batch.

        Arguments:
                recordBatchList {list} -- list of Batches
                recordType {class} -- batch record type
                tmpDir {str} -- path to temporary directory

        Returns:
                Batch
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

    def join(self, fjoin, **kwargs):
        """Joins SequenceCount batches.

        Arguments:
                fjoin {function} -- join function
                **kwargs {dict} -- join function keyword arguments
        """
        crawler = Crawler()
        crawler.doSmart = True
        crawler.desc = "Final joining..."
        for headers, seq in crawler.do_batch(self.collection):
            headers = list(chain(*headers))
            fjoin(headers, seq, **kwargs)

    @staticmethod
    def from_parent(parent, n_batches):
        assert KJoinerThreading == type(parent)
        return SeqCountBatcher(n_batches, parent.threads, tmp=parent.tmp)
