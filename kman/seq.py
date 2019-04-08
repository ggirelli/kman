
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
'''

from enum import Enum, unique
import gzip
import io
from itertools import chain
import oligo_melting as om
import re

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

	regexp = re.compile(''.join(['^(?P<ref>.+):',
		'(?P<start>[0-9]+)-(?P<end>[0-9]+):(?P<strand>[\\+-])$']))

	def __init__(self, ref, start, end, strand = STRAND.PLUS):
		super(SequenceCoords, self).__init__()
		assert start >= 0
		assert end >= 0
		assert strand in self.STRAND
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
	
	@staticmethod
	def rev(strand):
		"""Provides reverse strand.
		
		Arguments:
			strand {SequenceCoords.STRAND} -- current strand
		
		Returns:
			Sequence.Coords.STRAND -- reverse strand
		"""
		if SequenceCoords.STRAND.PLUS == strand:
			return SequenceCoords.STRAND.MINUS
		else:
			return SequenceCoords.STRAND.PLUS

	def __repr__(self):
		return "%s:%d-%d:%s" % (self.ref,
			self.start, self.end, self.strand.label)

	@staticmethod
	def from_str(s):
		"""Builds a SequenceCoords object from a string.
		
		Arguments:
			s {str} -- input string
		
		Returns:
			SequenceCoords
		"""
		ref, start, end, strand = SequenceCoords.regexp.search(s
			).group("ref", "start", "end", "strand")
		strand = list(SequenceCoords.STRAND)[[x.label
			for x in list(SequenceCoords.STRAND)].index(strand)]
		return SequenceCoords(ref, int(start), int(end), strand)

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

	def __init__(self, seq, t, name = None):
		super().__init__(seq, t, name)

	def kmers(self, k):
		"""Extract k-mers from Sequence.
		Args:
			k {int} -- substring length
		Returns:
			generator -- kmer generator
		"""
		return self.kmerator(self.text, k, self.natype, self.name,
			rc = self.doReverseComplement)

	def batches(self, k, batchSize):
		"""Split the sequence in batches ready to be fed to the kmerator.
		
		This is useful when parallelizing, due to the fact that generators
		cannot be pickled and can cause problems when using a process-based
		back-end. Instead, pass the sequences to the processes and run kmerator
		there.
		
		Arguments:
			seq {string} -- sequence to be batched
			k {int} -- length of substrings for kmerator
			batchSize {int} -- number of kmers per batch
		"""
		return self.batcher(self.text, k, batchSize)

	def kmers_batched(self, k, batchSize = 1):
		"""Extract batches of k-mers from Sequence.
		
		Arguments:
			k {int} -- substring length

		Keyword Arguments:
			batchSize {number} -- number of k-mers per batch (default: {1})

		Returns:
			generator -- k-mer batch generator
		"""
		assert batchSize >= 1
		if batchSize == 1:
			return self.kmers(k)
		else:
			return self.kmerator_batched(self.text, k, self.natype,
				batchSize, self.name, rc = self.doReverseComplement)

	@staticmethod
	def kmerator(seq, k, t, prefix = "ref", offset = 0,
		strand = SequenceCoords.STRAND.PLUS, rc = False):
		"""Extract k-mers from seq.
		
		Arguments:
			seq {string} -- input sequence
			k {int} -- substring length
			t {om.NATYPES} -- nucleic acid type
		
		Keyword Arguments:
			prefix {str} -- reference record name (default: {"ref"})
			offset {number} -- if this is a batch, current location for
			                   shifting (default: {0})
		"""
		def kmerGen(seq, prefix, k, t, offset, strand, rc):
			if rc:
				revStrand = SequenceCoords.rev(strand)
				for i in range(len(seq)-k+1):
					yield KMer(prefix, i+offset, i+offset+k, seq[i:i+k], t,
						strand = strand)
					yield KMer(prefix, i+offset, i+offset+k,
						Sequence.mkrc(seq[i:i+k], t), t, strand = revStrand)
			else:
				for i in range(len(seq)-k+1):
					yield KMer(prefix, i+offset, i+offset+k, seq[i:i+k], t,
						strand = strand)
		return kmerGen(seq, prefix, k, t, offset, strand, rc)

	@staticmethod
	def batcher(seq, k, batchSize):
		"""Split the sequence in batches ready to be fed to the kmerator.
		
		This is useful when parallelizing, due to the fact that generators
		cannot be pickled and can cause problems when using a process-based
		back-end. Instead, pass the sequences to the processes and run kmerator
		there.
		
		Arguments:
			seq {string} -- sequence to be batched
			k {int} -- length of substrings for kmerator
			batchSize {int} -- number of kmers per batch
		"""
		for i in range(0, len(seq)-k+1, batchSize):
			yield (seq[max(i, i-k+1):min(len(seq)-k+1, i+batchSize)], i)

	@staticmethod
	def kmerator_batched(seq, k, t, batchSize = 1, prefix = "ref", rc = False):
		"""Extract batches of k-mers from seq.
		
		Arguments:
			seq {string} -- input sequence
			k {int} -- substring length
			t {om.NATYPES} -- nucleic acid type
		
		Keyword Arguments:
			batchSize {number} -- number of kmers per batch (default: {1})
			prefix {str} -- reference record name (default: {"ref"})
		
		Returns:
			generator -- k-mer batch generator
		"""
		assert batchSize >= 1
		if batchSize == 1:
			return Sequence.kmerator(seq, k, t, prefix, rc = rc)
		else:
			for (seq2beKmered, i) in Sequence.batcher(seq, k, batchSize):
				yield Sequence.kmerator(seq2beKmered, k, t, prefix,
					offset = i, rc = rc)

class KMer(Sequence):
	"""K-mer object.
	
	Extends the Sequence class providing coordinate system and alphabet check.
	
	Extends:
		Sequence
	"""

	def __init__(self, chrom, start, end, seq,
		t = om.NATYPES.DNA, strand = SequenceCoords.STRAND.PLUS):
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

	@staticmethod
	def from_fasta(record, t = om.NATYPES.DNA):
		"""Reads a KMer from a Fasta record.
		
		Arguments:
			record {tuple} -- (header, seq)
		
		Keyword Arguments:
			t {om.NATYPES} -- nucleic acid type (default: {om.NATYPES.DNA})
		
		Returns:
			KMer
		"""
		coords = SequenceCoords.from_str(record[0])
		return KMer(coords.ref, coords.start, coords.end,
			record[1], t, strand = coords.strand)
	
	@staticmethod
	def from_file(*args, **kwargs):
		return KMer.from_fasta(*args, **kwargs)

	def as_fasta(self):
		"""Fasta-like representation."""
		return ">%s\n%s\n" % (self.header, self.seq)

	def __repr__(self):
		return "%s\t%s" % (self.header, self.seq)

	def is_ab_checked(self):
		"""Check if AB is fully respected.
		
		Checks if the KMer sequence respects the AB. I.e., it does not contain
		any foreign characters.
		
		Returns:
			bool -- whether AB is respected.
		"""
		for c in set(self.text):
			if not c in self.ab[0]:
				return False
		return True

class SequenceCount(Sequence):
	"""Sequence counting system.
	
	Retains a sequence and the headers of all the regions where it is present.
	
	Extends:
		Sequence
	
	Variables:
		__headers {list} -- list of headers.
	"""

	__headers = []

	def __init__(self, seq, headers, t = om.NATYPES.DNA):
		super().__init__(seq, t)
		self.__headers = headers

	@property
	def header(self):
		return self.__headers.copy()
	@property
	def seq(self):
		return self.text

	@staticmethod
	def from_text(line, t = om.NATYPES.DNA):
		"""Reads a KMer from a Fasta record.
		
		Arguments:
			record {tuple} -- (header, seq)
		
		Keyword Arguments:
			t {om.NATYPES} -- nucleic acid type (default: {om.NATYPES.DNA})
		
		Returns:
			KMer
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
