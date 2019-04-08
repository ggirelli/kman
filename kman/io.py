
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
'''

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import io

class SmartFastaParser(object):
	"""Fasta parser with minimally open buffer.
	
	Extends the SimpleFastaParser function to avoid 'too many open files' issue.
	The input buffer is opened only when needed, and the byte-based position in
	the Fasta file is kept in memory. This continuous opening/closing adds a bit
	of an overhead but is useful when many Fasta files must be parsed at the
	same time.

	Variables:
		__compressed {bool} -- whether the Fasta is compressed.
		__pos {number} -- byte-based location in the Fasta file.
	"""

	__compressed = False
	__pos = 0

	def __init__(self, FH):
		super(SmartFastaParser, self).__init__()
		if str == type(FH):
			if FH.endswith(".gz"):
				self.__compressed = True
				self.__FH = gzip.open(FH, "rt")
			else:
				self.__FH = open(FH, "r+")
		elif io.TextIOWrapper == type(FH):
			self.__FH = FH
			if self.__FH.name.endswith(".gz"):
				self.__compressed = True
		else:
			assert False, "type error."


	def __reopen(self):
		"""Re-open the buffer and seek the last recorded position."""
		if self.__FH.closed:
			if self.__compressed:
				self.__FH = gzip.open(self.__FH.name, "rt")
			else:
				self.__FH = open(self.__FH.name, "r+")
		self.__FH.seek(self.__pos)

	def parse(self):
		"""Iterate over Fasta records as string tuples.

		For each record a tuple of two strings is returned, the FASTA title
		line (without the leading '>' character), and the sequence (with any
		whitespace removed). The title line is not divided up into an
		identifier (the first word) and comment or description.

		Additionally, keep the Fasta handler open only when strictly necessary.

		>>> with open("Fasta/dups.fasta") as handle:
		...     for values in SimpleFastaParser(handle):
		...         print(values)
		...
		('alpha', 'ACGTA')
		('beta', 'CGTC')
		('gamma', 'CCGCC')
		('alpha (again - this is a dup entry to test the idx code)', 'ACGTA')
		('delta', 'CGCGC')
		
		"""
		# Skip any text before the first record (e.g. blank lines, comments)
		self.__reopen()

		while True:
			line = self.__FH.readline()
			self.__pos = self.__FH.tell()
			if line == "":
				return  # Premature end of file, or just empty?
			if line[0] == ">":
				break

		while True:
			self.__reopen()
			if type(None) == type(line):
				line = self.__FH.readline()

			if line[0] != ">":
				raise ValueError(
					"Records in Fasta files should start with '>' character")
			title = line[1:].rstrip()
			lines = []
			line = self.__FH.readline()
			while True:
				if not line:
					break
				if line[0] == ">":
					break
				else:
					self.__pos = self.__FH.tell()
				lines.append(line.rstrip())
				line = self.__FH.readline()

			# Remove trailing whitespace, and any internal spaces
			# (and any embedded \r which are possible in mangled files
			# when not opened in universal read lines mode)
			
			self.__FH.close()
			yield title, "".join(lines).replace(" ", "").replace("\r", "")

			if not line:
				return  # StopIteration
			line = None

		assert False, "Should not reach this line"
