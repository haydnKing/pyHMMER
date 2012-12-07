import unittest

from pyHMMER import HMMER, matchfile
from Bio import SeqIO

class jackHMMER(unittest.TestCase):

	def test_jackhmmer(self):
		seq = SeqIO.read('tests/data/jackhmmer_seq.fasta', 'fasta')
		seqdb = SeqIO.read('tests/data/matchtarget.fasta', 'fasta')

		j = HMMER.jackhmmer(seq, seqdb)

		#load the expected output
		seq_ = HMMER.wrap_seqrecords([seq,])
		seqdb_ = HMMER.wrap_seqrecords([seqdb,])
		m = matchfile.load('tests/data/jack_out', seq_, seqdb_)

		self.assertEqual(m, j.matches)
