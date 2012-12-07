import unittest, pickle

from pyHMMER import HMMER, matchfile
from Bio import SeqIO

class jackHMMER(unittest.TestCase):

	def test_jackhmmer(self):
		seq = SeqIO.read('tests/data/jackhmmer_seq.fasta', 'fasta')
		seqdb = SeqIO.read('tests/data/matchtarget.fasta', 'fasta')

		j = HMMER.jackhmmer(seq, seqdb)

		#load the expected output
		m = pickle.load(open('tests/data/jack_out.pickle'))

		self.assertEqual(m, j.matches)
