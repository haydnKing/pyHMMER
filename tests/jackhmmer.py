import unittest

from pyHMMER import HMMER, matchfile
from Bio import SeqIO, Alphabet

class jackHMMER(unittest.TestCase):

	def test_jackhmmer(self):
		seq = SeqIO.read('tests/data/jackhmmer_seq.fasta', 'fasta', 
				alphabet=Alphabet.generic_protein)
		seqdb = SeqIO.read('tests/data/matchtarget.fasta', 'fasta',
				alphabet=Alphabet.generic_protein)

		j = HMMER.jackhmmer(seq, seqdb)

		#load the expected output
		seq_ = HMMER.wrap_seqrecords([seq,])
		seqdb_ = HMMER.wrap_seqrecords([seqdb,])
		m = matchfile.load('tests/data/jack_out', seq_, seqdb_)

		self.assertEqual(m, j.matches)

	def test_jackhmmer_dna(self):
		seq = SeqIO.read('tests/data/jackhmmer_seq.fasta', 'fasta', 
				alphabet=Alphabet.generic_protein)
		seqdb = SeqIO.read('tests/data/dna_target.fasta', 'fasta',
				alphabet=Alphabet.generic_dna)

		j = HMMER.jackhmmer(seq, seqdb)


		#load the expected output
		seqdb_prot = SeqIO.read('tests/data/matchtarget.fasta', 'fasta',
				alphabet=Alphabet.generic_protein)
		seq_ = HMMER.wrap_seqrecords([seq,])
		seqdb_ = HMMER.wrap_seqrecords([seqdb_prot,])
		matches = matchfile.load('tests/data/jack_out', seq_, seqdb_)
		#scale the matches' locations to match the protein search
		for m in matches:
			m.scale(3)

		self.assertEqual([str(m) for m in matches], [str(m) for m in j.matches])

	def test_arguments(self):
		seq = SeqIO.read('tests/data/jackhmmer_seq.fasta', 'fasta',
				alphabet=Alphabet.generic_protein)
		seqdb = SeqIO.read('tests/data/matchtarget.fasta', 'fasta',
				alphabet=Alphabet.generic_protein)

		j = HMMER.jackhmmer(seq, seqdb)
		args = j.getArgs(max=True, E='something')

		self.assertEqual(args, ['--max', '-E', 'something'])
		
