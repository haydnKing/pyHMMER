from pyHMMER import hmmfile, matchfile, sequtils, HMMER
from Bio import SeqIO
import unittest, tempfile, re




from pyHMMER import HMMER

from matchfile import check_valid, hmm_seq, ali_seq

class Testhmmsearch(unittest.TestCase):

	def test_search(self):
		h = HMMER.hmmsearch('tests/data/valid.hmm', 'tests/data/matchtarget.fasta')
		self.assertEqual(len(h.matches), 17)
		check_valid(self, h.matches)

	def test_translation_search(self):
		h = HMMER.hmmsearch('tests/data/valid.hmm', 
				'tests/data/reversedmatchtarget.fasta')
		self.assertEqual(len(h.matches), 17)
		s = sequtils.translate(str(h.matches[0].getSequence()))
		self.assertEqual(s, hmm_seq)

class TestFiltering(unittest.TestCase):

	hmm = hmmfile.read('tests/data/valid.hmm')[0]
	
	target = []
	for r in SeqIO.parse('tests/data/matchtarget.fasta', 'fasta'):
		target.append(r)
	
	target = target[0]

	def setUp(self):
		#build a fake search
		self.hs = HMMER.hmmsearch()
		self.hs.hmms = [self.hmm,]
		self.hs.targets = [self.target,]

	def build_match(self, start, stop, score, frame=0):
		"""Build a match"""
		m = matchfile.Match()
		m.hmm_from = 0
		m.hmm_to = self.hmm.leng
		m.ali_from = m.env_from = start
		m.ali_to   = m.env_to   = stop
		m.score = score
		m.target = self.target
		m.query = self.hmm
		m.frame = frame
		return m

	def build_matches(self, spec):
		self.hs.matches = []
		for m in spec:
			self.hs.matches.append( self.build_match(*m))
	
	def assertMatches(self, spec):
		"""
			assert that the matches remaining correspond to spec and are in that
			order
		"""
		self.assertEqual(len(self.hs.matches), len(spec))
		
		for i,m in enumerate(self.hs.matches):
			s = spec[i]
			estr = "{0}: {1} != ({2.env_from}, {2.env_to}, {2.score})".format(i,s,m)
			self.assertEqual(s[0], m.env_from, estr)
			self.assertEqual(s[1], m.env_to, estr)
			self.assertEqual(s[2], m.score, estr)
			try:
				self.assertEqual(s[3], m.frame)
			except IndexError:
				pass

	def test_min_1(self):
		"""
			--<1>--|--<1>--|--<1>--  (d=0)  =>  --<1>--|--<1>--|--<1>--
		"""
		self.build_matches(((50, 60, 1),
												(60, 70, 1),
												(70, 80, 1),))
		self.hs.mindist(0, 'env')
		self.assertMatches(((50, 60, 1),
												(60, 70, 1),
												(70, 80, 1),))

	def test_min_2(self):
		"""
			--<1>--|--<2>--|--<3>--  (d=1)  =>  --<1>--|       |--<3>--
		"""
		self.build_matches(((50, 60, 1),
												(60, 70, 2),
												(70, 80, 3),))
		self.hs.mindist(1, 'env', True)
		self.assertMatches(((50, 60, 1),
												(70, 80, 3),))	
	
	def test_min_3(self):
		"""
			--<1>--|       |--<3>--  (d=0)  =>    --<2>--      --<3>--
			    --<2>--
		"""
		self.build_matches(((50, 60, 1),
												(55, 65, 2),
												(70, 80, 3),))
		self.hs.mindist(0, 'env')
		self.assertMatches(((55, 65, 2),
												(70, 80, 3),))	
		
	def test_min_4(self):
		"""
			--<2>--|--<2>--|  (d=0)  =>  --<2>--|--<2>--
					--<1>--|--<1>--
		"""
		self.build_matches(((50, 60, 2),
												(60, 70, 2),
												(55, 65, 1),
												(65, 75, 1),))
		self.hs.mindist(0, 'env')
		self.assertMatches(((50, 60, 2),
												(60, 70, 2),))
	
	def test_min_5(self):
		"""
			---------<1>--------------  (d=0)  =>  --<2>--|--<2>--|--<2>--
				--<2>--|--<2>--|--<2>--
		"""
		self.build_matches(((40,100, 1),
												(50, 60, 2),
												(60, 70, 2),
												(80, 90, 2),))
		self.hs.mindist(0, 'env')
		self.assertMatches(((50, 60, 2),
												(60, 70, 2),
												(80, 90, 2),))

	def test_min_4_frames(self):
		""" 
		Test 4, but with frames
			--<2>--|--<2>--|  (d=0)  =>  --<2>--|--<2>--
					--<1>--|--<1>--              --<1>--|--<1>--
		"""
		self.build_matches(((50, 60, 2, 1),
												(60, 70, 2, 1),
												(55, 65, 1, 2),
												(65, 75, 1, 2),))
		self.hs.mindist(0, 'env')
		self.assertMatches(((50, 60, 2, 1),
												(60, 70, 2, 1),
												(55, 65, 1, 2),
												(65, 75, 1, 2),))

	def test_max_1(self):
		"""
			--<1>--|       |--<2>--  (d=2)  =>
		"""
		self.build_matches(((50, 60, 1),
												(70, 80, 2),))
		self.hs.maxdist(2, 'env')
		self.assertMatches(tuple())

	def test_max_2(self):
		"""
			--<1>--|--<1>--               --<1>--  (d=5)  =>  --<1>--|--<1>--
		"""
		self.build_matches(((50, 60, 1),
												(60, 70, 1),
												(100,120, 1),))
		self.hs.maxdist(5, 'env')
		self.assertMatches(((50, 60, 1),
												(60, 70, 1),))

	def test_max_3(self):
		"""
			--<1>--     ...     --<1>--  (d=5)  =>  --<1>--   ...   --<1>--
		"""
		self.build_matches(((5, 25, 1),
												(665,685, 1),))
		self.hs.maxdist(200, 'env')
		self.assertMatches(((5, 25, 1),
												(665,685, 1),))
