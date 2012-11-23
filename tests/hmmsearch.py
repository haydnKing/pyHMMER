from pyHMMER import hmmfile, matchfile, sequtils
from Bio import SeqIO
import unittest, tempfile, re

test_data = {
	't_name': ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0',],
	't_accession': ['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-',], 
	't_len': [687,687,687,687,687,687,687,687,687,687,687,687,687,687,687,687,687,],
	'q_name': ['PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1','PPR_1',],
	'q_accession': ['PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2','PF12854.2',], 
	'q_len': [34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,], 
	'f_evalue': [4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,4.2e-183,], 
	'f_score': [579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,579.0,], 
	'f_bias': [31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,31.3,], 
	'num': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,], 
	'of': [17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,],
	'c_evalue': [  0.017, 3.3e-06, 1.3e-15, 6.2e-18, 7.1e-13, 6.5e-14, 1.6e-15,
		2.2e-14, 1.7e-14, 1.2e-17, 1.9e-11, 4.7e-15, 2.6e-13, 4.2e-17, 1.4e-19,
		4.3e-09,	 0.051,],
	'i_evalue': [0.017,3.3e-06,1.3e-15,6.2e-18,7.1e-13,6.5e-14,1.6e-15,2.2e-14,1.7e-14,1.2e-17,1.9e-11,4.7e-15,2.6e-13,4.2e-17,1.4e-19,4.3e-09,0.051,],
	'score': [1.2,13.1,43.2,50.6,34.4,37.7,42.9,39.3,39.6,49.6,29.8,41.4,35.8,47.9,55.8,22.3,-0.3,],
	'bias': [0.1	,	0.1	,	0.0	,	0.3	,	0.0	,	0.0	,	0.0	,	0.0	,	0.0	,	0.1	,	0.0	,	0.2	,	0.0	,	0.0	,	0.1	,	0.0	,	0.0	,		],
	'hmm_from': [23, 3, 2, 3, 1, 3, 1, 2, 1, 1, 2, 2, 1, 1, 3, 1, 2,		], 
	'hmm_to': [33,30,33,33,34,34,33,33,31,34,33,34,34,33,33,33,32,		], 
	'ali_from': [60,110,144,180,213,251,284,320,354,389,425,460,505,540,577,610,646,		], 
	'ali_to': [70,137,175,210,246,282,316,351,384,422,456,492,538,572,607,642,676,], 
	'env_from': [60,108,143,178,213,249,284,319,354,389,424,459,505,540,575,610,645,		], 
	'env_to': [71,140,176,211,246,282,317,352,387,422,457,492,538,573,608,643,678,		],
	'acc': [0.92,0.84,0.98,0.97,0.97,0.95,0.98,0.96,0.92,0.98,0.97,0.96,0.97,0.98,0.95,0.97,0.68,		], 
}

hmm_seq = "SCEAGFGGESLKLQSGFHEIKGLEDAIDLFSDML"
ali_seq = "LEDAIDLFSD"

def check_valid(self, matches):
	for i,m in enumerate(matches):
		for k in test_data.iterkeys():
			if re.search(r"_to$|_from$", k):
				v = getattr(m,k) + 1
			else:
				v = getattr(m,k)
			self.assertEqual(v, test_data[k][i], msg=
				"matches[{}].{} = {} != {}".format(i,k,getattr(m,k),test_data[k][i]))


class TestMatchRead(unittest.TestCase):

	hmms = hmmfile.read('tests/data/valid.hmm')

	targets = list()
	for r in SeqIO.parse('tests/data/matchtarget.fasta', 'fasta'):
		targets.append((r.id,r,))
	
	matches = matchfile.load('tests/data/matchfile', hmms, targets)

	def test_length(self):
		self.assertEqual(len(self.matches), 17)

	def test_values(self):
		check_valid(self, self.matches)	

	def test_writeread(self):
		f = tempfile.TemporaryFile('w+r')
		matchfile.save(self.matches, f)

		f.seek(0)

		m2 = matchfile.load(f, self.hmms, self.targets)

		check_valid(self, m2)

		f.close()

	def test_aliseq(self):
		self.assertEqual(str(self.matches[0].getSequence('ali')), ali_seq)

	def test_hmmseq(self):
		self.assertEqual(str(self.matches[0].getSequence('hmm')), hmm_seq)

from pyHMMER import HMMER

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
		m.hmm_to = self.hmm.LENG
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
