from pyHMMER import hmmfile
import unittest

class TestHMMRead(unittest.TestCase):

	parser = hmmfile.HMMParser()
	hmm = parser.read('tests/data/valid.hmm')[0]

	def test_name(self):
		self.assertEqual(self.hmm.NAME, 'PPR_1')

	def test_acc(self):
		self.assertEqual(self.hmm.ACC, 'PF12854.2')

	def test_desc(self):
		self.assertEqual(self.hmm.DESC, 'PPR repeat')

	def test_leng(self):
		self.assertEqual(self.hmm.LENG, 34)

	def test_alph(self):
		self.assertEqual(self.hmm.ALPH.lower(), 'amino')

	def test_rf(self):
		self.assertFalse(self.hmm.RF)

	def test_cs(self):
		self.assertFalse(self.hmm.CS)

	def test_map(self):
		self.assertTrue(self.hmm.MAP)

	def test_date(self):
		self.assertEqual(self.hmm.DATE, 'Tue Sep 27 21:08:33 2011')

	def test_nseq(self):
		self.assertEqual(self.hmm.NSEQ, 433)

	def test_effn(self):
		self.assertEqual(self.hmm.EFFN, 54.372765)

	def test_cksum(self):
		self.assertEqual(self.hmm.CKSUM, 570345876)

	def test_ga(self):
		self.assertEqual(self.hmm.GA, (23.8, 23.9))

	def test_tc(self):
		self.assertEqual(self.hmm.TC, (23.8, 23.7))

	def test_nc(self):
		self.assertEqual(self.hmm.NC, (23.7, 23.6))

	def test_stats(self):
		self.assertEqual(self.hmm.STATS, [
			('LOCAL', 'MSV', -7.2358,  0.71971,),
			('LOCAL', 'VITERBI', -7.4863, 0.71971,),
			('LOCAL', 'FORWARD', -4.4686, 0.71971,),
			])

	def test_warnings(self):
		parser = hmmfile.HMMParser()
		hmm = parser.read('tests/data/valid.hmm')
		self.assertEqual(len(parser.warnings), 2)

	def test_invalid_line(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/invalid_line.hmm'))
	
	def test_invalid_leng(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/invalid_leng.hmm'))

	def test_invalid_stats(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/invalid_stats.hmm'))
	
	def test_invalid_alph(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/invalid_alph.hmm'))

	def test_invalid_stat_type(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/invalid_stat_type.hmm'))
	
	def test_no_head(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/no_head.hmm'))
	
	def test_no_hmm(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/no_hmm.hmm'))
						
	def test_no_leng(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/no_leng.hmm'))
	
	def test_no_name(self):
		self.assertRaises(hmmfile.HMMFileException, self.parser.read, 
				('tests/data/no_name.hmm'))

	def test_state(self):
		s = self.hmm.states[33]
		self.assertEqual(s.num, 33)
		self.assertEqual(s.me,[3.79753,7.35139,8.03144,7.44063,
					4.52454,7.32717,7.48394,3.90706,7.25453,3.18850,0.14720,7.47555,
					7.49919,7.18151,7.18775,6.67656,6.28879,4.64358,4.98721,4.09848,])
		self.assertEqual(s.MAP, 33)
		self.assertEqual(s.RF, '-')
		self.assertEqual(s.CS, '-')
		self.assertEqual(s.ie,[2.68618,4.42225,2.77519,2.73123,
			3.46354,2.40513,3.72494,3.29354,2.67741,2.69355,4.24690,2.90347,
			2.73739,3.18146,2.89801,2.37887,2.77519,2.98518,4.58477,3.61503,])
		self.assertEqual(s.tr, 
			[0.00075,7.59383,8.31617,0.61958,0.77255,0.48576,0.95510,])

import tempfile

def save_load(hmm):
	f = tempfile.TemporaryFile('w+r')
	hmmfile.write(hmm, f)
	f.seek(0)
	parser = hmmfile.HMMParser()
	hmmout = parser.read(f)[0]
	f.close()
	return hmmout


class TestHMMReadWrite(unittest.TestCase):

	parser = hmmfile.HMMParser()
	hmm = parser.read('tests/data/valid.hmm')[0]

	def test_options(self):
		hmm = save_load(self.hmm)
		for o in hmmfile.OPTIONS:
			if hasattr(self.hmm, o):
				self.assertEqual(getattr(self.hmm, o), getattr(hmm, o))

	def test_model(self):
		hmm = save_load(self.hmm)
		if hasattr(self.hmm, 'COMPO'):
			self.assertEqual(self.hmm.COMPO, hmm.COMPO)

		for i,s in enumerate(self.hmm.states):
			s2 = hmm.states[i]
			self.assertEqual(s.me, s2.me)
			self.assertEqual(s.ie, s2.ie)
			self.assertEqual(s.tr, s2.tr)
			self.assertEqual(s.MAP, s2.MAP)
			self.assertEqual(s.RF, s2.RF)
			self.assertEqual(s.CS, s2.CS)


		
