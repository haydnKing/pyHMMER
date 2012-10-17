from pyHMMER import hmmfile
import unittest

class TestHMMRead(unittest.TestCase):

	hmm = hmmfile.HMM('tests/data/valid.hmm')
	
	def test_valid(self):
		self.assertTrue(self.hmm.isValid())

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
		self.assertEqual(len(self.hmm.warnings), 2)









		
