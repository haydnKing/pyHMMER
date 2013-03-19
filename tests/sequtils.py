from pyHMMER import sequtils
import unittest

test = ("TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCC" +
 "CACCGCATCACCAACAGCGTCGCCGACGGATTATcATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGA" +
 "AGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG")

translation ="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

six_frame = (
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	"FSYCLPHRITNSVADGFSYCLPHRITNSVADGLS**LPQRITKRVAEGLS*WLPQRMTKRVAE",
	"FLIVFLIVLLIVLLMASPTASPTASPTASPTDYHNDYHNE*QKE*QKGCRSGCRSG*RRGWRR",
	"PSATLFVIRCGSHYDNPSATLFVIRCGSHYDNPSATLLVMRWGRQ*EKPSATLLVIR*GRQ*EK",
	"PPPPSSSSAAAATTTTLLLLFLLFVVVVIMIIRRRRCW*CGGGGSRRSHQQHY**YDEEDNKK",
	"LRHPLRHPLRQPLRQPFCYSFCYSLW*SL**SVGDAVGDAVGEAVGEAISNTISNTMRKTIRK",
)


class TestConversion(unittest.TestCase):

	def test_translation(self):
		s = sequtils.translate(test)
		self.assertEqual(s,translation)

	def test_sixFrame(self):
		c = sequtils.sixFrameTranslation(test)
		self.assertEqual(tuple(c), six_frame)

	def test_minLength(self):
		self.assertRaises(ValueError, sequtils.translate, 'aa')
		self.assertRaises(ValueError, sequtils.sixFrameTranslation, 'aa')

	def test_invalidChar(self):
		self.assertRaises(ValueError, sequtils.translate, 'aaaaaf')
		c = sequtils.sixFrameTranslation('aaaaaaaaaaajk')
		self.assertRaises(ValueError, c.next)

	
class TestItentification(unittest.TestCase):

	def test_dna(self):
		self.assertEqual(sequtils.seq_type('ATGCAtGA'), 'DNA')

	def test_rna(self):
		self.assertEqual(sequtils.seq_type('AUGCAuGA'), 'RNA')

	def test_amino(self):
		self.assertEqual(sequtils.seq_type('LEDaIDLFSD'), 'AMINO')
