from pyHMMER import seqConvert
import unittest

test = ("TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCC" +
 "CACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGA" +
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
		s = seqConvert.translate(test)
		self.assertEqual(s,translation)

	def test_sixFrame(self):
		c = seqConvert.sixFrameTranslation(test)
		self.assertEqual(c, six_frame)

	def test_minLength(self):
		self.assertRaises(ValueError, seqConvert.translate, 'aa')
		self.assertRaises(ValueError, seqConvert.sixFrameTranslation, 'aa')

	def test_invalidChar(self):
		self.assertRaises(ValueError, seqConvert.translate, 'aaaaaf')
		self.assertRaises(ValueError, seqConvert.sixFrameTranslation, 'aaaaaajk')

	
