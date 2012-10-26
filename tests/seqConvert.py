from pyHMMER import seqConvert
import unittest

test = ("TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCC" +
 "CACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGA" +
 "AGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG")

six_frame = (
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	"FSYCLPHRITNSVADGFSYCLPHRITNSVADGLS**LPQRITKRVAEGLS*WLPQRMTKRVAE",
	"FLIVFLIVLLIVLLMASPTASPTASPTASPTDYHNDYHNE*QKE*QKGCRSGCRSG*RRGWRR",
	"PSATLFVIRCGSHYDNPSATLFVIRCGSHYDNPSATLLVMRWGRQ*EKPSATLLVIR*GRQ*EK",
	"PPPPSSSSAAAATTTTLLLLFLLFVVVVIMIIRRRRCW*CGGGGSRRSHQQHY**YDEEDNKK",
	"LRHPLRHPLRQPLRQPFCYSFCYSLW*SL**SVGDAVGDAVGEAVGEAISNTISNTMRKTIRK",
)

class TestConversion(unittest.TestCase):

	def test_sixFrame(self):
		c = seqConvert.sixFrameTranslation(test)
		self.assertEqual(c, six_frame)

	
