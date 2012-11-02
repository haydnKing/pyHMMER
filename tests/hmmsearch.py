from pyHMMER import hmmfile, matchfile, sequtils
from Bio import SeqIO
import unittest, tempfile, re

test_data = {
	't_name': ['gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|','gi|157931526|gb|ABW04887.1|',],
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
	'desc': ['PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]','PPR [Raphanus sativus]',],
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
		targets.append(r)
	
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
