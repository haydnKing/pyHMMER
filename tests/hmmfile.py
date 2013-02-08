from pyHMMER import hmmfile
from StringIO import StringIO
import unittest

class TestHMMRead(unittest.TestCase):

	parser = hmmfile.HMMParser()
	hmm = parser.read('tests/data/valid.hmm')[0]

	def test_name(self):
		self.assertEqual(self.hmm.name, 'PPR_1')

	def test_acc(self):
		self.assertEqual(self.hmm.acc, 'PF12854.2')

	def test_desc(self):
		self.assertEqual(self.hmm.desc, 'PPR repeat')

	def test_leng(self):
		self.assertEqual(self.hmm.leng, 34)

	def test_alph(self):
		self.assertEqual(self.hmm.alph.lower(), 'amino')

	def test_rf(self):
		self.assertFalse(self.hmm.rf)

	def test_cs(self):
		self.assertFalse(self.hmm.cs)

	def test_map(self):
		self.assertTrue(self.hmm.map)

	def test_date(self):
		self.assertEqual(self.hmm.date, 'Tue Sep 27 21:08:33 2011')

	def test_nseq(self):
		self.assertEqual(self.hmm.nseq, 433)

	def test_effn(self):
		self.assertEqual(self.hmm.effn, 54.372765)

	def test_cksum(self):
		self.assertEqual(self.hmm.cksum, 570345876)

	def test_ga(self):
		self.assertEqual(self.hmm.ga, (23.8, 23.9))

	def test_tc(self):
		self.assertEqual(self.hmm.tc, (23.8, 23.7))

	def test_nc(self):
		self.assertEqual(self.hmm.nc, (23.7, 23.6))

	def test_stats(self):
		self.assertEqual(self.hmm.stats, [
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
		self.assertEqual(s.map, 33)
		self.assertEqual(s.rf, '-')
		self.assertEqual(s.cs, '-')
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
		if hasattr(self.hmm, 'compo'):
			self.assertEqual(self.hmm.compo, hmm.compo)

		for i,s in enumerate(self.hmm.states):
			s2 = hmm.states[i]
			self.assertEqual(s.me, s2.me)
			self.assertEqual(s.ie, s2.ie)
			self.assertEqual(s.tr, s2.tr)
			self.assertEqual(s.map, s2.map)
			self.assertEqual(s.rf, s2.rf)
			self.assertEqual(s.cs, s2.cs)

build_data = [{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':4},
		"emission":{'A': 9, 'C': 2, 'T': 4, 'G': 27},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 1, 'C': 1, 'T': 1, 'G': 1},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 0, 'C': 0, 'T': 0, 'G': 3},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 1, 'C': 1, 'T': 1, 'G': 1},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 2, 'C': 0, 'T': 0, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 7, 'C': 22, 'T': 68, 'G': 7},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 6, 'C': 46, 'T': 31, 'G': 3},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 6, 'C': 46, 'T': 31, 'G': 3},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 6, 'C': 46, 'T': 31, 'G': 3},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 1, 'C': 12, 'T': 3, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 0, 'C': 5, 'T': 3, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 6, 'C': 46, 'T': 31, 'G': 3},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 1, 'C': 1, 'T': 1, 'G': 1},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{		
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 0, 'C': 0, 'T': 3, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 7, 'C': 22, 'T': 68, 'G': 7},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 1, 'C': 1, 'T': 1, 'G': 1},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 0, 'C': 5, 'T': 3, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
	{
		"transition":{'MD': 0, 'DM': 1, 'MM': 7, 'DD': 0, 'MI': 1, 'II': 6, 'IM':
			4},
		"emission":{'A': 0, 'C': 0, 'T': 2, 'G': 0},
		"insert_emission":None,
		"MAP":-1,
		"CS":"-",
		"RF":"-"
		},
]

class TestHMMBuild(unittest.TestCase):

	em = {'a':2, 'C':1, 'D':1,}
	iem = {'D':2, 'E':1, 'F':1,}
	tr = {'MM':1,'MI':2,'MD':1,'IM':2,'II':8,'DM':1,'DD':1,}

	e_em = [0.5,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]
	e_iem= [0,0,0.5,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]

	e_tr = [[0.25,0.5,0.25,0.2,0.8,1,0,],
					[0.25,0.5,0.25,0.2,0.8,0.5,0.5,],
					[1.0/3,2.0/3,0,0.2,0.8,1,0,],]
	
	def test_build(self):
		hmm = hmmfile.HMM(alphabet='DNA')

		#build the HMM
		for i in build_data:
			hmm.addState(**i)
		hmm.clean()

		#there are the correct number of states
		self.assertEqual(len(hmm.states), len(build_data), 
			"{} states built - should be {}".format(len(hmm.states),
				len(build_data)))

		s = StringIO()

		hmmfile.write(hmm, s)

		s.seek(0)
		hmm2 = hmmfile.read(s)[0]

		for o in hmmfile.OPTIONS:
			if hasattr(hmm, o):
				self.assertEqual(getattr(hmm, o), getattr(hmm2, o))

		def r(l):
			ret = []
			for i in l:
				if type(i) == float:
					ret.append(round(i,5))
				else:
					ret.append(i)
			return ret

		for i,(e,g) in enumerate(zip(hmm.states, hmm2.states)):
			msg = "State {}: {} != {}"
			self.assertEqual(r(e.me), g.me, msg.format(i,e.me,g.me))
			self.assertEqual(r(e.ie), g.ie, msg.format(i,e.ie,g.ie))
			self.assertEqual(r(e.tr), g.tr, msg.format(i,e.tr,g.tr))
			self.assertEqual(e.map,g.map, msg.format(i,e.map,g.map))
			self.assertEqual(e.rf, g.rf, msg.format(i,e.rf,g.rf))
			self.assertEqual(e.cs, g.cs, msg.format(i,e.cs,g.cs))


	def test_blank(self):
		hmm = hmmfile.HMM(alphabet = 'DNA')

		for i in range(5):
			hmm.addState()
		hmm.clean()

		self.assertEqual(len(hmm.states), 5)
		#inner states should be equal, except for num
		for i in range(2,4):
			hmm.states[i].num = 1
			self.assertEqual(hmm.states[1], hmm.states[i])

	def test_duplicates(self):
		hmm = hmmfile.HMM(alphabet='AMINO')

		self.assertRaises(ValueError, hmm.addState, 
				transition={'MM':1,'MI':2,'MD':1,'IM':2,'II':8,'DM':1,'DD':1,'mm':0,},
				emission=self.em, 
				insert_emission=self.iem)
		
		self.assertRaises(ValueError, hmm.addState, 
				transition=self.tr, 
				emission={'a':2, 'C':1, 'D':1,'A':2}, 
				insert_emission=self.iem)
		
		self.assertRaises(ValueError, hmm.addState, 
				transition=self.tr, 
				emission=self.em, 
				insert_emission={'D':2, 'E':1, 'F':1, 'd':2,})



