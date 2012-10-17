"""Read and Write HMMER's .hmm format"""

import re

class HMMFileException(Exception):
	pass

OPTIONS = [	
						'NAME',
						'ACC',
						'DESC',
						'LENG',
						'ALPH',
						'RF',
						'CS',
						'MAP',
						'DATE',
						'COM',
						'NSEQ',
						'EFFN',
						'CKSUM',
						'GA',
						'TC',
						'NC',
						'STATS',
					]
REQUIRED = [
						'NAME',
						'LENG',
						'ALPH',
						]

class HMM:
	"""A Hidden Markov Model"""
	def __init__(self, fname):
		#open the file and read the file
		f = open(fname, 'r')
		lines = list(enumerate([line.strip() for line in f]))
		f.close()
		
		#first line must start with HMMER3/b
		if not re.match(r"^HMMER3/b", lines[0][1]):
			raise HMMFileException("Invalid File:"
					"file must start with	\'HMMER3/b\'")

		#split into header and model
		header = []
		model = []
		found = False
		for i, l in lines:
			if re.match(r'^HMM\s', l):
				header = lines[1:i]
				model = lines[i:]
				found = True
				break
		if not found:
			raise HMMFileException('Invalid File: no \'HMM\' line found')
		
		self.errors = []
		self.warnings = []

		#parse the header section
		self._read_hdr(header)
		
		#parse the model section
		#self._read_mdl(model)
		
	def printErrors(self):
		if len(self.errors):
			print "Error parsing file:"
			for error in self.errors:
				print "\t%s" % error

	def printWarnings(self):
		if len(self.warnings):
			print "Warning:"
			for warn in self.warnings:
				print "\t%s" % warn

	def isValid(self):
		return (len(self.errors) == 0)

	def _read_hdr(self, lines):
		"""parse the header"""
		self.COM = []
		self.STATS = []
		r = re.compile(r'^(?P<key>\w+)\s+(?P<value>.+)$')
		for line in lines:
			m = r.match(line[1])
			if not m:
				self.errors.append('Line %s: invalid line'% line[0] + 1)
				continue
			key = m.group('key').upper()
			val = m.group('value')
			if not (key in OPTIONS):
				self.warnings.append('Line %s: ignoring unknown option \'%s\'' % 
						(line[0] + 1,	key))
				continue

			#simple strings
			if key in ['NAME', 'ACC', 'DESC', 'ALPH', 'DATE',]:
				setattr(self, key, val)
			#integers
			elif key in ['LENG', 'NSEQ', 'CKSUM',]:
				try:
					val = int(val)
					if val < 0:
						raise ValueError
					setattr(self, key, val)
				except ValueError:
					self.errors.append('Line %s: Must be a positive integer')
					continue
			#bools
			elif key in ['RF', 'CS', 'MAP',]:
				b = {'yes': True, 'no': False,};
				val = val.lower()
				if val in b:
					setattr(self, key, b[val])
				else:
					self.errors.append('Line %s: value must be \'yes\' or \'no\'' % line[0]+1)
					continue
			#COM
			elif key == 'COM':
				m2 = re.match(r'(\d+)\s+(\S+)$', val)
				if m2:
					try:
						self.COM.append( (int(m2.group(1)), m2.group(2)))
					except ValueError:
						self.errors.append('Line %s: Invalid command number'% line[0]+1)
						continue
				else:
					self.COM.append( (len(self.COM)+1, value))
			elif key == 'EFFN':
				try:
					val = float(val)
					if val < 0:
						raise ValueError
					self.EFFN = val
				except ValueError:
					self.errors.append('Line %s: EFFN must be a positive real'% line[0]+1)
					continue
			#pairs of floats
			elif key in ['GA', 'TC', 'NC']:
				m2 = re.match(r'^([\d\.]+)\s+([\d\.]+);?$', val)
				if m2:
					try:
						v1 = float(m2.group(1))
						v2 = float(m2.group(2))
						if v1 < 0 or v2 < 0:
							raise ValueError
						setattr(self, key, (v1, v2))
					except ValueError:
						self.errors.append('Line %s: %s must be two positive reals' % 
								(line[0]+1, key))
						continue
				else:
					self.errors.append('Line %s: %s must be two positive reals' % 
								(line[0]+1, key))
			#STATS
			elif key == 'STATS':
				m2 = re.match(r'^(\w+)\s+(\w+)\s+([\d\.-]+)\s+([\d\.-]+)$', val)
				if m2:
					try:
						s1 = m2.group(1).upper()
						s2 = m2.group(2).upper()
						f1 = float(m2.group(3))
						f2 = float(m2.group(4))
						if s1 == 'LOCAL' and s2 in ['MSV', 'VITERBI', 'FORWARD']:
							self.STATS.append((s1, s2, f1, f2))
					except ValueError:
						self.errors.append('Line %s: STATS <s1> <s2> <f1> <f2>' % line[0]+1)
						continue
				else:
					self.errors.append('Line %s: STATS <s1> <s2> <f1> <f2>' % line[0]+1)
					continue
		
		#Parsed each line
		#check all the required options were present
		for o in REQUIRED:
			try:
				getattr(self, o)
			except AttributeError:
				self.errors.append('Option \'%s\' is required')
									


