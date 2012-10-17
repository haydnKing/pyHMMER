"""Read and Write HMMER's .hmm format"""

import re

class HMMFileException(Exception):
	def __init__(self, errors, warnings):
		self.errors = errors
		self.warnings = warnings

	def __unicode__(self):
		ret = ''
		if(len(self.errors) > 0):
			if(len(self.errors) == 1):
				ret += 'Error:\n'
			else:
				ret += 'Errors:\n'
			for e in self.errors:
				ret += '\t%s\n' % e

		if(len(self.warnings) > 0):
			if(len(self.warnings) == 1):
				ret += 'Warning:\n'
			else:
				ret += 'Warnings:\n'
			for w in self.warnings:
				ret += '\t%s\n' % w
		return ret

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

ALPHABETS = {
		'DNA': 'ACGT',
		'RNA': 'ACGU',
		'AMINO': 'ACDEFGHIKLMNPQRSTVWY',
		}

class HMM:
	"""A Hidden Markov Model"""
	def __init__(self, fname):
		#open the file and read the file
		f = open(fname, 'r')
		lines = list(enumerate([line.strip() for line in f]))
		f.close()
		
		#first line must start with HMMER3/b
		if not re.match(r"^HMMER3/b", lines[0][1]):
			raise HMMFileException(["Invalid File: file must start with "
				"\'HMMER3/b\'"], [])

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
			raise HMMFileException(['Invalid File: no \'HMM\' line found',], [])
		
		self.errors = []
		self.warnings = []

		#parse the header section
		self._read_hdr(header)
		
		#parse the model section
		#self._read_mdl(model)

		if not self.isValid():
			raise HMMFileException(self.errors, self.warnings)
		
	def getWarnings(self):
		ret = ''
		if len(self.warnings):
			if len(self.warnings) == 1:
				ret += "Warning parsing file:"
			else:
				ret += "Warnings parsing file:"
			for error in self.errors:
				if error[0] > 0:
					ret += "\tLine %d: %s" % (error[0], error[1])
				else:
					ret += "\t%s" % error[1]
		return ret

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
				self._addError(line[0], 'invalid line')
				continue
			key = m.group('key').upper()
			val = m.group('value')
			if not (key in OPTIONS):
				self._addWarning(line[0], 'ignoring unknown option \'%s\'' % key)
				continue

			#simple strings
			if key in ['NAME', 'ACC', 'DESC', 'ALPH', 'DATE',]:
				setattr(self, key, val)
				if key == 'ALPH':
					if self.ALPH.upper() in ALPHABETS:
						self.SYMBOLS = ALPHABETS[self.ALPH.upper()]
						self.K = len(self.SYMBOLS)
					else:
						self._addError(line[0], 'ALPH must be \'DNA\', \'RNA\' or \'AMINO\'')
			#integers
			elif key in ['LENG', 'NSEQ', 'CKSUM',]:
				try:
					val = int(val)
					if val < 0:
						raise ValueError
					setattr(self, key, val)
				except ValueError:
					self._addError(line[0], 'Must be a positive integer')
					continue
			#bools
			elif key in ['RF', 'CS', 'MAP',]:
				b = {'yes': True, 'no': False,};
				val = val.lower()
				if val in b:
					setattr(self, key, b[val])
				else:
					self._addError(line[0], 'Line %s: value must be \'yes\' or \'no\'')
					continue
			#COM
			elif key == 'COM':
				m2 = re.match(r'(\d+)\s+(\S+)$', val)
				if m2:
					try:
						self.COM.append( (int(m2.group(1)), m2.group(2)))
					except ValueError:
						self.COM.append( (len(self.COM)+1, m2.group(0)))
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
					self._addError(line[0], 'EFFN must be a positive real')
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
						self._addError(line[0], '%s must be two positive reals' % key)
						continue
				else:
					self._addError(line[0], '%s must be two positive reals' % key)
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
						else:
							if s1 != 'LOCAL':
								self._addError(line[0], 's1 must equal \'LOCAL\'')
							else:
								self._addError(line[0], 
										's2 must be \'MSV\', \'VITERBI\' or \'FORWARD\'')
					except ValueError:
						self._addError(line[0], 'STATS <s1> <s2> <f1> <f2>')
						continue
				else:
					self._addError(line[0], 'STATS <s1> <s2> <f1> <f2>')
					continue
		
		#Parsed each line
		#check all the required options were present
		for o in REQUIRED:
			try:
				getattr(self, o)
			except AttributeError:
				self._addError(-1, 'Option \'%s\' is required')


									
	def _addError(self, line, what):
		self.errors.append((line, what))

	def _addWarning(self, line, what):
		self.warnings.append((line, what))


