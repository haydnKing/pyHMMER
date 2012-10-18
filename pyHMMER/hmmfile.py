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
				if e[0] >= 0:
					ret += '\tLine %s: %s\n' % (e[0]+1, e[1])
				else:
					ret += '\t%s\n' % e[1]

		if(len(self.warnings) > 0):
			if(len(self.warnings) == 1):
				ret += 'Warning:\n'
			else:
				ret += 'Warnings:\n'
			for w in self.warnings:
				if w[0] >= 0:
					ret += '\tLine %s: %s\n' % (w[0]+1, w[1])
				else:
					ret += '\t%s\n' % w[1]
		return ret

	def __str__(self):
		return unicode(self).encode('utf-8')

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

class State:
	"""A state within an HMM"""
	def __init__(self, num, MAP_annot, CS_annot, RF_annot, 
						match_emission, insert_emission, transition):
		self.num = num
		self.MAP_annot = MAP_annot
		self.RF_annot = RF_annot
		self.CS_annot = CS_annot
		self.match_emission = match_emission
		self.insert_emission = insert_emission
		self.transition = transition

class HMM:
	"""A Hidden Markov Model"""

	def __init__(self, fname=''):
		self.errors = []
		self.warnings = []
		self.line = -1

		self.states = []
		self.COM = []
		self.STATS = []
		#assume these to be false if not present
		self.RM = self.CS = self.MAP = False
		self.K = 0

		if fname:
			self.loadFile(fname)
	
	def loadFile(self, fname):
		"""load a file from disk"""
		#open the file and read the file
		f = open(fname, 'r')
		self.load(f)
		f.close()

	def load(self, f):
		"""load from a file object"""
		lines = [line.strip() for line in f]

		#Should I try to read the model
		self._continue = True			

		#parse the header section
		self._read_hdr(lines)
		
		if self._continue:
			#parse the model section
			self._read_mdl(lines)

		if not self.isValid():
			raise HMMFileException(self.errors, self.warnings)
		
	def isValid(self):
		return (len(self.errors) == 0)

	def writeToFile(self, fname):
		f = open(fname, 'w')
		self.write(f)
		f.close()

	def write(self, f):
		f.writelines(self.getLines())

	def getLines(self):
		ret = [];
		ret.append('HMMER3/b [pyHMMER | 2012]\n')
		#write header
		for o in ['NAME', 'ACC', 'DESC', 'LENG', 'ALPH', 'DATE', 'NSEQ', 'EFFN',
				'CKSUM',]:
			try:
				ret.append('{:<5s} {}\n'.format(o, getattr(self, o)))
			except AttributeError:
				pass
		for o in ['RF', 'CS', 'MAP']:
			try:
				if getattr(self, o):
					ret.append('{:<5s} yes\n'.format(o))
				else:
					ret.append('{:<5s} no\n'.format(o))
			except AttributeError:
				pass

		for c in self.COM:
			ret.append('COM   %s %s\n' % c)

		for o in ['GA', 'TC', 'NC',]:
			try:
				out = getattr(self, o)
				ret.append('{:<5s} {:.2f}  {:.2f}\n'.format(o, out[0], out[1]))
			except AttributeError:
				pass

		for s in self.STATS:
			ret.append('STATS {:<5s} {:<10s} {:8f} {:8f}\n'.format(*s))

		#write HMM lines
		ret.append('HMM ')
		ret.append( '%9s'*len(self.SYMBOLS) % tuple(self.SYMBOLS))
		ret.append(
			'\n            m->m     m->i     m->d'
			'     i->m     i->i     d->m     d->d\n')

		def ff(f, l=9, p=5):
			try:
				return ('%'+str(l)+'.'+str(p)+'f') % f
			except TypeError:
				return ' '*(l-1) + '*'
		#write CMPO, if it exists. Should probably calculate it...
		try:
			ret.append(("%7s" + ("%9s" * self.K) + "\n") % 
					(('COMPO',) + tuple((ff(f) for f in self.COMPO))))
		except AttributeError:
			pass

		#write model
		me_fmt = "%7i" + ("%9s" * self.K) + "%7i %1s %1s\n"
		ie_fmt = (" " * 7) + ("%9s" * self.K) + "\n"
		tr_fmt = (" " * 7) + ("%9s" * 7) + "\n"

		#0th state gets special treatment
		if len(self.states) < 1:
			return ret

		ret.append(ie_fmt % tuple((ff(f) for f in self.states[0].insert_emission)))
		ret.append(tr_fmt % tuple((ff(f) for f in self.states[0].transition)))

		#write model
		for i,s in enumerate(self.states[1:], 1):
			ret.append(me_fmt % ((i,) + tuple((ff(f) for f in s.match_emission)) + 
					(s.MAP_annot, s.RF_annot, s.CS_annot)))
			ret.append(ie_fmt % tuple((ff(f) for f in s.insert_emission)))
			ret.append(tr_fmt % tuple((ff(f) for f in s.transition)))

		#end with //
		ret.append('      //\n')
		return ret

	def _read_hdr(self, lines):
		"""parse the header"""


		#have we see the HMM line?
		seen_HMM = False

		#look for the format specification
		#first line must start with HMMER3/b
		self._next_line(lines);
		if not re.match(r"^HMMER3/b", lines[self.line]):
			raise HMMFileException(["Invalid File: file must start with "
				"\'HMMER3/b\'"], [])

		r = re.compile(r'^(?P<key>\w+)\s+(?P<value>.+)$')
		while self._next_line(lines) < len(lines):
			if re.match(r'^HMM\s', lines[self.line]):
				seen_HMM = True
				break
			
			m = r.match(lines[self.line])
			if not m:
				self._addError('invalid line')
				continue
			
			key = m.group('key').upper()
			val = m.group('value')
			#check if the key is valid
			if not (key in OPTIONS):
				self._addWarning('ignoring unknown option \'%s\'' % key)
				continue

			#simple strings
			if key in ['NAME', 'ACC', 'DESC', 'ALPH', 'DATE',]:
				setattr(self, key, val)
				if key == 'ALPH':
					if self.ALPH.upper() in ALPHABETS:
						self.SYMBOLS = ALPHABETS[self.ALPH.upper()]
						self.K = len(self.SYMBOLS)
					else:
						self._addError('ALPH must be \'DNA\', \'RNA\' or \'AMINO\'')
			#integers
			elif key in ['LENG', 'NSEQ', 'CKSUM',]:
				try:
					val = int(val)
					if val < 0:
						raise ValueError
					setattr(self, key, val)
				except ValueError:
					self._addError('Must be a positive integer')
					continue
			#bools
			elif key in ['RF', 'CS', 'MAP',]:
				b = {'yes': True, 'no': False,};
				val = val.lower()
				if val in b:
					setattr(self, key, b[val])
				else:
					self._addError('value must be \'yes\' or \'no\'')
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
					self._addError('EFFN must be a positive real')
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
						self._addError('%s must be two positive reals' % key)
						continue
				else:
					self._addError('%s must be two positive reals' % key)
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
								self._addError('s1 must equal \'LOCAL\'')
							else:
								self._addError('s2 must be \'MSV\', \'VITERBI\' or \'FORWARD\'')
					except ValueError:
						self._addError('STATS <s1> <s2> <f1> <f2>')
						continue
				else:
					self._addError('STATS <s1> <s2> <f1> <f2>')
					continue
		
		#Parsed each line
		#check all the required options were present
		for o in REQUIRED:
			try:
				getattr(self, o)
			except AttributeError:
				self._addError('Option \'%s\' is required' % o, False)
				if o == 'LENG':
					setattr(self, o, 0)
				else:
					setattr(self, o, '')

		if not seen_HMM:
			self._addError('No HMM line found', False)
			self._continue = False

	def _read_mdl(self, lines):
		"""read the model"""
		#ignore the line immediately after the HMM line
		self._next_line(lines)
		#parse the COMPO line, if it's present
		if re.match(r'^COMPO\s', lines[self.line+1]):
			self._next_line(lines)
			self.COMPO = self._parse_prob(lines[self.line].split()[1:], self.K)

		#variables to collect for each model
		num = 0
		MAP_annot = 0
		RF_annot = '-'
		CS_annot = '-'
		match_emission = []
		insert_emission = ['*' for i in range(0,self.K)]
		transition = []

		#state = 0 - expect Match Emission line or //
		#state = 1 - expect Insert Emission line
		#state = 2 - expect State Transision line
		#state = 4 - we've seen //
		state = 1
		while self._next_line(lines) < len(lines):
			if state == 0:
				if lines[self.line] == '//':
					state = 4
					break
				#parse the emission line
				l = lines[self.line].split()
				if len(l) != (self.K + 4):
					self._addError('Malformed Match Emission line')
					return
				try:
					num = int(l[0])
					if num != len(self.states):
						self._addError('Expected state number %s' % len(self.states))
				except ValueError:
					self._addError(self.parse_line, 'node number must be a positive integer')
				match_emission = self._parse_prob(l[1:self.K+1])
				#MAP number
				try:
					if(l[self.K+1] != '-'):
						MAP_annot = int(l[self.K+1])
					elif self.MAP:
						self._addWarning('Map annotation is \'-\', even though MAP is \'yes\'')
				except ValueError:
					self._addError('Map Annotation must be an integer or \'-\'')
				#RF annotation
				RF_annot = l[self.K+2]
				if len(RF_annot) != 1:
					self._addError('RF annotation must be a single character')
				#CS annotation
				CS_annot = l[self.K+3]
				if len(CS_annot) != 1:
					self._addError('CS annotation must be a single character')
				#we're now expecting an IE line
				state = 1
			
			elif state == 1:
				insert_emission = self._parse_prob(lines[self.line].split(), self.K)
				state = 2

			elif state == 2:
				transition = self._parse_prob(lines[self.line].split(), 7)
				state = 0
				#add the node
				self.states.append(State(num, MAP_annot, CS_annot, RF_annot, 
						match_emission, insert_emission, transition))

		if len(self.states) != (self.LENG + 1):
			self._addError('Expected %s states (incl. BEGIN) but found %s' % 
					(self.LENG+1, len(self.states)))
		if state != 4:
			self._addError('No \'//\' found at end of file')

	def _next_line(self, lines):
		self.line = self.line + 1
		try:
			while (len(lines[self.line].strip()) == 0 
					or lines[self.line].strip()[0] == '#'):
				self.line = self.line + 1
		except IndexError:
			pass

		return self.line

	def _parse_prob(self, l, expected=-1):
		ret = []
		for v in l:
			if str(v) == "*":
				ret.append(v)
			else:
				try:
					if float(v) < 0:
						self._addError('log probability was less than zero (%s)' % v)
						ret.append('*')
					else:
						ret.append(float(v))
				except ValueError:
					self._addError('probability should be a positive float (%s)' % v)
		if expected > 0 and len(ret) != expected:
			self._addError('expected %s floats' % expected)
		return ret

	def _addError(self, what, show_line=True):
		if(show_line):
			self.errors.append((self.line, what))
		else:
			self.errors.append((-1, what))

	def _addWarning(self, what, show_line=True):
		if(show_line):
			self.warnings.append((self.line, what))
		else:
			self.warnings.append((-1, what))


