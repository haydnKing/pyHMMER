"""Read and Write HMMER's .hmm format"""

import re, itertools, collections

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

class HMMParser:
	"""Parse a Hidden Markov Model in HMMER3's .hmm format"""

	#States
	S_default = 0
	S_header = 1
	S_model_me = 2
	S_model_ie = 3
	S_model_st = 4

	#useful regexes
	hdr_re = re.compile(r'^(?P<key>\w+)\s+(?P<value>.+)$')

	def __init__(self):
		self.lineNum = -1
		self.errors = []
		self.warnings = []

	def read(self, f):
		"""Read all the HMMs in the file, either a fileobject or a string"""

		#if f is a string, try to open the file
		if isinstance(f, basestring):
			f = open(f, 'r')
		
		ret = []
		self.errors = []
		self.warnings = []

		def _cleanLine(l):
			"""remove any comments and strip whitespace"""
			i = l.find('#')
			if i:
				l = l[:i]
			return l.strip()

		def consume(iterator, n):
			'''Advance the iterator n-steps ahead. If n is none, consume
					entirely.'''
			collections.deque(itertools.islice(iterator, n), maxlen=0)

		#read in all the lines, cleaning each one
		lines = [_cleanLine(l) for l in f]

		state = self.S_default

		iterator = enumerate(lines).__iter__()
		for self.lineNum, line in iterator:
			#ignore blank lines
			if not line:
				continue

			#find the start of an HMM definition
			if state == self.S_default:
				if not re.match(r"^HMMER3/b", line):
					self._addError("Invalid File: file must start with \'HMMER3/b\'")
					#can't continue
					self._raiseErrors()
				else:
					#make a new model
					hmm = HMM()
					state = self.S_header
					continue

			#have we found the end of the HMM
			if re.match(r'^//', line):
				#if this isn't where it's supposed to be
				if state != self.S_model_me:
					#provide a helpful error
					ltype = 'None'
					if state == self.S_header:
						ltype = 'header'
					elif state == self.S_model_ie:
						ltype = 'insert emission'
					elif state == self.S_model_st:
						ltype = 'state transition'
					self._addError('Unexpected \'//\', expecting {} line'.format(ltype))

				#check if we've found the advertised number of states 
				#			nb. hmm.LENG doesn't include the begin state """
				if len(hmm.states)-1 < hmm.LENG:
					self._addError('Too few states in model ({:d} < {:d})'
							.format(len(hmm.states)-1, hmm.LENG))
				if len(hmm.states)-1 > hmm.LENG:
					self._addError('Too many states in model ({:d} > {:d})'
							.format(len(hmm.states)-1, hmm.LENG))
				#add the HMM to the list
				ret.append(hmm)
				state = self.S_default
				continue

			#parse a header line
			if state == self.S_header:
				if re.match(r'^HMM\s', line):
					#Parsed each header line
					#check all the required options were present
					for o in REQUIRED:
						try:
							if not getattr(hmm, o):
								self._addError('Option \'{}\' is required'.format(o), False)
						except AttributeError:
							self._addError('Option \'{}\' is required'.format(o), False)
							
							if o == 'LENG':
								setattr(hmm, o, 0)
							else:
								setattr(hmm, o, '')	
					
					#check if the line after next is a COMPO line
					if re.match(r'^COMPO\s', lines[self.lineNum+2]):
						hmm.COMPO = self._parse_prob(lines[self.lineNum+2].split()[1:], hmm.K)
						#drop the next two lines
						consume(iterator, 2)
					else:
						#otherwise drop only one line
						consume(iterator, 1)
					
					model_state = State()
					state = self.S_model_ie
					continue
				
				m = self.hdr_re.match(line)
				if not m:
					self._addError('Invalid header line')
					continue
				
				key = m.group('key').upper()
				val = m.group('value')
				#check if the key is valid
				if not (key in OPTIONS):
					self._addWarning('Ignoring unknown option \'%s\'' % key)
					continue

				#simple strings
				if key in ['NAME', 'ACC', 'DESC', 'ALPH', 'DATE',]:
					setattr(hmm, key, val)
					if key == 'ALPH':
						if hmm.ALPH.upper() in ALPHABETS:
							hmm.SYMBOLS = ALPHABETS[hmm.ALPH.upper()]
							hmm.K = len(hmm.SYMBOLS)
						else:
							self._addError('ALPH must be \'DNA\', \'RNA\' or \'AMINO\'')
							continue
				#integers
				elif key in ['LENG', 'NSEQ', 'CKSUM',]:
					try:
						val = int(val)
						if val < 0:
							raise ValueError
						setattr(hmm, key, val)
					except ValueError:
						self._addError('{} must be a positive integer'.format(key))
						continue
				#bools
				elif key in ['RF', 'CS', 'MAP',]:
					b = {'yes': True, 'no': False,}
					val = val.lower()
					if val in b:
						setattr(hmm, key, b[val])
					else:
						self._addError('{} must be \'yes\' or \'no\''.format(key))
						continue
				#COM
				elif key == 'COM':
					m2 = re.match(r'(\d+)\s+(\S+)$', val)
					if m2:
						try:
							hmm.COM.append( (int(m2.group(1)), m2.group(2)))
						except ValueError:
							hmm.COM.append( (len(self.COM)+1, m2.group(0)))
							continue
					else:
						hmm.COM.append( (len(self.COM)+1, value))
				elif key == 'EFFN':
					try:
						val = float(val)
						if val < 0:
							raise ValueError
						hmm.EFFN = val
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
							setattr(hmm, key, (v1, v2))
						except ValueError:
							self._addError('{} must be two positive reals'.format(key))
							continue
					else:
						self._addError('{} must be two positive reals'.format(key))
						continue
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
								hmm.STATS.append((s1, s2, f1, f2))
							else:
								if s1 != 'LOCAL':
									self._addError('s1 must equal \'LOCAL\'')
								else:
									self._addError('s2 must be \'MSV\', \'VITERBI\' or \'FORWARD\'')
								continue
						except ValueError:
							self._addError('STATS <s1> <s2> <f1> <f2>')
							continue
					else:
						self._addError('STATS <s1> <s2> <f1> <f2>')
						continue
			
			if state == self.S_model_me:
				#parse the ME line
				l = line.split()
				
				model_state = State()

				try:
					model_state.num = int(l[0])
					if model_state.num != len(hmm.states):
						self._addError('Expected state number %s' % len(hmm.states))
				except ValueError:
					self._addError('Node number must be a positive integer')

				model_state.me = self._parse_prob(l[1:hmm.K+1], hmm.K)
				
				#MAP number
				try:
					if(l[hmm.K+1] != '-'):
						model_state.MAP = int(l[hmm.K+1])
					elif hmm.MAP:
						self._addWarning('Map annotation is \'-\', even though MAP is \'yes\'')
				except ValueError:
					self._addError('Map Annotation must be an integer or \'-\'')
				except IndexError:
					self._addError('No Map annotation provided')
				#RF annotation
				try:
					model_state.RF = l[hmm.K+2]
					if len(model_state.RF) != 1:
						self._addError('RF annotation must be a single character')
				except IndexError:
					self._addError('No RF annotation provided')
				#CS annotation
				try:
					model_state.CS = l[hmm.K+3]
					if len(model_state.CS) != 1:
						self._addError('CS annotation must be a single character')
				except IndexError:
					self._addError('No CS annotation provided')
				#we're now expecting an IE line
				state = self.S_model_ie
				continue

			if state == self.S_model_ie:
				model_state.ie = self._parse_prob(line.split(), hmm.K)
				state = self.S_model_st
				continue

			if state == self.S_model_st:
				model_state.tr = self._parse_prob(line.split(), 7)
				#add the state to the current hmm
				hmm.states.append(model_state)
				model_state = State()
				state = self.S_model_me
				continue

		#parsed every line in the file
		#if there were errors, raise them. Otherwise keep quiet
		if len(self.errors):
			self._raiseErrors()

		return ret
	
	def _raiseErrors(self):
		raise HMMFileException(self.errors, self.warnings)

	def _addError(self, what, show_line=True):
		if(show_line):
			self.errors.append((self.lineNum, what))
		else:
			self.errors.append((-1, what))

	def _addWarning(self, what, show_line=True):
		if(show_line):
			self.warnings.append((self.lineNum, what))
		else:
			self.warnings.append((-1, what))


	def _parse_prob(self, l, expected=-1):
		"""Parse a list of log probabilities"""
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

def read(f):
	"""Read in all the HMMs in f"""
	p = HMMParser()
	return p.read(f)

def write(hmms, f):
	"""Write out the HMMs in hmms to f, either a fileObject or a string"""	
	if isinstance(f, basestring):
		f = open(f, 'w')

	if not isinstance(hmms, collections.Iterable):
		hmms = [hmms,]

	for hmm in hmms:
		f.write('HMMER3/b [pyHMMER | 2012]\n')
		#write header
		for o in ['NAME', 'ACC', 'DESC', 'LENG', 'ALPH', 'DATE', 'NSEQ', 'EFFN',
				'CKSUM',]:
			try:
				f.write('{:<5s} {}\n'.format(o, getattr(hmm, o)))
			except AttributeError:
				pass
		for o in ['RF', 'CS', 'MAP']:
			try:
				if getattr(hmm, o):
					f.write('{:<5s} yes\n'.format(o))
				else:
					f.write('{:<5s} no\n'.format(o))
			except AttributeError:
				pass

		for c in hmm.COM:
			f.write('COM   %s %s\n' % c)

		for o in ['GA', 'TC', 'NC',]:
			try:
				out = getattr(hmm, o)
				f.write('{:<5s} {:.2f}  {:.2f}\n'.format(o, out[0], out[1]))
			except AttributeError:
				pass

		for s in hmm.STATS:
			f.write('STATS {:<5s} {:<10s} {:8f} {:8f}\n'.format(*s))

		#write HMM lines
		f.write('HMM ')
		f.write( '%9s'*len(hmm.SYMBOLS) % tuple(hmm.SYMBOLS))
		f.write(
			'\n            m->m     m->i     m->d'
			'     i->m     i->i     d->m     d->d\n')

		def ff(f, l=9, p=5):
			try:
				return ('%'+str(l)+'.'+str(p)+'f') % f
			except TypeError:
				return ' '*(l-1) + '*'
		#write CMPO, if it exists. Should probably calculate it...
		try:
			f.write(("%7s" + ("%9s" * hmm.K) + "\n") % 
					(('COMPO',) + tuple((ff(f) for f in hmm.COMPO))))
		except AttributeError:
			pass

		#write model
		me_fmt = "%7i" + ("%9s" * hmm.K) + "%7i %1s %1s\n"
		ie_fmt = (" " * 7) + ("%9s" * hmm.K) + "\n"
		tr_fmt = (" " * 7) + ("%9s" * 7) + "\n"

		#0th state gets special treatment
		if len(hmm.states) < 1:
			f.write('      //\n')
			continue

		f.write(ie_fmt % tuple((ff(f) for f in hmm.states[0].ie)))
		f.write(tr_fmt % tuple((ff(f) for f in hmm.states[0].tr)))

		#write model
		for i,s in enumerate(hmm.states[1:], 1):
			f.write(me_fmt % ((i,) + tuple((ff(f) for f in s.me)) + 
				(s.MAP, s.RF, s.CS)))
			f.write(ie_fmt % tuple((ff(f) for f in s.ie)))
			f.write(tr_fmt % tuple((ff(f) for f in s.tr)))

		#end with //
		f.write('      //\n')
		




class State:
	"""A state within an HMM"""
	def __init__(self):
		self.num = -1
		self.MAP = -1
		self.RF = '-'
		self.CS = '-'
		self.me = []
		self.ie = []
		self.tr = []

class HMM:
	"""A Hidden Markov Model"""
	def __init__(self):
		self.NAME = ''
		self.LENG = 0

		self.states = []
		self.COM = []
		self.STATS = []
		#assume these to be false if not present
		self.RM = self.CS = self.MAP = False
		self.K = 0
		

