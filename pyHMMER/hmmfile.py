"""Read and Write HMMER's .hmm format"""

import re, itertools, collections, math

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
				if len(hmm.states)-1 < len(hmm):
					self._addError('Too few states in model ({:d} < {:d})'
							.format(len(hmm.states)-1, len(hmm)))
				if len(hmm.states)-1 > len(hmm):
					self._addError('Too many states in model ({:d} > {:d})'
							.format(len(hmm.states)-1, len(hmm)))
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
							if not getattr(hmm, o.lower()):
								self._addError('Option \'{}\' is required'.format(o), False)
						except AttributeError:
							self._addError('Option \'{}\' is required'.format(o), False)
							
							if o == 'LENG':
								setattr(hmm, o.lower(), 0)
							else:
								setattr(hmm, o.lower(), '')	
					
					#check if the line after next is a COMPO line
					if re.match(r'^COMPO\s', lines[self.lineNum+2]):
						hmm.compo = self._parse_prob(lines[self.lineNum+2].split()[1:], hmm.K)
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
					setattr(hmm, key.lower(), val)
					if key == 'ALPH':
						if hmm.alph.upper() in ALPHABETS:
							hmm.symbols = ALPHABETS[hmm.alph.upper()]
							hmm.K = len(hmm.symbols)
							hmm.alpha = hmm.alph.upper()
						else:
							self._addError('ALPH must be \'DNA\', \'RNA\' or \'AMINO\'')
							continue
				#integers
				elif key in ['LENG', 'NSEQ', 'CKSUM',]:
					try:
						val = int(val)
						if val < 0:
							raise ValueError
						setattr(hmm, key.lower(), val)
					except ValueError:
						self._addError('{} must be a positive integer'.format(key))
						continue
				#bools
				elif key in ['RF', 'CS', 'MAP',]:
					b = {'yes': True, 'no': False,}
					val = val.lower()
					if val in b:
						setattr(hmm, key.lower(), b[val])
					else:
						self._addError('{} must be \'yes\' or \'no\''.format(key))
						continue
				#COM
				elif key == 'COM':
					m2 = re.match(r'(\d+)\s+(\S+)$', val)
					if m2:
						try:
							hmm.com.append( (int(m2.group(1)), m2.group(2)))
						except ValueError:
							hmm.com.append( (len(hmm.com)+1, m2.group(0)))
							continue
					else:
						hmm.com.append( (len(hmm.com)+1, val))
				elif key == 'EFFN':
					try:
						val = float(val)
						if val < 0:
							raise ValueError
						hmm.effn = val
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
							setattr(hmm, key.lower(), (v1, v2))
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
								hmm.stats.append((s1, s2, f1, f2))
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
						model_state.map = int(l[hmm.K+1])
					elif hmm.map:
						self._addWarning('Map annotation is \'-\', even though MAP is \'yes\'')
				except ValueError:
					self._addError('Map Annotation must be an integer or \'-\'')
				except IndexError:
					self._addError('No Map annotation provided')
				#RF annotation
				try:
					model_state.rf = l[hmm.K+2]
					if len(model_state.rf) != 1:
						self._addError('RF annotation must be a single character')
				except IndexError:
					self._addError('No RF annotation provided')
				#CS annotation
				try:
					model_state.cs = l[hmm.K+3]
					if len(model_state.cs) != 1:
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

		f.write('{:<5s} {}\n'.format('NAME', hmm.name if hmm.name else 
			"<Untitled HMM>"))
		f.write('{:<5s} {}\n'.format('LENG', len(hmm.states) - 1))
		
		for o in ['ACC', 'DESC', 'ALPH', 'DATE', 'NSEQ', 'EFFN',
				'CKSUM',]:
			try:
				f.write('{:<5s} {}\n'.format(o, getattr(hmm, o.lower())))
			except AttributeError:
				pass
		for o in ['RF', 'CS', 'MAP']:
			try:
				if getattr(hmm, o.lower()):
					f.write('{:<5s} yes\n'.format(o))
				else:
					f.write('{:<5s} no\n'.format(o))
			except AttributeError:
				pass

		for c in hmm.com:
			f.write('COM   %s %s\n' % c)

		for o in ['GA', 'TC', 'NC',]:
			try:
				out = getattr(hmm, o.lower())
				f.write('{:<5s} {:.2f}  {:.2f}\n'.format(o, out[0], out[1]))
			except AttributeError:
				pass

		for s in hmm.stats:
			f.write('STATS {:<5s} {:<10s} {:8f} {:8f}\n'.format(*s))

		#write HMM lines
		f.write('HMM ')
		f.write( '%9s'*len(hmm.symbols) % tuple(hmm.symbols))
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
					(('COMPO',) + tuple((ff(f) for f in hmm.compo))))
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
				(s.map, s.rf, s.cs)))
			f.write(ie_fmt % tuple((ff(f) for f in s.ie)))
			f.write(tr_fmt % tuple((ff(f) for f in s.tr)))

		#end with //
		f.write('      //\n')
		




class State:
	"""A state within an HMM"""
	def __init__(self):
		self.num = -1
		self.map = -1
		self.rf = '-'
		self.cs = '-'
		self.me = []
		self.ie = []
		self.tr = []

	def __eq__(self, rhs):
		return (self.num == rhs.num and
						self.map == rhs.map and
						self.rf  == rhs.rf  and
						self.cs  == rhs.cs  and
						self.me  == rhs.me  and
						self.ie  == rhs.ie  and
						self.tr  == rhs.tr)
	
	def __ne__(self, rhs):
		return not self.__eq__(rhs)

	def __unicode__(self):
		return ("HMM State num: {}" + ("\n  {}: {}"*6)).format(
				self.num, 'map', self.map, 'rf', self.rf, 'cs', self.cs, 
				'me', self.me, 'ie', self.ie, 'tr', self.tr)

	def __str__(self):
		return unicode(self).encode('utf-8')

	def __repr__(self):
		return "HMM State({})".format(self.num)

class HMM:
	"""A Hidden Markov Model"""
	_trorder = ['MM','MI','MD','IM','II','DM','DD']

	def __init__(self, name='', alphabet=None):
		self.name = name
		self.leng = 0

		self.states = []
		self.com = []
		self.stats = []
		#assume these to be false if not present
		self.rm = self.cs = self.map = False
		self.alph = self.alpha = alphabet
		if self.alph:
			self.symbols = ALPHABETS[self.alpha.upper()]
			self.K = len(self.symbols)
		else:
			self.symbols = ''
			self.K = 0

	def clean(self):
		"""Makes sure the states adhere to the requirements placed on the first and
		last states - call this after the model is built"""
		#State0: (DM, DD) = (1, 0)
		self.states[0].tr[5:7] = [0.0, '*',]
		#State0 is mute
		self.states[0].me = []

		#FinalState: (DM, DD) = (1, 0)
		s = self.states[-1]
		s.tr[5:7] = [0.0, '*',]
		# MD = 0
		s.tr[2] = '*'
		s.tr[0:2] = logodds(self._norm(expodds(s.tr[0:2])))
	
	def addState(self, transition = None, emission = None, insert_emission = None,
								MAP=-1, CS='-', RF='-'):
		"""Add a new state with the properties specified
			
			transtision: optional dict containing the optional keys 
					MM, MI, MD, IM, II, DM, DD
				or lower case equivalents

			emission: dictionary containing the match emission probabilities for each
				symbol in the alphabet. Omitted symbols are given p=0.
			
			insert_emission: dictionary containing the insert emission probabilities 
				for each symbol in the alphabet. Omitted symbols are given p=0.

			MAP: MAP annotation for this node, integer

			RF: Reference Annotation for this node, char

			CS: Consensus Structure annotation for this node
		"""
		#print ("{{\n\"transition\"={},\n\"emission\"={},\n\"insert_emission\"={},\n"
		#	+ "\"MAP\"={},\n\"CS\"=\"{}\",\n\"RF\"=\"{}\"\n}},").format(transition, emission,
		#		insert_emission, MAP, CS, RF)
		
		#defaults
		if not transition:
			transition = {'MM': 1.0,}
		if not emission:
			emission = {}
			for letter in self.symbols:
				emission[letter] = 1
		if not insert_emission:
			insert_emission = {}
			for letter in self.symbols:
				insert_emission[letter] = 1

		#transitions
		tr = self._gettr(transition)

		#emissions
		em = self._getem(emission, 'match emissions')
		iem= self._getem(insert_emission, 'insert emissions')

		s = State()
		s.tr=tr
		s.me=em
		s.ie=iem
		s.num=len(self.states)
		s.map=MAP
		s.cs=CS
		s.rf=RF

		self.states.append(s)

	def _gettr(self, transition):
		#get the values
		tr = [-1,]*7
		for (k,v) in transition.iteritems():
			if v < 0:
				raise ValueError('Value of \'{}\' is less than zero in transision'
						.format(k))
			try:
				i = self._trorder.index(k.upper())
			except ValueError:
				raise ValueError('Unknown key \'{}\' in transition'.format(k))
			if tr[i] >= 0:
				raise ValueError('Duplicate key \'{}\' in transition'.format(k))
			tr[i] = v

		#make sure no special values persist
		for i in range(len(tr)):
			if tr[i] < 0:
				tr[i] = 0

		#set defaults and normalise
		if sum(tr[0:3]) == 0:
			tr[0:3] = [1,0,0,]
		tr[0:3] = self._norm(tr[0:3])

		if sum(tr[3:5]) == 0:
			tr[3:5] = [1,0,]
		tr[3:5] = self._norm(tr[3:5])

		if sum(tr[5:7]) == 0:
			tr[5:7] = [1,0,]
		tr[5:7] = self._norm(tr[5:7])

		return logodds(tr)				

	def _getem(self, emission, name):
		r = [0,] * self.K

		for (k,v) in emission.iteritems():
			if v < 0:
				raise ValueError("Value \'{}\' less than zero in {}".format(k,name))
			if emission.has_key(k.lower()) and emission.has_key(k.upper()):
				raise ValueError("Duplicate key \'{}\' in {}".format(k,name))
			i = ALPHABETS[self.alpha].find(k.upper())
			if i < 0:
				raise ValueError("Unknown symbol \'{}\' not in alphabet ({}[\'{}\']) when"
						" parsing {}".format(k,self.alpha, ALPHABETS[self.alpha], name))
			r[i] = v

		if sum(r) == 0:
			r = [1,]*self.K

		return logodds(self._norm(r))

	def _norm(self, l):
		t = sum(l)
		if t == 0:
			t = 1
		r = []
		for v in l:
			r.append(float(v) / t)
		return r

	def __len__(self):
		return self.leng

	def __nonzero__(self):
		return True

	def __unicode__(self):
		s = "HMM with {} {}.\n\t".format(len(self.states), 
				"states" if len(self.states)>1 else "state")
		for state in self.states:
			if state.me:
				s += ALPHABETS[self.alpha][state.me.index(max(state.me))]
		return s

	def __str__(self):
		return self.__unicode__().encode('utf-8')

def logodds(p):
	"""Takes a single or an interable of probabilities and returns -log(p)"""
	if hasattr(p, '__iter__'):
		r = []
		for i in p:
			if i <= 0:
				r.append('*')
			else:
				r.append(-math.log(i))
		return r

	if p <= 0:
		return '*'
	return -math.log(p)

def expodds(p):
	if hasattr(p, '__iter__'):
		r = []
		for i in p:
			if isinstance(i, basestring):
				r.append(0.0)
			else:
				r.append(math.exp(-i))
		return r

	if isinstance(p, basestring):
		return 0.0
	return math.exp(-p)

