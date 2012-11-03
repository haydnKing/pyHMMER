"""Read and Write HMMER's match format"""
import re, sequtils

class Match:
	"""Represents an HMM match"""
	#Format strings for output
	fmt = ( 
		"{:20s} {:20s} {:4d} {:7g} {:7g} {:7g} {:9d} {:9d} {:9d} {:9d} {:9d} {:9d}"
			+ " {:s}"
		)
	hdr = (
		("{:20s} {:20s} {:4s} {:7s} {:7s} {:7s} {:9s} {:9s} {:9s} {:9s} {:9s} {:9s}"
			+ " {:s}")
		.format("Target", "Query", "#", "c-Evalue", "i-Evalue", "Score",
			"hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
			"description")
		)
	
	def __init__(self):
		self.target = None
		self.query = None
		self.number = 0
		self.c_evalue = 0
		self.i_evalue = 0
		self.score = 0

		self.hmm_from = 0
		self.hmm_to = 0

		self.ali_from = 0
		self.ali_to = 0

		self.env_from = 0
		self.env_to = 0

		self.description = ""

		self.translation = {'query': 'DNA', 'target': 'DNA',}
		self.frame = 0

	def getTarget(self):
		if self.target:
			return self.target
		return "<Unknown>"

	def getFrame(self):
		if self.frame and self.frame in range(-3,3):
			return self.frame
		return 0

	def getQuery(self):
		if self.query:
			return self.query
		return "<Unknown>"

	def getHMMPos(self):
		"""Get the match's position withing the HMM"""
		return (self.hmm_from, self.hmm_to)

	def getAliPos(self):
		return self._map_position((self.ali_from, self.ali_to))

	def getEnvPos(self):
		return self._map_position((self.env_from, self.env_to))

	def getFrameSpan(self, mode='hmm'):
		"""get the raw coordinates of the match"""
		mode = mode.lower()
		if mode == 'hmm':
			start = self.ali_from - self.hmm_from
			end = self.ali_to + (self.query.LENG - self.hmm_to)
		elif mode == 'env':
			start = self.env_from
			end = self.env_to
		elif mode == 'ali':
			start = self.ali_from
			end = self.ali_to
		else:
			raise ValueError("Invalid mode \'{}\', " +
				"must be \'hmm\', \'env\' or \'ali\'".format(mode))
		return (start,end)

	def getTargetSpan(self, mode='hmm'):
		"""Get the span of the match mapped onto the target sequence"""	
		return self._map_position(self.getFrameSpan())

	def isTranslation(self):
		"""return true if the target and hmm have different alphabets"""
		return self.translation['query'].lower() != self.translation['target'].lower()

	def getSequence(self, mode='hmm'):
		"""get the sequence of the match
			mode: either 'hmm', 'ali' or 'env' depending on which region of the 
				sequence is required
		"""
		span = self.getTargetSpan(mode)
		
		if span[0] > span[1]:
			return self.target.seq[span[1]:span[0]].reverse_complement()
		return self.target.seq[span[0]:span[1]]

	def _map_position(self, pos):
		"""Map the position given onto the target"""
		if self.isTranslation():
			if (self.translation['target'].upper() in ['DNA', 'RNA'] and 
					self.translation['query'].upper() == 'AMINO'):
				#get the length of the target sequence
				l = len(self.target.seq)
				#and find my position within it
				ret = (3 * pos[0], 3* pos[1])
				if self.frame not in [1,2,3,-1,-2,-3]:
					raise ValueError("Nonsensical Frame \'{}\'".format(self.frame))
				if self.frame > 0:
					ret = tuple(x+self.frame-1 for x in ret)
				elif self.frame < 0:
					ret = tuple(l-x for x in ret)
				return ret
			else:
				raise ValueError("Unhandled Translation {query} to {target}"
						.format(**self.translation))
		else:
			return pos
	



	def __unicode__(self):
		return self.fmt.format(self.target.name, self.query.NAME, self.number,
				self.c_evalue, self.i_evalue, self.score, self.hmm_from, self.hmm_to,
				self.ali_from, self.ali_to, self.env_from, self.env_to, 
				self.desc)

	def __str__(self):
		return unicode(self).encode('utf-8')

format_list = [
		{'name': 't_name',			'len': 19, 'just': '<', 'type': 's', 
				'lname':'target name',},
		{'name': 't_accession', 'len': 10, 'just': '<', 'type': 's',
				'lname':'accession',},
		{'name': 't_len',				'len':  5, 'just': '>', 'type': 'd',
				'lname':'tlen',},
		{'name': 'q_name',			'len': 20, 'just': '<', 'type': 's',
				'lname':'query name',},
		{'name': 'q_accession',	'len': 10, 'just': '<', 'type': 's',
				'lname':'accession',},
		{'name': 'q_len',				'len':  5, 'just': '>', 'type': 'd',
				'lname':'qlen',},
		{'name': 'f_evalue',		'len':  9, 'just': '>', 'type': 'g', 'prec':2,
				'lname':'E-value',},
		{'name': 'f_score',			'len':  6, 'just': '>', 'type': 'f', 'prec':1,
				'lname':'score',},
		{'name': 'f_bias',			'len':  5, 'just': '>', 'type': 'f', 'prec':1,
				'lname':'bias',},
		{'name': 'num',					'len':  3, 'just': '>', 'type': 'd',
				'lname':'#',},
		{'name': 'of',					'len':  3, 'just': '>', 'type': 'd',
				'lname':'of',},
		{'name': 'c_evalue',		'len':  9, 'just': '>', 'type': 'g', 'prec':2,
				'lname':'c-Evalue',},
		{'name': 'i_evalue',		'len':  9, 'just': '>', 'type': 'g', 'prec':2,
				'lname':'i-Evalue',},
		{'name': 'score',				'len':  6, 'just': '>', 'type': 'f', 'prec':1,
				'lname':'score',},
		{'name': 'bias',				'len':  5, 'just': '>', 'type': 'f', 'prec':1,
				'lname':'bias',},
		{'name': 'hmm_from',		'len':  5, 'just': '>', 'type': 'd',
				'lname':'from',},
		{'name': 'hmm_to',			'len':  5, 'just': '>', 'type': 'd',
				'lname':'to',},
		{'name': 'ali_from',		'len':  5, 'just': '>', 'type': 'd',
				'lname':'from',},
		{'name': 'ali_to',			'len':  5, 'just': '>', 'type': 'd',
				'lname':'to',},
		{'name': 'env_from',		'len':  5, 'just': '>', 'type': 'd',
				'lname':'from',},
		{'name': 'env_to',			'len':  5, 'just': '>', 'type': 'd',
				'lname':'to',},
		{'name': 'acc',					'len':  4, 'just': '>', 'type': 'f', 'prec':2,
				'lname':'acc',},
		{'name': 'desc',				'len': 21, 'just': '<', 'type': 's',
				'lname':'description of target',},
	]


def save(matches, f):
	"""Save the match or matches in matches to the file f, either a string
		filename of a file object"""
	#should I close the file object?
	should_close = False
	if isinstance(f, basestring):
		f = open(f, 'w')
 		should_close = True

	#for each output match
	for o in matches:
		#and for each item
		for fmt in format_list:
			if fmt.has_key('prec'):
				sfmt = "{{.{name}:{just}.{prec}{type}}}".format(**fmt)
			else:
				sfmt = "{{.{name}:{just}{type}}}".format(**fmt)
			l = len(sfmt.format(o))
			#increase the length if the formatted length is longer
			fmt['len'] = max(l, fmt['len'])

	#Actual output

	#write header
	gap = sum((fmt['len']+1 for fmt in format_list[0:6]))
	hfmt = "#" + (" " * (gap-1))
	gap = sum((fmt['len']+1 for fmt in format_list[6:9])) -1
	hfmt += "{{:-^{}s}} ".format(gap)
	gap = sum((fmt['len']+1 for fmt in format_list[9:15])) -1
	hfmt += "{{:-^{}s}} ".format(gap)

	gap = sum((fmt['len']+1 for fmt in format_list[15:17])) -1
	hfmt += "{{:^{}s}} ".format(gap)
	gap = sum((fmt['len']+1 for fmt in format_list[17:19])) -1
	hfmt += "{{:^{}s}} ".format(gap)
	gap = sum((fmt['len']+1 for fmt in format_list[19:21])) -1
	hfmt += "{{:^{}s}}\n".format(gap)
	
	f.write(hfmt.format(" full sequence ", " this domain ", "hmm coord",
		"ali coord", "env coord"))

	dash = "#" + ("-" * (format_list[0]['len']-1))
	hfmt = "# {{:{}s}} ".format(format_list[0]['len'] - 2)
	for fmt in format_list[1:]:
		hfmt += "{{:{}{}s}} ".format(fmt['just'], fmt['len'])
		dash += ' ' + ('-' * fmt['len'])
	hfmt = hfmt.rstrip() + "\n"
	dash += "\n"

	f.write(hfmt.format( *tuple( fmt['lname'] for fmt in format_list) ))
	f.write(dash)

	#build the format for the main output
	sfmt = ""
	for fmt in format_list:
		if fmt.has_key('prec'):
			sfmt += "{{{name}:{just}{len:d}.{prec:d}{type}}}".format(**fmt)
		else:
			sfmt += "{{{name}:{just}{len:d}{type}}}".format(**fmt)
		if fmt != format_list[-1]:
			sfmt += " "
		else:
			sfmt += "\n"

	#convert to 1-offset and write
	out = []
	for o in matches:
		o_ = dict()
		for fmt in format_list:
			v = getattr(o, fmt['name'])
			if re.search(r"_from$|_to$", fmt['name']):
				v = int(v) + 1
			o_[fmt['name']] = v
		out.append(o_)


	for o in out:
		f.write(sfmt.format(**o))
				
	if should_close:
		f.close()


def load(f, hmms, targets):
	"""Load the results of the search stored in f (filename or file object)
		which was conducted by searching hmms against targets"""
	#make sure the arguments are iterable
	if not hasattr(hmms, "__iter__"):
		hmms = [hmms,]
	if not hasattr(targets, "__iter__"):
		targets = [targets,]

	#should I close the file object?
	should_close = False
	if isinstance(f, basestring):
		f = open(f, 'r')
		should_close = True

	#get the HMM alphabet
	if hmms:
		hmm_alpha = hmms[0].ALPH
		for h in hmms:
			if h.ALPH != hmm_alpha:
				raise ValueError("Not all HMMS have the same alphabet")

	#do the loading
	matches = []
	for line in f:
		#skip comments
		if line.lstrip()[0] == '#':
			continue

		#skip lines which aren't long enough
		l = line.split()
		if len(l) < len(format_list):
			continue

		match = Match()

		#attempt to set target
		name = l[0]
		for t in targets:
			if t.name == name:
				match.target = t
				break
		#query
		for q in hmms:
			if q.NAME == l[3]:
				match.query = q
				break

		for i,fmt in enumerate(format_list[0:-1]):
			v = l[i]
			if fmt['type'] in ['f','g',]:
				setattr(match, fmt['name'], float(v))
			elif fmt['type'] == 'd':
				if re.search(r"_from$|_to$", fmt['name']):
					v = int(v) - 1
				setattr(match, fmt['name'], int(v))
			else:
				setattr(match, fmt['name'], str(v))

		match.desc = " ".join(l[22:])

		#is there frame information attached?
		m = re.search(r"frame:\s(?P<frame>[+-]?\d)", match.desc)
		if m:
			match.frame = int(m.group("frame"))
			match.translation['target'] = 'DNA'
		if match.target:
			match.translation['target'] = sequtils.seq_type(str(match.target.seq))
		if hmm_alpha:
			match.translation['query'] = hmm_alpha

		matches.append(match)

	if should_close:
		f.close()

	return matches




