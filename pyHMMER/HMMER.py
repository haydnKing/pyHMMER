"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet, Alphabet

import cStringIO as StringIO
import hmmfile, tools, tempfile, subprocess, os, re
from subprocess import Popen, PIPE, STDOUT

class Match:
	"""Represents an HMM match"""
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
		self.cEvalue = 0
		self.iEvalue = 0
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

	def getSpan(self):
		start = self.env_from - self.hmm_from
		end = self.env_to + (self.query.LENG - self.hmm_to)
		return self._map_position((start,end))

	def isTranslation(self):
		"""return true if the target and hmm have different alphabets"""
		return self.translation['query'].lower() != self.translation['target'].lower()

	def _map_position(self, pos):
		"""Map the position given onto the target"""
		if self.isTranslation():
			if (self.translation['target'].upper() == 'DNA' and 
					self.translation['query'].upper() == 'AMINO'):
				#get the length of the target sequence
				l = len(self.target.seq)
				#and find my position within it
				ret = (3 * pos[0], 3* pos[1])
				if self.frame not in [1,2,3,-1,-2,-3]:
					raise ValueError("Nonsensical Frame \'{}\'" % self.frame)
				if self.frame > 0:
					ret = tuple(x+self.frame-1 for x in ret)
				elif self.frame < 0:
					ret = tuple(l-1-x for x in ret)
				return ret
			else:
				raise ValueError("Unhandled Translation {query} to {target}"
						.format(**self.translation))
		else:
			return pos

	def __unicode__(self):
		return self.fmt.format(self.target.name, self.query.NAME, self.number,
				self.cEvalue, self.iEvalue, self.score, self.hmm_from, self.hmm_to,
				self.ali_from, self.ali_to, self.env_from, self.env_to, 
				self.description)

	def __str__(self):
		return unicode(self).encode('utf-8')
	
		

class hmmsearch:
	"""Search for an HMM in a database and collect the results"""
	def __init__(self, hmm, targets):
		"""Perform the search
				hmm: a file name or an HMM object which has been loaded from a file
				targets: the sequences to search - a fasta filename or one or 
					more Bio.SeqRecord

				If the hmm performs searches on Amino Acids and and of the inputs are
				DNA sequences, 6-frame translations will be produced automatically
				Reverse translations (from Amino Acid to DNA) are not supported
		"""
		# Load the HMM(s)
		if not hasattr(hmm, "__iter__"):
			hmm = [hmm,]
		#load the file if h is not an HMM object
		self.hmm = list()
		for h in hmm:
			if not isinstance(h, hmmfile.HMM):
				self.hmm = self.hmm + hmmfile.read(h)
			else:
				self.hmm.append(h)

		#make sure targets is iterable
		if not hasattr(targets, '__iter__'):
			targets = [targets,]
		self.targets = list()
		#load a fasta file if t is not a seqRecord
		for t in targets:
			if not isinstance(t, SeqRecord):
				for req in SeqIO.parse(t, 'fasta'):
					self.targets.append(req)
			else:
				self.targets.append(t)

		ttargets = []
		hmm_alpha = self.hmm[0].ALPH.upper()
		for h in self.hmm:
			if h.ALPH.upper() != hmm_alpha:
				raise ValueError("The HMMs don't all have the same alphabet")
		#Translate targets if necessary
		for t in self.targets:
			if hmm_alpha == "DNA":
				if (isinstance(t.seq.alphabet, DNAAlphabet) or 
						isinstance(t.seq.alphabet, Alphabet)):
					ttargets.append(t)
				else:
					raise ValueError("No translation available for %s to %s" % 
							(hmm_alpha, t.seq.alphabet))
			elif hmm_alpha == "RNA":
				raise ValueError("RNA HMMs are not supported")
			elif hmm_alpha == "AMINO":
				if isinstance(t.seq.alphabet, ProteinAlphabet):
					ttargets.append(t)
				else:
					ttargets = ttargets + tools.getSixFrameTranslation(t)

		#write the HMM to a temporary file
		hmm_file = tempfile.NamedTemporaryFile()
		target_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()

		hmmfile.write(self.hmm, hmm_file)
		hmm_file.flush()

		SeqIO.write(ttargets, target_file, 'fasta')
		target_file.flush()
		del ttargets

		p = Popen(['hmmsearch', '--tformat', 'fasta', 
			'--domtblout', out_file.name, hmm_file.name, target_file.name,], 
				stdout=PIPE, stdin=PIPE, stderr=PIPE)
		out = p.communicate()


		#TODO: Check stdout for errors

		self.matches = []
		for line in out_file:
			#skip comments
			if line.lstrip()[0] == '#':
				continue

			#skip lines which aren't long enough
			l = line.split()
			if len(l) < 21:
				continue

			match = Match()
			#target
			name = l[0]
			

			for t in self.targets:
				if t.name == name:
					match.target = t
					break
			#query
			for q in self.hmm:
				if q.NAME == l[3]:
					match.query = q
					break

			match.number = int(l[9])
			match.cEvalue = float(l[11])
			match.iEvalue = float(l[12])
			match.score = float(l[13])

			match.hmm_from = int(l[15])
			match.hmm_to = int(l[16])
			match.ali_from = int(l[17])
			match.ali_to = int(l[18])
			match.env_from = int(l[19])
			match.env_to = int(l[20])

			match.description = " ".join(l[22:])

			m = re.search(r"frame:\s(?P<frame>[+-]?\d)", match.description)
			if m:
				print "m.group(0) = %s" % m.group(0)
				match.frame = int(m.group("frame"))
				match.translation['target'] = 'DNA'
			match.translation['query'] = hmm_alpha

			self.matches.append(match)

		out_file.close()
		hmm_file.close()
		target_file.close()

	def discardOverlaps(self):
		"""Remove overlaping sequences from the results
		"""
		#construct a list of ranges for each target frame
		l = dict()
		for t in self.targets:
			l[t.name] = list()

		class minimatch:
			def __init__(self, m):
				self.match = m
				self.span = m.getSpan()
				self.span = (min(self.span), max(self.span),)
				self.score = m.score

			def overlaps(self, m):
				# either m starts in me OR m ends in me OR m contains me
				return ( (self.span[0] < m.span[0] and self.span[1] > m.span[0]) or
								 (self.span[0] < m.span[1] and self.span[1] > m.span[1]) or
								 (self.span[0] > m.span[0] and self.span[1] < m.span[1]))

		#add each match to the list
		for m in self.matches:
			l[m.target.name].append(minimatch(m))

		for m in l.itervalues():
			#sort ascending start positions
			m.sort(key=lambda item: item.span[0])
			while len(m):
				#get all the matches in the first conflicting group
				conflict = [m[0],]
				for m_ in m[1:]:
					#if m_ starts after the conflict ends
					if m_.span[0] > conflict[-1].span[1]:
						break
					#otherwise add it to the conflict
					conflict.append(m_)
				
				#deal with all the items in this conflicting group
				for c in conflict:
					m.remove(c)

				#special case where there is only one item
				if len(conflict) == 1:
					#accept the item
					continue
				
				#otherwise, sort by score
				conflict.sort(key=lambda item: item.score)
				conflict.reverse()
				#continue until there aren't any more conflicts
				while len(conflict):
					#remove any items which overlap with the highest scoring item
					for m_ in conflict[1:]:
						if conflict[0].overlaps(m_):
							#m_ overlaps a higher scoring match
							self.matches.remove(m_.match)
							conflict.remove(m_)
					#conflict[0] no longer conflicts with anything, so accept it
					conflict.pop(0)


	def extractSequences(self):
		"""Extract the sequences from the remaining matches from the database"""
		pass

	def __unicode__(self):
		ret = "Found {:d} matches\n".format(len(self.matches))
		if self.matches:
			ret += self.matches[0].hdr + '\n'
			for m in self.matches:
				ret += str(m) + '\n'
		return ret

	def __str__(self):
		return unicode(self).encode('utf-8')

	def __repr__(self):
		return self.__str__()

	def __getitem__(self, i):
		return self.matches[i]


def test_deps():
	"""Test that all the required binaries are present"""
	bins = [('hmmalign',   'HMMER', 'hmmer.janelia.org'), 
					('hmmbuild',   'HMMER', 'hmmer.janelia.org'), 
					('hmmconvert', 'HMMER', 'hmmer.janelia.org'), 
					('hmmemit',    'HMMER', 'hmmer.janelia.org'), 
					('hmmfetch',   'HMMER', 'hmmer.janelia.org'), 
					('hmmpress',   'HMMER', 'hmmer.janelia.org'), 
					('hmmscan',		 'HMMER', 'hmmer.janelia.org'), 
					('hmmsearch',  'HMMER', 'hmmer.janelia.org'), 
					('hmmsim',		 'HMMER', 'hmmer.janelia.org'), 
					('hmmstat',		 'HMMER', 'hmmer.janelia.org'), 
					('jackhmmer',  'HMMER', 'hmmer.janelia.org'), 
					('phmmer',	   'HMMER', 'hmmer.janelia.org'), 
					('clustalo',   'Clustal Omega', 'www.clustal.org/omega/'),
				 ]
	failed = []

	def which(program):
		import os
		def is_exe(fpath):
			return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

		fpath, fname = os.path.split(program)
		if fpath:
			if is_exe(program):
				return program
		else:
			for path in os.environ["PATH"].split(os.pathsep):
				exe_file =	os.path.join(path, program)
				if is_exe(exe_file):
					return exe_file

		return None

	for b in bins:
		if not which(b[0]):
			print b
			failed.append('Failed to find \'{}\' from package \'{}\'. ' 
					'Please install it from {}'.format(*b))

	if failed:
		raise ImportError('\n'.join(failed))

test_deps()
