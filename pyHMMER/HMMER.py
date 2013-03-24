"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet, Alphabet, IUPAC

import cStringIO as StringIO
import hmmfile, tools, tempfile, subprocess, os, re, resource, psutil
from subprocess import Popen, PIPE, STDOUT

import matchfile
import sequtils	

class wrap(object):
	"""A class which wraps anothe class, overriding specific values"""
	def __init__(self, wrapped, overrides):
		self.wrapped = wrapped
		self.overrides = overrides

	def __getattr__(self, item):
		#if the item starts with an underscore, return the unwrapped version
		if item[0] == '_' and self.overrides.has_key(item[1:]):
			return getattr(self.wrapped, item[1:])
		if self.overrides.has_key(item):
			return self.overrides[item]
		return getattr(self.wrapped, item)

	def __len__(self):
		return len(self.wrapped)

def wrap_seqrecords(records, alpha = None):
	return [wrap(r, {'id': str(i), 
		'alpha': alpha or sequtils.seq_type(str(r.seq))}) 
				for i,r in enumerate(records)]

def wrap_hmms(hmms):
	return [wrap(h, {'name': str(i), 'alpha': h.alph.upper(),}) 
					for i,h in enumerate(hmms)]

def set_limits():
	"""Set process limits in the child thread"""
	rsc = resource.RLIMIT_AS
	soft, hard = resource.getrlimit(rsc)
	#limit to the total physical memory
	resource.setrlimit(rsc, (psutil.virtual_memory().total, hard))

class hmmertool:
	"""Class that each tool inherits from"""
	SWITCHES = []
	ARGS = []
	CMD = ''

	def getArgs(self, **kwargs):
		args = []
		for k,v in kwargs.iteritems():
			fmt = "--{}" if (len(k) > 1) else "-{}"
			if k in self.ARGS:
				args += [fmt.format(k), str(v)]
			elif k in self.SWITCHES:
				args += [fmt.format(k)]
			else:
				raise ValueError("Unknown {} argument \'{}\'"
						.format(self.__class__.__name__, k))
		return args

class jackhmmer(hmmertool):
	"""Iteratively seach a protein database with a protein sequence

	Attributes:
		- matches: matches found by jackhmmer (matchfile.Match)
		- hmms: list of the HMMs used during jackhmmer iterations
	"""

	SWITCHES = ['max', 'nobias', 'fast', 'hand','wpb','wgsc','wblosum','wnone',
			'eent', 'enone','mpi']
	ARGS = ['popen', 'pextend', 'mxfile', 'E', 'T', 'Z', 'domE', 'domT', 'domZ',
			'incE', 'incT', 'incdomE', 'incdomT', 'F1', 'F2', 'F3', 'symfrac',
			'fragthresh','wid','eset','ere','esigma','eid','EmL','EmN','EvL','EvN',
			'EfL','EfM','Eft','seed','cpu']

	def __init__(self, seq, seqdb, verbose=False, **kwargs):
		"""
			seq: the sequence to search with
				a single or a list of SeqRecords

			seqdb: the sequence database to search
				a single or a list of SeqRecords

			keyword arguments: other arguments to jackhmmer - see HMMER docs
		"""
		#turn everything into lists
		if isinstance(seq, SeqRecord):
			seq = [seq,]
		if isinstance(seqdb, SeqRecord):
			seqdb = [seqdb,]

		args = self.getArgs(**kwargs)

		#apply unique ids to the targets
		self.seq = wrap_seqrecords(seq, alpha='AMINO')
		self.seqdb=wrap_seqrecords(seqdb, alpha='AMINO')

		seq_file = tempfile.NamedTemporaryFile()
		seqdb_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()
		hmm_file = tempfile.NamedTemporaryFile()

		SeqIO.write(self.seq, seq_file, 'fasta')
		SeqIO.write(self.seqdb, seqdb_file, 'fasta')
		seq_file.flush()
		seqdb_file.flush()

		p = Popen(['jackhmmer', '--qformat', 'fasta', '--tformat', 'fasta', 
			'--chkhmm', hmm_file.name, '--domtblout', out_file.name,] + args + 
			[seq_file.name, seqdb_file.name,], 
				stdout=PIPE, stdin=PIPE, stderr=PIPE)
		out = p.communicate()

		self.matches = matchfile.load(out_file, self.seq, self.seqdb)

		#load the hmms
		self.hmms = []
		try:
			i = 1
			while True:
				f = "{}-{}.hmm".format(hmm_file.name, i)
				self.hmms.append(hmmfile.read(f)[0])
				os.remove(f)
				i += 1
		except IOError:
			pass

		out_file.close()
		seq_file.close()
		seqdb_file.close()

	def annotate(mode='hmm', type=None):
		"""Annotate the target seqdbs with the discovered features

			- mode: the size of the feature
							'hmm': the entire HMM
							'ali': alignment values
							'env': environment values

			- type: type to apply to the SeqFeatures, defaults to name of query seq
		"""
		for m in self.matches:
			m.getTarget().features.append(m.asSeqFeature(mode=mode,type=type))


class hmmsearch(hmmertool):
	"""Search for an HMM in a database and collect the results"""
	SWITCHES = ['cut_ga','cut_nc','cut_tc','max','nobias','nonull2',
			'stall','mpi',]
	ARGS = ['E','T','domE','domT','incE','incT','incdomE','incdomT',
			'F1','F2','F3','Z','domZ','seed','tformat','cpu',]
	
	def __init__(self, hmm = None, targets = None, **kwargs):
		"""Initialise - search if hmm and targets have been provided"""
		if hmm and targets:
			self.search(hmm, targets, **kwargs)

	def search(self, hmm, targets, **kwargs):
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
		if not hasattr(targets, '__iter__') or isinstance(targets, SeqRecord):
			targets = [targets,]
		self.targets = list()
		#load a fasta file if t is not a seqRecord
		for t in targets:
			if not isinstance(t, SeqRecord):
				for req in SeqIO.parse(t, 'genbank'):
					self.targets.append(req)
				for req in SeqIO.parse(t, 'fasta'):
					self.targets.append(req)
			else:
				self.targets.append(t)

		#apply unique ids
		self.targets = wrap_seqrecords(self.targets) 
		self.hmm = wrap_hmms(self.hmm)

		ttargets = []
		hmm_alpha = self.hmm[0].alph.upper()
		for h in self.hmm:
			if h.alph.upper() != hmm_alpha:
				raise ValueError("The HMMs don't all have the same alphabet")

		#get the arguments for HMMER
		args = self.getArgs(**kwargs)
		#clear the matches
		self.matches = []

		#Translate targets if necessary
		for t in self.targets:
			t_alpha = sequtils.seq_type(str(t.seq))
			
			#if target and hmms have same alphabet
			if hmm_alpha == t_alpha:
				self.matches += self._do_search(self.hmm, t, args)
			#else try translating
			else:
				if hmm_alpha == "AMINO" and t_alpha == "DNA":
					#looks like we have to convert
					for tt in tools.getSixFrameTranslation(t):
						self.matches += self._do_search(self.hmm, tt, args)
				else:
					raise ValueError('Unknown Translation \'{}\' to \'{}\''
							.format(t_alpha, hmm_alpha))

	def _do_search(self, hmm, target, args):
		#write the HMM to a temporary file
		hmm_file = tempfile.NamedTemporaryFile()
		target_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()

		hmmfile.write(hmm, hmm_file)
		hmm_file.flush()

		SeqIO.write([target,], target_file, 'fasta')
		target_file.flush()

		p = Popen(['hmmsearch',] + args + ['--tformat', 'fasta', 
			'--domtblout', out_file.name, hmm_file.name, target_file.name,]
			 , stdout=PIPE, stdin=PIPE, stderr=PIPE, preexec_fn=set_limits)
		out = p.communicate()

		if p.returncode != 0:
			if p.returncode == -6 and "alloc" in out[1].lower():
				seq = target.seq
				split = len(seq) / 2
				target.seq = seq[0:split]
				m_left = self._do_search(hmm, target, args)
				target.seq = seq[split:]
				m_right = self._do_search(hmm, target, args)
				for m in m_right:
					m.offset(split)
				target.seq = seq
				return m_left + m_right

			raise RuntimeError('hmmsearch error ({}) : {}'.format(
				p.returncode, out[1]))
		else:
			matches = matchfile.load(out_file, hmm, self.targets)

		out_file.close()
		hmm_file.close()
		target_file.close()
		return matches

	def mindist(self, dist=0, mode='hmm', verb=False):
		"""
			Filter the results so that each match has at least dist symbols
			between it and the next match in the same frame
		"""
		#construct a list of minimatches for each target frame
		l = self._get_by_frame(mode)

		for target in l.itervalues():
			for frame in target:
				#sort ascending start positions
				frame.sort(key=lambda item: item.span[0])
				while len(frame):
					#get all the matches in the first conflicting group
					conflict = [frame[0],]
					for m_ in frame[1:]:
						#if m_ starts after the conflict ends
						if m_.span[0] > (dist + conflict[-1].span[1]):
							break
						#otherwise add it to the conflict
						conflict.append(m_)

					#deal with all the items in this conflicting group
					for c in conflict:
						frame.remove(c)

					#special case where there is only one item
					if len(conflict) == 1:
						#accept the item
						continue
					
					#otherwise, sort by score
					conflict.sort(key=lambda item: item.match.score)
					conflict.reverse()
					#continue until there aren't any more conflicts
					while len(conflict):
						#remove any items which overlap with the highest scoring item
						for m_ in conflict[1:]:
							if conflict[0].overlaps(m_, dist):
								#m_ overlaps a higher scoring match
								self.matches.remove(m_.match)
								conflict.remove(m_)
						#conflict[0] no longer conflicts with anything, so accept it
						conflict.pop(0)

	def maxdist(self, dist=0, mode='hmm'):
		"""
			Filter the results to remove all matches which are more than dist from
			andother match in its frame
		"""
		l = self._get_by_frame(mode)
		to_remove = list()

		def d(lst):
			if len(lst) == 0:
				return
			if len(lst) == 1:
				yield (None, lst[0], None)
				return

			for i in range(0,len(lst)):
				yield (lst[i].dist(lst[i-1]), 
									lst[i], 
							lst[i].dist(lst[(i+1)%len(lst)]))

		for target in l.itervalues():
			for frame in target:
				frame.sort(key=lambda item: item.span[0])
				for a,b,c in d(frame):
					if (a > dist) and (c > dist):
						to_remove.append(b.match)

		for r in to_remove:
			try:
				self.matches.remove(r)
			except ValueError:
				pass

	def filter(self, minlen=0, maxlen=None, minscore=0):
		"""Filter the matches as defined by the arguments
				minlen: minimum length of a match envelope
				maxlen: maximum length of a match envelope
				minscore: minimum score of a match
		"""
		to_remove = []
		for i,m in enumerate(self.matches):
			l = abs(m.env_to - m.env_from)
			if minlen:
				if l < minlen:
					to_remove.append(m)
			if maxlen:
				if l > maxlen:
					to_remove.append(m)
			if minscore:
				if m.score < minscore:
					to_remove.append(m)
		
		#remove all the matches which have failed
		for m in to_remove:
			try:
				self.matches.remove(m)
			except ValueError:
				pass
	
	def getFeatures(self, target, mode='hmm'):
		"""Return a list of feature found in a specific target"""
		ret = []
		for match in self.matches:
			if match.target == target:
				ret.append(match.asSeqFeature(mode=mode))
		return ret

	def annotate(self, targets, mode='hmm'):
		"""Annotate a target or list of targets given in target"""
		if isinstance(targets, SeqRecord):
			targets = [targets,]

		for t in targets:
			t.features.extend(self.getFeatures(t))

	class _minimatch:
		def __init__(self, m, mode='hmm'):
			self.match = m
			self.span = m.getFrameSpan(mode)

		def overlaps(self, m, dist=0):
			# either m starts in me OR m ends in me OR m contains me
			s1 = self.span[0] 
			e1 = self.span[1] + dist
			s2 = m.span[0]
			e2 = m.span[1] + dist

			return ( (s1 < s2 and e1 > s2) or #m starts in me
							 (s1 < e2 and e1 > e2) or #m ends in me
							 (s1 > s2 and e1 < e2)) #m contains me

		def dist(self, m):
			"""
				return the closest distance to m or None if the target or the frame
				are not equal			
			"""
			if (m.match.getTarget() != self.match.getTarget() or
					m.match.getFrame() != self.match.getFrame()):
				return None
			t = m.match.getTarget()
			if isinstance(t, basestring):
				return None

			if self.span[0] < m.span[0]:
				d1 = m.span[0] - self.span[1]
				d2 = len(t.seq) - (m.span[1] - self.span[0])
			else:
				d1 = self.span[0] - m.span[1]
				d2 = len(t.seq) - (self.span[1] - m.span[0])

			if abs(d1) < abs(d2):
				return d1
			return d2


	def _get_by_frame(self, mode='hmm'):
		r = dict()
		for t in self.targets:
			r[t.name] = [[], [], [], [], [], [], [],]

		for m in self.matches:
			r[m.target.name][m.frame].append(self._minimatch(m, mode))

		return r

	def __unicode__(self):
		def _seq(seq):
			if len(seq) < 60:
				return seq
			return "{} ... {}".format(seq[0:25], seq[-25:-1])

		ret = "hmmsearch:\n\t{}:\n".format(
				"Query" if len(self.hmm)==1 else "Queries")
		for i,hmm in enumerate(self.hmm,1):
			ret += "\t\t{}) {} [{} nodes],\n".format(i,hmm._name,len(hmm))
		ret += "\t{}:\n".format("Target" if len(self.targets)==1 else "Target")
		for i,t in enumerate(self.targets,1):
			ret += "\t\t{}) \'{}\': {}\n".format(i,t.name,_seq(t.seq))
		ret += "\nFound {:d} matches:\n".format(len(self.matches))
		
		if self.matches:
			for i,m in enumerate(self.matches,1):
				ret += "{})  {}\n".format(i, m)
		return ret

	def __str__(self):
		return unicode(self).encode('utf-8')

	def __repr__(self):
		return "[hmmsearch: {:d} matches]".format(len(self.matches))

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
					#('clustalo',   'Clustal Omega', 'www.clustal.org/omega/'),
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
			failed.append('Failed to find \'{}\' from package \'{}\'. ' 
					'Please install it from {}'.format(*b))

	if failed:
		raise ImportError('\n'.join(failed))

test_deps()
