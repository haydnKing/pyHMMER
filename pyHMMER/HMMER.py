"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet, Alphabet, IUPAC

import cStringIO as StringIO
import hmmfile, tools, tempfile, subprocess, os, re
from subprocess import Popen, PIPE, STDOUT

import matchfile
import sequtils	

class wrap(object):
	"""A class which wraps anothe class, overriding specific values"""
	def __init__(self, wrapped, overrides):
		self.wrapped = wrapped
		self.overrides = overrides

	def __getattr__(self, item):
		if self.overrides.has_key(item):
			return self.overrides[item]
		return getattr(self.wrapped, item)

def wrap_seqrecords(records):
	return [wrap(r, {'id': str(i), 'alpha': sequtils.seq_type(str(r.seq))}) 
				for i,r in enumerate(records)]

def wrap_hmms(hmms):
	return [wrap(h, {'name': str(i), 'alpha': h.alph.upper(),}) 
					for i,h in enumerate(hmms)]

class jackhmmer:
	"""Iteratively seach a protein database with a protein sequence

	Attributes:
		- matches: matches found by jackhmmer (matchfile.Match)
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

		args = []
		for k,v in kwargs.iteritems():
			if k in self.ARGS:
				args += ["--{}".format(k), v]
			elif k in self.SWITCHES:
				args += ['--{}'.format(k)]
			else:
				raise ValueError("Unknown jackhmmer argument \'{}\'".format(k))

		#apply unique ids to the targets
		self.seq = wrap_seqrecords(seq)
		self.seqdb=wrap_seqrecords(seqdb)

		seq_file = tempfile.NamedTemporaryFile()
		seqdb_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()

		SeqIO.write(self.seq, seq_file, 'fasta')
		SeqIO.write(self.seqdb, seqdb_file, 'fasta')
		seq_file.flush()
		seqdb_file.flush()

		p = Popen(['jackhmmer', '--qformat', 'fasta', '--tformat', 'fasta', 
			'--domtblout', out_file.name,] + args + [seq_file.name, seqdb_file.name,], 
				stdout=PIPE, stdin=PIPE, stderr=PIPE)
		out = p.communicate()

		self.matches = matchfile.load(out_file, self.seq, self.seqdb)

		if verbose:
			print "Closing files..."

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


class hmmsearch:
	"""Search for an HMM in a database and collect the results"""

	def __init__(self, hmm = None, targets = None, **kwargs):
		"""Initialise - search if hmm and targets have been provided"""
		if hmm and targets:
			self.search(hmm, targets, **kwargs)

	def search(self, hmm, targets, verbose=False):
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
		
		if verbose:
			print "Loading HMMs and Targets..."

		#apply unique ids
		self.targets = wrap_seqrecords(self.targets) 
		self.hmm = wrap_hmms(self.hmm)

		ttargets = []
		hmm_alpha = self.hmm[0].alph.upper()
		for h in self.hmm:
			if h.alph.upper() != hmm_alpha:
				raise ValueError("The HMMs don't all have the same alphabet")

		#Translate targets if necessary
		for t in self.targets:
			t_alpha = sequtils.seq_type(str(t.seq))
			
			#if target and hmms have same alphabet
			if hmm_alpha == t_alpha:
					ttargets.append(t)
			#else try translating
			else:
				if hmm_alpha == "AMINO" and t_alpha == "DNA":
					try:
						#looks like we might have to convert
						ttargets = ttargets + tools.getSixFrameTranslation(t)
					except ValueError:
						#probably was a protein after all
						ttargets.append(t)
				else:
					raise ValueError('Unknown Translation \'{}\' to \'{}\''
							.format(t_alpha, hmm_alpha))

		if verbose:
			print "Writing Temporary Files"

		#write the HMM to a temporary file
		hmm_file = tempfile.NamedTemporaryFile()
		target_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()

		hmmfile.write(self.hmm, hmm_file)
		hmm_file.flush()

		SeqIO.write(ttargets, target_file, 'fasta')
		target_file.flush()
		del ttargets

		if verbose:
			print "Calling hmmsearch..."

		p = Popen(['hmmsearch', '--tformat', 'fasta', 
			'--domtblout', out_file.name, hmm_file.name, target_file.name,], 
				stdout=PIPE, stdin=PIPE, stderr=PIPE)
		out = p.communicate()


		#TODO: Check stdout for errors
		if verbose:
			print "Reading in matches..."

		self.matches = matchfile.load(out_file, self.hmm, self.targets)

		if verbose:
			print "Closing files..."

		out_file.close()
		hmm_file.close()
		target_file.close()

		if verbose:
			print "Found {} matches in {} targets".format(len(self.matches), 
					len(self.targets))

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
	
	def chain(self, mingap=0, maxgap=0, minlen=1):
		"""
			Extract a chain of in-frame matches
				mingap: minimum gap between matches (0)
				maxgap: maximum gap between matches (0)
				minlen: minimum chain length for reporting (1)
		"""
		l = self._get_by_frame('hmm')

		chains = []

		def add_chain(chain):
			if len(chain) >= minlen:
				chains.append([m.match for m in chain])

		#for each target
		for target in l.itervalues():
			for frame in target:
				#ignore empty frames
				if not frame:
					continue
				chain = [frame[0],]
				for m in frame[1:]:
					dist = chain[-1].dist(m)
					#add m to the chain if it's within range
					if dist >= mingap and dist <= maxgap:
						chain.append(m)
					#or if m is too far away, start a new chain
					elif dist > maxgap:
						#add chain to chains if it's long enough
						add_chain(chain)
						#start a new chain
						chain = [m,]
					#ignore matches that are too close

				add_chain(chain)

		print "Found {} chains".format(len(chains))

		return chains

	def annotate(self, target, mode='hmm'):
		"""Annotate the seqRecord given by target"""
		target = target or self.targets[0]

		for match in self.matches:
			if match.target == target:
				target.features.append(match.asSeqFeature(mode=mode))

		return target

	def getProteins(self, **kwargs):
		#prepare the arguments
		chain_args = kwargs.copy()
		for p in ['max_5_prime', 'max_3_prime', 'mode']:
			try:
				chain_args.pop(p)
			except KeyError:
				pass
		prot_args = kwargs.copy()
		for p in ['mingap', 'maxgap', 'minlen',]:
			try:
				prot_args.pop(p)
			except KeyError:
				pass
		

		chains = self.chain(**chain_args)
		ret = []
		for c in chains:
			p = self.getProtein(c, **prot_args)
			if p:
				ret.append(p)

		print "Found {} proteins".format(len(ret))
		return ret
					
	def getProtein(self, chain, mode='hmm', max_5_prime=None, max_3_prime=None):
		"""
			Extract the sequence of the chain, extending backwards to the start codon
			and forwards to the stop codon
			All matches in the chain must have the same target, be in the same 
			frame and the target must be a DNA alphabet
		"""
		if not chain:
			return Seq('')

		#sort the chain by start point
		chain.sort(key=lambda m: m.getFrameSpan(mode)[0])

		target = chain[0].target
		query = chain[0].query
		frame = chain[0].frame
		start = chain[0].getTargetSpan(mode)[0]
		end = chain[-1].getTargetSpan(mode)[1]

		if frame == 0:
			frame = 1
		if frame > 0:
			step = +3
			fwd = +1
		elif frame < 0:
			step = -3
			fwd = -1

		prot = [start, end]
		#move back until we find a start codon
		while True:
			if frame > 0:
				codon = str(target.seq[prot[0]:prot[0]+step]).upper()
			else:
				codon = str(
						target.seq[prot[0]+step:prot[0]].reverse_complement()).upper()
			if codon == 'ATG':
				break
			prot[0] = prot[0] - step
			if prot[0] < 0 or prot[0] > len(target.seq):
				prot[0] = prot[0] + step
				break

		#move end forward until we meet a stop
		while True:
			if frame > 0:
				codon = str(target.seq[prot[1]:prot[1]+step]).upper()
			else:
				codon = str(
						target.seq[prot[1]+step:prot[1]].reverse_complement()).upper()
			if codon in ['TAG', 'TAA', 'TGA',]:
				prot[1] = prot[1] + step
				break
			prot[1] = prot[1] + step
			if prot[1] < 0 or prot[1] > len(target.seq):
				prot[1] = prot[1] - step
				break

		#check if the 5' and 3' ends are within the limits
		if max_5_prime and max_5_prime < abs(start - prot[0]):
			return None
		if max_3_prime and max_3_prime < abs(end - prot[1]):
			return None

		if frame < 0:
			seq = target.seq[prot[1]:prot[0]].reverse_complement()
			s_loc = FeatureLocation(prot[0], prot[1], strand=-1)
		else:
			seq = target.seq[prot[0]:prot[1]]
			s_loc = FeatureLocation(prot[0], prot[1], strand=1)

		seq_type = sequtils.seq_type(str(seq))
		if seq_type == 'DNA':
			seq.alphabet = IUPAC.ambiguous_dna
		elif seq_type == 'RNA':
			seq.alphabet = IUPAC.ambiguous_rna
		elif seq_type == 'AMINO':
			seq.alphabet = IUPAC.protein

		#Create the features
		feats = []
		for m in chain:
			feats.append(m.asSeqFeature(mode=mode, offset=prot[0]))
			
		return SeqRecord(seq, name=query.name, 
				description="{} containing protein".format(query.name), 
				features=feats,
				annotations={'source': target,
										 'sourceloc': s_loc,})



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
		ret = "Found {:d} matches\n".format(len(self.matches))
		if self.matches:
			for m in self.matches:
				ret += str(m) + '\n'
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
			failed.append('Failed to find \'{}\' from package \'{}\'. ' 
					'Please install it from {}'.format(*b))

	if failed:
		raise ImportError('\n'.join(failed))

test_deps()
