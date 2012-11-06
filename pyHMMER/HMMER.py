"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet, Alphabet

import cStringIO as StringIO
import hmmfile, tools, tempfile, subprocess, os, re
from subprocess import Popen, PIPE, STDOUT

import matchfile
import sequtils	

class hmmsearch:
	"""Search for an HMM in a database and collect the results"""

	def __init__(self, hmm = None, targets = None):
		"""Initialise - search if hmm and targets have been provided"""
		if hmm and targets:
			self.search(hmm, targets)

	def search(self, hmm, targets):
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

		print "Loading HMMs and Targets..."

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
					try:
						#looks like we might have to convert
						ttargets = ttargets + tools.getSixFrameTranslation(t)
					except ValueError:
						#probably was a protein after all
						ttargets.append(t)

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

		print "Calling hmmsearch..."

		p = Popen(['hmmsearch', '--tformat', 'fasta', 
			'--domtblout', out_file.name, hmm_file.name, target_file.name,], 
				stdout=PIPE, stdin=PIPE, stderr=PIPE)
		out = p.communicate()


		#TODO: Check stdout for errors
		print "Reading in matches..."

		self.matches = matchfile.load(out_file, self.hmm, self.targets)

		print "Closing files..."

		out_file.close()
		hmm_file.close()
		target_file.close()

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

		return chains
					
	def extractProtein(self, chain):
		"""
			Extract the sequence of the chain, extending backwards to the start codon
			and forwards to the stop codon
			All matches in the chain must have the same target, be in the same 
			frame and the target must be a DNA alphabet
		"""
		if not chain:
			return Seq('')

		#sort the chain by start point
		chain.sort(key=lambda m: m.getTargetSpan()[0])

		target = chain[0].target
		frame = chain[0].frame
		start = chain[0].getTargetSpan()[0]
		end = chain[-1].getTargetSpan()[1]

		if frame > 0:
			step = +3
		elif frame < 0:
			step = -3

		prot = [start, end]
		#move back until we find a start codon
		while True:
			codon = str(target.seq[prot[0]:prot[0]+step]).upper()
			if ((frame > 0 and codon == 'ATG') or
					(frame < 0 and codon == 'CAT')):
				break
			prot[0] = prot[0] - step
			if prot[0] < 0 or prot[0] > len(target.seq):
				prot[0] = prot[0] + step
				break

		#move end forward until we meet a stop
		while True:
			codon = str(target.seq[prot[1]:prot[1]+step]).upper()
			if ((frame > 0 and codon in ['TAG', 'TAA', 'TGA',]) or
					(frame < 0 and codon in ['CTA', 'TTA', 'TCA',])):
				prot[1] = prot[1] + step
				break
			prot[1] = prot[1] + step
			if prot[1] < 0 or prot[1] > len(target.seq):
				prot[1] = prot[1] - step
				break

		ret = target.seq[prot[0]:prot[1]]
		if frame < 0:
			return ret.reverse_complement()
		return ret
		

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
