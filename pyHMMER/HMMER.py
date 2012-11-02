"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet, Alphabet

import cStringIO as StringIO
import hmmfile, tools, tempfile, subprocess, os, re
from subprocess import Popen, PIPE, STDOUT

import matchfile
		

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

		self.matches = matchfile.load(out_file, self.hmm, self.targets)

		out_file.close()
		hmm_file.close()
		target_file.close()

	def mindist(self, dist=0):
		"""
			Filter the results so that each match has a minimum number of dist
			symbols between it and the next match in the same frame
		"""
		#construct a list of ranges for each target frame
		l = dict()
		for t in self.targets:
			l[t.name] = list()


		#add each match to the list
		for m in self.matches:
			l[m.target.name].append(m)

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

	def filter(self, mindist=0, maxdist=None, minlen=0, maxlen=None, minscore=0):
		"""Filter the matches as defined by the arguments
				Nb. The order in which the filtering takes place is not guarenteed - 
						If you want to filter by multiple things (e.g. length then
						distance) you should make multiple calls to filter
			mindist: minimum distance to another match
			maxdist: maximum distance to another match
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
		
		
		#if we need to filter based distance
		if mindist or maxdist:
			#sort based on target and frame
			l = dict()
			for t in self.targets:
				l[t.name] = [ [],[],[],[],[],[],[],]
			for m in self.matches:
				if m.frame:
					l[m.target.name][m.getFrame()].append(self.minimatch(m))

			def group(lst):
				for i in range(0,len(lst)):
					yield (lst[i-1],lst[i],lst[(i+1)%len(lst)])

			print "l = %s" % l
			#for each target
			for frames in l.itervalues():
				print "frames = %s" % frames
				#for each frame in the target
				for frame in frames:
					print "frame = %s" % frame
					#sort the frame on ascending start position
					frame.sort(key=lambda m: m.span[0])
					#for each match in the frame
					for a,b,c in group(frame):
						#find the distance to the closest match
						dist = min(b.dist(a), b.dist(c))
						if mindist and dist < mindist:
							to_remove.append(b)
						if maxdist and dist > maxdist:
							to_remove.append(b)
		
		#remove all the matches which have failed
		for m in to_remove:
			try:
				self.matches.remove(m)
			except ValueError:
				pass
			
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
			failed.append('Failed to find \'{}\' from package \'{}\'. ' 
					'Please install it from {}'.format(*b))

	if failed:
		raise ImportError('\n'.join(failed))

test_deps()
