"""Implements some of the core features of HMMER in an object orientated
	manner"""

from Bio import SeqIO, SeqRecord

import cStringIO as StringIO
import hmmfile, tempfile, subprocess, os
from subprocess import Popen, PIPE, STDOUT

class Match:
	"""Represents an HMM match"""
	fmt = ( 
		"{:20s} {:20s} {:4d} {:7g} {:7g} {:7g} {:9d} {:9d} {:9d} {:9d} {:9d} {:9d}"
		)
	hdr = (
		"{:20s} {:20s} {:4s} {:7s} {:7s} {:7s} {:9s} {:9s} {:9s} {:9s} {:9s} {:9s}"
		.format("Target", "Query", "#", "c-Evalue", "i-Evalue", "Score",
			"hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to")
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

	def __unicode__(self):
		return self.fmt.format(self.target.name, self.query.NAME, self.number,
				self.cEvalue, self.iEvalue, self.score, self.hmm_from, self.hmm_to,
				self.ali_from, self.ali_to, self.env_from, self.env_to)

	def __str__(self):
		return unicode(self).encode('utf-8')
	
		

class hmmsearch:
	"""Search for an HMM in a database and collect the results"""
	def __init__(self, hmm, targets):
		"""Perform the search
				hmm: a file name or an HMM object which has been loaded from a file
				targets: the sequences to search - one or more Bio.SeqRecord

				If the hmm performs searches on Amino Acids and and of the inputs are
				DNA sequences, 6-frame translations will be produced automatically
				Reverse translations (from Amino Acid to DNA) are not supported
		"""
		# Load the HMM
		#if the hmm is not an HMM object
		if not isinstance(hmm, hmmfile.HMM):
			hmm = hmmfile.read(hmm)

		#make sure targets is iterable
		if not hasattr(targets, '__iter__'):
			targets = [targets,]

		#TODO: Translate targets if necessary

		#write the HMM to a temporary file
		hmm_file = tempfile.NamedTemporaryFile()
		target_file = tempfile.NamedTemporaryFile()
		out_file = tempfile.NamedTemporaryFile()

		hmmfile.write(hmm, hmm_file)
		hmm_file.flush()

		SeqIO.write(targets, target_file, 'fasta')
		target_file.flush()

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

			match = Match()
			l = line.split()
			#target
			for t in targets:
				if t.name == l[0]:
					match.target = t
					break
			#query
			for q in hmm:
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

			self.matches.append(match)

		out_file.close()
		hmm_file.close()
		target_file.close()

	def discardOverlaps(self, choice_fn):
		"""Remove overlaping sequences from the results, using choice_fn to decide
			which result to remove and which to keep
		"""
		pass

	def expandResults(self):
		"""Expand the results to cover the entire HMM, even if they are only
			partial"""
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
