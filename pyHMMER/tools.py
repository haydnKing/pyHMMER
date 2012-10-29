"""Some helpful functions"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

import seqConvert

def getSixFrameTranslation(rec):
	"""Return the six frame translations of SeqRecord rec"""
	frames = ["+1", "+2", "+3", "-1", "-2", "-3",]
	c = seqConvert.sixFrameTranslation(str(rec.seq))
	r = []
	for i,s in enumerate(c):
		r.append(SeqRecord(Seq(s, generic_protein),
			name=rec.name, id=rec.id,
			description= "%s; frame: %s" % (rec.description, frames[i]), 
			))
	return r

