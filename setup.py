import os
from setuptools import setup, Extension

def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()


seqConvert = Extension('pyHMMER.seqConvert',
										sources = ['pyHMMER/seqConvert.c'])

setup(
		name = "pyHMMER",
		version = "0.0.1",
		author = "Haydn King",
		author_email = "hjking734@gmail.com",
		description = ("A python wrapper for HMMER"),
		license = "BSD",
		keywords = "wrapper HMMER HMM bioinformatics",
		packages=['pyHMMER', 'tests'],
		long_description=read('README'),
		test_suite='tests',
		ext_modules=[seqConvert,],
		)
