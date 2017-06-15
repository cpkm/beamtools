#from distutils.core import setup
from setuptools import setup
setup(
  name = 'beamtools',
  packages = ['beamtools'], # this must be the same as the name above
  version = '0.3',
  platforms=['any'],
  description = 'A suite of tools for the analysis of optical beams',
  author = 'Kyle Manchee',
  author_email = 'cpkmanchee@gmail.com',
  url = 'https://github.com/cpkm/beamtools',
  download_url = 'https://github.com/cpkm/beamtools/archive/v0.3.tar.gz',
  keywords = ['optics', 'analysis', 'beams', 'spectrum'],
  classifiers = [],
)
