import os
import sys
import re
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

included_dirs = []
try:
	import numpy as np
	included_dirs.append(np.get_include())
except:
	pass

#Test if numpy is installed
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy as np
       	self.include_dirs.append(np.get_include())
       	included_dirs.append(np.get_include())

#Add cython modules depending on the availability of cython
try:
	from Cython.Distutils import build_ext
except ImportError:
	use_cython = False
else:
	use_cython = True

if use_cython:
	ext_modules = [Extension("tobias.utils.ngs", ["tobias/utils/ngs.pyx"], include_dirs=included_dirs), #,sinclude[np.get_include()]),
					Extension("tobias.utils.sequences", ["tobias/utils/sequences.pyx"], include_dirs=included_dirs), #include_dirs=[np.get_include()]),
					Extension("tobias.utils.signals", ["tobias/utils/signals.pyx"], include_dirs=included_dirs)] #, include_dirs=[np.get_include()])]

else:
	ext_modules = [Extension("tobias.utils.ngs", ["tobias/utils/ngs.c"], include_dirs=included_dirs), #, include_dirs=[np.get_include()]),
					Extension("tobias.utils.sequences", ["tobias/utils/sequences.c"], include_dirs=included_dirs), #, include_dirs=[np.get_include()]),
					Extension("tobias.utils.signals", ["tobias/utils/signals.c"], include_dirs=included_dirs)] #, include_dirs=[np.get_include()])]

cmdclass = {'build_ext': build_ext}

#Path of setup file to establish version
setupdir = os.path.abspath(os.path.dirname(__file__))

def find_version(init_file):
	version_file = open(init_file).read()
	version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
	if version_match:
		return version_match.group(1)
	else:
		raise RuntimeError("Unable to find version string.")

def readme():
	with open('README.md') as f:
		return f.read()

setup(name='tobias',
		version=find_version(os.path.join(setupdir, "tobias", "__init__.py")),	#get version from __init__.py
		description='Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal',
		long_description=readme(),
		url='https://github.molgen.mpg.de/loosolab/TOBIAS',
		author='Mette Bentsen',
		author_email='mette.bentsen@mpi-bn.mpg.de',
		license='MIT',
		packages=['tobias', 'tobias.footprinting', 'tobias.plotting', 'tobias.motifs', 'tobias.misc', 'tobias.utils'],
		entry_points={
			'console_scripts': ['TOBIAS=tobias.TOBIAS:main']
		},
		ext_modules=ext_modules,
		cmdclass=cmdclass,
		setup_requires=["numpy"],
		include_dirs=included_dirs,
		#dependency_links=['https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz#egg=MOODS-python-1.9.3'],	
		install_requires=[
			'numpy',
			'scipy',
			'pysam',
			'pybedtools',
			'matplotlib>=2',
			'scikit-learn',
			'pandas',
			'pypdf2',
			'xlsxwriter',
			'adjustText',
			'pyBigWig',
		],

		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Intended Audience :: Science/Research',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			'Programming Language :: Python :: 3'
		],
		zip_safe=True,
		)
