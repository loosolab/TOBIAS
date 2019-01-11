from setuptools import setup, Extension
import numpy as np

def readme():
    with open('README.md') as f:
        return f.read()

ext_modules = [Extension("tobias.utils.ngs", ["tobias/utils/ngs.pyx"], include_dirs=[np.get_include()]),
              Extension("tobias.utils.sequences", ["tobias/utils/sequences.pyx"], include_dirs=[np.get_include()]),
               Extension("tobias.utils.signals", ["tobias/utils/signals.pyx"], include_dirs=[np.get_include()])]

setup(name='tobias',
      version='0.1',
      description='Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal',
      long_description=readme(),
      url='https://github.molgen.mpg.de/loosolab/TOBIAS',
      author='Mette Bentsen',
      author_email='mette.bentsen@mpi-bn.mpg.de',
      license='MIT',
      packages=['tobias', 'tobias.footprinting', 'tobias.utils', 'tobias.plotting', 'tobias.motifs'],
      entry_points = {
        'console_scripts': ['TOBIAS=tobias.TOBIAS:main']
      },
      install_requires=[
        'setuptools_cython',
        'numpy',
        'scipy',
        'pyBigWig',
        'pysam',
        'pybedtools',
        'matplotlib>=2',
        'scikit-learn',
        'pandas',
        'pypdf2',
        'xlsxwriter',
        'adjustText',
      ],
      #dependency_links=['https://github.com/jhkorhonen/MOODS/tarball/master'],
      classifiers = [
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3'
      ],
      zip_safe=False,
      include_package_data=True,
      ext_modules = ext_modules,
      scripts=["tobias/utils/peak_annotation.sh"]
      )
