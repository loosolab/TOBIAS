import os
import sys
from setuptools import setup, Extension

# numpy and cython are guaranteed present at build time via
# pyproject.toml [build-system].requires (PEP 517 build isolation).
import numpy as np

# Compile from .pyx if Cython is available, otherwise from pre-generated .c
try:
	from Cython.Distutils import build_ext
	cmdclass = {'build_ext': build_ext}
	source_ext = "pyx"
except ImportError:
	cmdclass = {}
	source_ext = "c"
	if not os.path.exists("tobias/utils/ngs.c"):  # only happens when installing directly from source
		sys.exit("Cython is needed to compile TOBIAS.")

ext_modules = [Extension("tobias.utils.ngs", ["tobias/utils/ngs." + source_ext], include_dirs=[np.get_include()]),
			Extension("tobias.utils.sequences", ["tobias/utils/sequences." + source_ext], include_dirs=[np.get_include()]),
			Extension("tobias.utils.signals", ["tobias/utils/signals." + source_ext], include_dirs=[np.get_include()])]

# All package metadata now lives in pyproject.toml; only the C-extension build stays here.
setup(ext_modules=ext_modules, cmdclass=cmdclass)
