from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

examples_extension = Extension(
    name="pymultiply", # This is the name you will 'import'
    sources=["pymultiply.pyx"], # Cython interface file
    libraries=["multiply"], # Links to library named 'lib[name].a
    library_dirs=[".", '/anaconda/pkgs/mpc-1.0.3-0/lib/'], # Could be omitted, since everything in one dir
    include_dirs=[np.get_include()] # Use numpy arrays so we need header files
)
setup(
    ext_modules=cythonize([examples_extension])
)
