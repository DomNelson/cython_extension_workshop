from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

examples_extension = Extension(
    name="pymultiply", # Arbitrary (?) name
    sources=["pymultiply.pyx"], # Cython interface file
    libraries=["multiply"], # Name of compiled C library
    library_dirs=["."],
    include_dirs=[".",
                    np.get_include()] # If numpy is ever used
)
setup(
    name="PyMultiply",
    ext_modules=cythonize([examples_extension])
)
