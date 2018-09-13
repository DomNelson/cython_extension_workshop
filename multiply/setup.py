from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

examples_extension = Extension(
    name="pymultiply", # This is what you will call 'import' on
    sources=["pymultiply.pyx"], # Cython interface file
    libraries=["multiply"], # Links to library named 'lib[name].a
    library_dirs=["."], # Could be omitted, since everything in one dir
    include_dirs=[np.get_include()] # Need numpy header files
)
setup(
    ext_modules=cythonize([examples_extension])
)
