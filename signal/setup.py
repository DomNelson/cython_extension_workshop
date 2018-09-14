from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

examples_extension = Extension(
    name="pysignal",
    sources=["pysignal.pyx"],
    libraries=["signal", "gsl"],
    library_dirs=["."],
    include_dirs=[np.get_include()]
)
setup(
    name="pysignal",
    ext_modules=cythonize([examples_extension], gdb_debug=True)
)
