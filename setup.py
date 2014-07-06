from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "searchtools",
    ext_modules = cythonize('searchtools.pyx')
)
