from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
os.environ['CC'] = 'clang'
os.environ['CCFLAGS'] = "-fopenmp=libomp"
os.environ['CXX'] = 'clang++'
os.environ['CXXFLAGS'] = "-fopenmp=libomp"

ext_modules=[
    Extension("wrapPlanet",
    sources=["wrapPlanet.pyx"],
    language="c++"
    )]

setup(
    name="wrapPlanet",
    ext_modules=cythonize(ext_modules)
)
