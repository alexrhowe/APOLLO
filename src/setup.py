from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
os.environ['CC'] = 'gcc'
os.environ['CXX'] = 'g++'

ext_modules=[
    Extension("wrapPlanet",
    sources=["wrapPlanet.pyx"],
    language="c++"
    )]

setup(
    name="wrapPlanet",
    ext_modules=cythonize(ext_modules)
)
