from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
os.environ['CC'] = 'gcc-mp-9'
os.environ['CXX'] = 'g++-mp-9'

ext_modules=[
    Extension("wrapPlanet_layer",
    sources=["wrapPlanet_layer.pyx"],
    language="c++"
    )]

setup(
    name="wrapPlanet_layer",
    ext_modules=cythonize(ext_modules)
)

ext_modules=[
    Extension("wrapPlanet_auto",
    sources=["wrapPlanet_auto.pyx"],
    language="c++"
    )]

setup(
    name="wrapPlanet_auto",
    ext_modules=cythonize(ext_modules)
)
