from distutils.core import setup
from Cython.Build import cythonize

setup(name="wrapPlanet_layer",
      ext_modules=cythonize("wrapPlanet_layer.pyx"))

setup(name="wrapPlanet_auto",
      ext_modules=cythonize("wrapPlanet_auto.pyx"))
