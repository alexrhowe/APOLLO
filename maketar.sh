#!/bin/bash

tar -cvf Apollo0115.tar --exclude=src/w*cpp Apollo.py Apollo.pbs Apollo.run FilterList.dat Filters/*.dat Makefile Plot.Apollo.py Plot.Converge.py Plot.Haze.py Plot.Opac.py README.txt examples/*.dat maketar.sh modelspectra/example.spectrum.dat modelspectra/example.binned.dat modelspectra/example.fullres.dat plots/example.TP.png plots/example.basic.png plots/example.fullres.fit.png plots/example.binned.fit.png plots/example.gases.png samples/example.dat src/*.h src/*.cpp src/*.py src/*.pxd src/*.pyx
