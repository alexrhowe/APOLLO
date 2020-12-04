CXXFLAGS=-g
all: MakeHaze
MakeHaze : src/MakeHaze.o src/Atmosphere.haze.o src/constants.o
	g++ -o MakeHaze -g src/MakeHaze.o src/Atmosphere.haze.o src/constants.o
MakeHaze.o Atmosphere.haze.o constants.o : src/Atmosphere.haze.h src/constants.h
