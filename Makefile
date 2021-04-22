CXXFLAGS=-g
all: MakeHaze
MakeHaze : MakeHaze.o Atmosphere.haze.o constants.o
	g++ -o MakeHaze -g MakeHaze.o Atmosphere.haze.o constants.o
MakeHaze.o Atmosphere.haze.o constants.o : Atmosphere.haze.h constants.h
