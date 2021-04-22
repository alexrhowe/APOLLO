from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Planet_auto.h":
    cdef cppclass Planet:
        Planet(int modenum, vector[double] wavens, int cloudmod, int hazetype, int wrange, vector[int] mollist, string opacdir) except +
        void setParams(vector[double] plparams, vector[double] abund, vector[double] tpprofile)
        vector[double] getSpectrum(int streams)
        void readopac(vector[int] mollist, vector[double] wavens, string opacdir)
        void setWave(int npoints, double rxsec, vector[double] wavens, vector[double] abund)
        vector[double] getFlux(vector[double] wavens)
        vector[double] transFlux(double rs, vector[double] wavens)