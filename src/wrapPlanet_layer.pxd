from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Planet_layer.h":
    cdef cppclass Planet:
        Planet(vector[int] switches, vector[double] wavens, vector[double] wavenslo, vector[int] mollist, string opacname, string hir, string lor) except +
        void setParams(vector[double] plparams, vector[double] abund, vector[double] tpprofile)
        double getTeff()
        vector[double] getSpectrum()
        void readopac(vector[int] mollist, vector[double] wavens, string table, string opacdir)
        void setWave(int npoints, double rxsec, vector[double] wavens, vector[double] abund)
        vector[double] getFlux(vector[double] wavens, string table)
        vector[double] transFlux(double rs, vector[double] wavens, string table)