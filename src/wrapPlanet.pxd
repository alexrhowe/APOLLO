from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Planet.h":
    cdef cppclass Planet:
        Planet(vector[int] switches, vector[double] wavens, vector[double] wavenslo, vector[int] mollist, string opacname, string hir, string lor) except +
        void setParams(vector[double] plparams, vector[double] abund, vector[double] rxsecs, vector[double] tpprofile)
        double getTeff()
        vector[double] getSpectrum()
        vector[double] getClearSpectrum()
        void readopac(vector[int] mollist, vector[double] wavens, string table, string opacdir)
        void setWave(int npoints, vector[double] rxsecs, vector[double] wavens, vector[double] abund)
        vector[vector[double]] getContribution()
        vector[vector[double]] getCloudContribution()
        vector[vector[double]] getGasContribution()
        vector[vector[vector[double]]] getSpeciesContribution()
        vector[double] getFlux(vector[double] wavens, string table)
        vector[double] transFlux(double rs, vector[double] wavens, string table)