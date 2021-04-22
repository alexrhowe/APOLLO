from libcpp.vector cimport vector

cdef extern from "Planet_layer.h":
    cdef cppclass Planet:
        Planet(int modenum, vector[double] wavens, int cloudmod, int hazetype, vector[int] mollist) except +
        void setParams(vector[double] plparams, vector[double] abund, vector[double] tpprofile)
        double getTeff()
        vector[double] getSpectrum(int streams)
        void readopac(vector[int] mollist)
        void setWave(int npoints, double rxsec, vector[double] wavens, vector[double] abund)
        vector[double] getFlux(vector[double] wavens)
        vector[double] transFlux(double rs, vector[double] wavens)