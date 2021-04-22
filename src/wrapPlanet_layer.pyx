# distutils: language = c++
# distutils: sources = [Planet_layer.cpp,Atmosphere.cpp,constants.cpp]

cdef class PyPlanet:
    cdef Planet *thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def MakePlanet(self, modenum, wavens, cloudnum, hazetype, wrange, mollist, opacdir):
        self.thisptr = new Planet(modenum, wavens, cloudnum, hazetype, wrange, mollist, opacdir)

    def set_Params(self, plparams, abund, tpprofile):
        return self.thisptr.setParams(plparams, abund, tpprofile)

    def get_Spectrum(self, streams):
        return self.thisptr.getSpectrum(streams)

    def readopac(self, mollist, wavens, opacdir):
        return self.thisptr.readopac(mollist, wavens, opacdir)

    def setWave(self, npoints, rxsec, wavens, abund):
        return self.thisptr.setWave(npoints, rxsec, wavens, abund)

    def getFlux(self, wavens):
        return self.thisptr.getFlux(wavens)

    def transFlux(self, rs, wavens):
        return self.thisptr.transFlux(rs, wavens)