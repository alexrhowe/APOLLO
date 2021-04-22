# distutils: language = c++
# distutils: sources = [Planet_auto.cpp,Atmosphere.cpp,constants.cpp]

cdef class PyPlanet:
    cdef Planet *thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def MakePlanet(self, switches, wavens, wavenslo, mollist, opacdir, hires, lores):
        self.thisptr = new Planet(switches, wavens, wavenslo, mollist, opacdir, hires, lores)
    
    def get_Teff(self):
        return self.thisptr.getTeff()

    def set_Params(self, plparams, abund, tpprofile):
        return self.thisptr.setParams(plparams, abund, tpprofile)

    def get_Spectrum(self):
        return self.thisptr.getSpectrum()
    
    def readopac(self, mollist, wavens, table, opacdir):
        return self.thisptr.readopac(mollist, wavens, table, opacdir)

    def setWave(self, npoints, rxsec, wavens, abund):
        return self.thisptr.setWave(npoints, rxsec, wavens, abund)

    def getFlux(self, wavens, table):
        return self.thisptr.getFlux(wavens, table)

    def transFlux(self, rs, wavens, table):
        return self.thisptr.transFlux(rs, wavens, table)