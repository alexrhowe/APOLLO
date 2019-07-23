# Functions that return thruput for individual spectroscopic modes
# for the JWST pipeline.

import numpy as np
from scipy import interpolate

class Filter():
    def __init__(self):
        self.delwave = np.zeros(2048)
        self.calwave = np.zeros(2048)
        self.resolve = 1.0
        self.slit = 1.0
        self.omeg_back = 1.e-6
        self.scale = 1.0
        self.effic = 1.0
        self.npix1 = 1
        self.npix2 = 2
        self.npix3 = 3
        self.npix4 = 4

    def NIRSpec(self,fnum):
        fnum = int(fnum)
        fFilter = open('FilterList.dat','r')
        flines = fFilter.readlines()
        flength = len(flines)
        filter = flines[fnum].split()[0]
        fmin = float(flines[fnum].split()[1])
        fmax = float(flines[fnum].split()[2])

        if fnum<4:
            self.delwave = (fmax-fmin)/1399.0
            self.calwave = np.zeros(1400)
        if 4<=fnum<8:
            self.delwave = (fmax-fmin)/511.0
            self.calwave = np.zeros(512)
        if fnum==8:
            self.delwave = (fmax-fmin)/255.0
            self.calwave = np.zeros(256)
        
        for i in range(0,len(self.calwave)):
            self.calwave[i] = fmin + self.delwave*i
        
        fin = open(filter,'r')
        lines = fin.readlines()
        length = len(lines)
        microns = np.zeros(length)
        thruput = np.zeros(length)
        for i in range(0,length):
            microns[i] = lines[i].split()[0]
            thruput[i] = lines[i].split()[1]
        
        f = interpolate.interp1d(microns,thruput,kind='linear')
        xthruput = f(self.calwave)

        theta = 10.08  # half-angle of beam to the detector in degrees
        self.scale = 0.065  # arcsec per pixel
        self.slit = -999.   # IFU = slitless
        self.resolve = 2.0  # 2 column resolution
        
        ngr = 4
        self.effic = (ngr-1.)/(ngr+1.)
        return xthruput
    
    def NIRSpecG395H(self):
        self.delwave = (4.9-2.9)/2047.
        self.calwave = np.zeros(2048)
        for i in range(0,len(self.calwave)):
            self.calwave[i] = 2.9 + self.delwave*i

        # read in *.txt response function (photons to electrons) to xmicrons and xthruput       
        fin = open('Filters/NIRSpec_G395H.dat','r')
        lines = fin.readlines()
        length = len(lines)
        microns = np.zeros(length)
        thruput = np.zeros(length)
        for i in range(0,length):
            microns[i] = lines[i].split()[0]
            thruput[i] = lines[i].split()[1]

        f = interpolate.interp1d(microns,thruput,kind='linear')
        xthruput = f(self.calwave)

        theta = 10.08       # half-angle of beam to the detector in degrees
        self.scale = 0.065  # arcsec per pixel
        self.slit = 24.6    # entrance slit size in columns, this is the 1.6x1.6 arcsec case
        self.resolve = 2.0  # 2 column resolution
        
        ngr = 4
        self.effic = (ngr-1.)/(ngr+1.)
        return xthruput

    def NIRCamF322W2(self):
        self.effic = 0.8
        npix = 1675
        self.delwave = (5.0-2.4)/(npix-1.)
        self.calwave = np.zeros(npix)
        for i in range(0,len(self.calwave)):
            self.calwave[i] = 2.4 + self.delwave*i

        theta = 0.0
        self.scale = 0.065
        self.slit = -999.   # no slit
        self.resolve = 2.0
        tel = 0.94     # telescope efficiency
        optics = 0.541 # instrument optics + filter
        
        gwave = np.zeros(29)  # grism wavelengths
        for i in range(0,len(gwave)):
            gwave[i] = 2.30 + 0.1*i
            # grism thruput
        grespon = [0.24, 0.24, 0.306,0.374,0.437,0.494,0.543,0.585,0.62, 0.648,0.67, 0.686,0.696,0.702,0.705,0.703,0.702,0.694,0.685,0.674,0.661,0.649,0.636,0.621,0.609,0.593,0.579,0.566,0.566]
            # total instrument thruput
        trespon2 = [0.10,0.10,0.13,0.16,0.19,0.22,0.24,0.26,0.29,0.30,0.32,0.34,0.36,0.37,0.37,0.37,0.36,0.35,0.35,0.34,0.32,0.30,0.28,0.26,0.24,0.23,0.20,0.10,0.10]
            # detector QE
        qe = 0.934871364 + 0.051540672*gwave - 0.281664062*gwave*gwave + 0.243866616*pow(gwave,3) - 0.086009299*pow(gwave,4) + 0.014508849*pow(gwave,5) - 0.001000108*pow(gwave,6)

        qe = qe*0.88
        trespon = np.zeros(27)
        for i in range(0,27): trespon = tel*optics*grespon[i]*qe[i]

        f = interpolate.interp1d(gwave,trespon2,kind='linear')
        xthruput = f(self.calwave)
            
        return xthruput

    def MIRILRS(self,slitin):
        self.effic = 0.8
        npix = 320   # number of detector pixels used by the spectrum, from a Kendrew figure
        self.delwave = (14.0-5.0)/(npix-1.)
        self.calwave = np.zeros(npix)
        for i in range(0,len(self.calwave)):
            self.calwave[i] = 5.0 + self.delwave*i

        theta = 0.0
        self.scale = 0.11
        if(slitin==True):
            self.slit = 4.6
            self.resolve = 4.6
            optics = 0.6
        if(slitin==False):
            self.slit = -999.0
            self.resolve = 2.7
            optics = 1.0

        tel = 0.94
        gwave = np.zeros(19)
        grespon = [0.1049,0.1409,0.2292,0.2524,0.2945,0.2589,0.2631,0.3123,0.3047,0.2805,0.2602,0.2352,0.2156,0.1550,0.1132,0.0771,0.0502,0.0251,0.0251]
        trespon = np.zeros(19)
        for i in range(0,len(gwave)):
            gwave[i] = 5.0 + 0.5*i
            trespon[i] = tel*optics*grespon[i]

        f = interpolate.interp1d(gwave,trespon,kind='linear')
        xthruput = f(self.calwave)
        return xthruput

    def NIRISS(self):
        # Alternative throughput function
        #fin = open('Filters/NIRISS_Throughput.dat','r')
        fin = open('Filters/NIRISS_G700XD.dat','r')
        lines = fin.readlines()
        nwave = np.zeros(len(lines))
        nresponse = np.zeros(len(lines))
        for i in range(0,len(lines)):
            nwave[i] = lines[i].split()[0]
            nresponse[i] = lines[i].split()[1]

        self.delwave = (2.5-0.602)/2047.0
        self.calwave = np.zeros(2048)
        for i in range(0,len(self.calwave)):
            self.calwave[i] = 0.602 + self.delwave*i

        theta = 10.08
        self.scale = 0.0658
        self.slit = -999.0
        self.resolv = 2.0
        self.effic = 0.8

        f = interpolate.interp1d(nwave,nresponse,kind='linear')
        xthruput = f(self.calwave)
        return xthruput

    def MIRIMRS(self,channel):
        self.effic = 0.8  # arbitrary for now
        
        if channel == 'A':
            lam1 = 4.87
            lam2 = 5.82
            r1 = 3320.0
            r2 = 3710.0
        if channel == 'B':            
            lam1 = 5.62
            lam2 = 6.73
            r1 = 3190.0
            r2 = 3750.0
        if channel == 'C':            
            lam1 = 6.49
            lam2 = 7.76
            r1 = 3100.0
            r2 = 3610.0

        w1 = 1.405
        w2 = 1.791
        w = (w1+w2)/2.0   # slice width in pixels (width = dispersion direction)
        r = (r1+r2)/2.0   # approximate the resolving power at mid-channel
        self.npix1 = np.modf(w*(lam2-lam1)*r + 0.5)[1]
        self.delwave = (lam2-lam1)/self.npix1
        calwave1 = np.zeros(self.npix1)
        for i in range(0,len(calwave1)):
            calwave1[i] = lam1 + self.delwave*i
        cwave11 = np.zeros(self.npix1)
        cwave12 = np.zeros(self.npix1)  # wavelengths at the boundaries of the pixels

        for i in range(0,self.npix1-2):
            cwave12[i] = (calwave1[i] + calwave1[i+1])/2.0
        for i in range(0,self.npix1-1):
            cwave11[i] = (calwave1[i-1] + calwave1[i])/2.0
        cwave11[0] = calwave1[0] - (cwave11[1] - calwave1[0])
        cwave12[self.npix1-1] = calwave1[self.npix1-1] + (calwave1[self.npix1-1] - cwave12[self.npix1-2])
        
        if channel == 'A':
            lam1 = 7.45
            lam2 = 8.90
            r1 = 2990.0
            r2 = 3110.0
        if channel == 'B':
            lam1 = 8.61
            lam2 = 10.28
            r1 = 2750.0
            r2 = 3170.0
        if channel == 'C':
            lam1 = 9.91
            lam2 = 11.87
            r1 = 2860.0
            r2 = 3300.0

        w1 = 1.452
        w2 = 1.821
        w = (w1+w2)/2.0   # slice width in pixels (width = dispersion direction)
        r = (r1+r2)/2.0   # approximate the resolving power at mid-channel
        self.npix2 = np.modf(w*(lam2-lam1)*r + 0.5)[1]
        self.delwave = (lam2-lam1)/self.npix2
        calwave2 = np.zeros(self.npix2)
        for i in range(0,len(calwave2)):
            calwave2[i] = lam1 + self.delwave*i
        cwave21 = np.zeros(self.npix2)
        cwave22 = np.zeros(self.npix2)  # wavelengths at the boundaries of the pixels

        for i in range(0,self.npix2-2):
            cwave22[i] = (calwave2[i] + calwave2[i+1])/2.0
        for i in range(0,self.npix2-1):
            cwave21[i] = (calwave2[i-1] + calwave2[i])/2.0
        cwave21[0] = calwave2[0] - (cwave21[1] - calwave2[0])
        cwave22[self.npix2-1] = calwave2[self.npix2-1] + (calwave2[self.npix2-1] - cwave22[self.npix2-2])
        
        if channel == 'A':
            lam1 = 11.47
            lam2 = 13.67
            r1 = 2530.0
            r2 = 2880.0
        if channel == 'B':
            lam1 = 13.25
            lam2 = 15.80
            r1 = 1790.0
            r2 = 2640.0
        if channel == 'C':
            lam1 = 15.30
            lam2 = 18.24
            r1 = 1980.0
            r2 = 2790.0

        w1 = 1.629
        w2 = 2.043
        w = (w1+w2)/2.0   # slice width in pixels (width = dispersion direction)
        r = (r1+r2)/2.0   # approximate the resolving power at mid-channel
        self.npix3 = np.modf(w*(lam2-lam1)*r + 0.5)[1]
        self.delwave = (lam2-lam1)/self.npix3
        calwave3 = np.zeros(self.npix3)
        for i in range(0,len(calwave3)):
            calwave3[i] = lam1 + self.delwave*i
        cwave31 = np.zeros(self.npix3)
        cwave32 = np.zeros(self.npix3)  # wavelengths at the boundaries of the pixels

        for i in range(0,self.npix3-2):
            cwave32[i] = (calwave3[i] + calwave3[i+1])/2.0
        for i in range(0,self.npix3-1):
            cwave31[i] = (calwave3[i-1] + calwave3[i])/2.0
        cwave31[0] = calwave3[0] - (cwave31[1] - calwave3[0])
        cwave32[self.npix3-1] = calwave3[self.npix3-1] + (calwave3[self.npix3-1] - cwave32[self.npix3-2])
        
        if channel == 'A':
            lam1 = 17.54
            lam2 = 21.10
            r1 = 1460.0
            r2 = 1930.0
        if channel == 'B':
            lam1 = 20.44
            lam2 = 24.72
            r1 = 1680.0
            r2 = 1770.0
        if channel == 'C':
            lam1 = 23.84
            lam2 = 28.82
            r1 = 1630.0
            r2 = 1330.0

        w1 = 2.253
        w2 = 2.824
        w = (w1+w2)/2.0   # slice width in pixels (width = dispersion direction)
        r = (r1+r2)/2.0   # approximate the resolving power at mid-channel
        self.npix4 = np.modf(w*(lam2-lam1)*r + 0.5)[1]
        self.delwave = (lam2-lam1)/self.npix4
        calwave4 = np.zeros(self.npix4)
        for i in range(0,len(calwave4)):
            calwave4[i] = lam1 + self.delwave*i
        cwave41 = np.zeros(self.npix4)
        cwave42 = np.zeros(self.npix4)  # wavelengths at the boundaries of the pixels

        for i in range(0,self.npix4-2):
            cwave42[i] = (calwave4[i] + calwave4[i+1])/2.0
        for i in range(0,self.npix4-1):
            cwave41[i] = (calwave4[i-1] + calwave4[i])/2.0
        cwave41[0] = calwave4[0] - (cwave41[1] - calwave4[0])
        cwave42[self.npix4-1] = calwave4[self.npix4-1] + (calwave4[self.npix4-1] - cwave42[self.npix4-2])

        self.scale = 0.20
        self.slit = 1.0
        self.resolv = 2.0
        
        xthruput = np.zeros(self.npix1+self.npix2+self.npix3+self.npix4)
        for i in range(0,self.npix1):
            if channel == 'A': xthruput[i] = 0.117
            if channel == 'B': xthruput[i] = 0.117
            if channel == 'C': xthruput[i] = 0.138
        for i in range(0,self.npix2):
            if channel == 'A': xthruput[self.npix1+i] = 0.114
            if channel == 'B': xthruput[self.npix1+i] = 0.130
            if channel == 'C': xthruput[self.npix1+i] = 0.134
        for i in range(0,self.npix3):
            if channel == 'A': xthruput[self.npix1+self.npix2+i] = 0.115
            if channel == 'B': xthruput[self.npix1+self.npix2+i] = 0.108
            if channel == 'C': xthruput[self.npix1+self.npix2+i] = 0.114
        for i in range(0,self.npix4):
            if channel == 'A': xthruput[self.npix1+self.npix2+self.npix3+i] = 0.028
            if channel == 'B': xthruput[self.npix1+self.npix2+self.npix3+i] = 0.020
            if channel == 'C': xthruput[self.npix1+self.npix2+self.npix3+i] = 0.011

        self.calwave = np.zeros(self.npix1+self.npix2+self.npix3+self.npix4)
        calwave12 = np.concatenate((calwave1,calwave2),axis=0)
        calwave123 = np.concatenate((calwave12,calwave3),axis=0)
        self.calwave = np.concatenate((calwave123,calwave4),axis=0)

        return xthruput
