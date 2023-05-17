from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import matplotlib.pyplot as plt
from . import Filter
import sys
from matplotlib import rc
from scipy import interpolate
rc('text',usetex=True)

def addNoise(mode,obstype,mod_wave,mod_flux,noise_params,starspec = '',fname = ''):
    # mod_flux = emission from the target in (erg s^-1 cm^-2 Hz^-1)

    if mode==-1 and fname=='':
        print('Error: no filter file provided.')
        return 0
    
    # Fundamental constants
    G = 6.67259e-8
    c = 2.9979e10
    h = 6.626075e-27
    kB = 1.380758e-16
    sig = 5.67e-5
    u = 1.66054e-24
    Na = 6.02214e23
    pi = 3.14159265

    # Model Parameters
    Rp = noise_params[0]
    Rs = noise_params[1]
    dist = noise_params[2]
    RA = noise_params[3]
    Dec = noise_params[4]
    Teff = noise_params[5]
    duration = noise_params[6]    # exposure time in hours
    ntran = 1         # number of transits if applicable

    # compute solid angles of star and planet
    theta_star = 6.96e10*Rs/(3.086e18*dist)
    theta_planet = 6.371e8*Rp/(3.086e18*dist)
    omega_star = pi*theta_star*theta_star
    omega_planet = pi*theta_planet*theta_planet

    length2 = len(mod_wave)
    planet_flux = np.zeros(length2)
    bnu = np.zeros(length2)
    star_flux = np.zeros(length2)

    j=length2-1
    k=0

    uwave = mod_wave
    # mod_wave = 10000./mod_wave

    # planet_flux = mod_flux/np.pi * omega_planet
    for i in range(0,len(mod_flux)):
        planet_flux[i] = (float)(mod_flux[i])/pi*omega_planet  # multiply by solid angle

    print("Flux\n")
    print(mod_flux[-10:])
    print(omega_planet)
    print(planet_flux)
    received_flux = planet_flux                  # actual flux density (erg/s/cm^2/Hz) at the telescope
    planet_flux = planet_flux/c/uwave/uwave      # convert cm^-1 to Hz^-1
    planet_flux = planet_flux/(h*c*uwave)        # divide by photon energy
    planet_flux = planet_flux*2.5e5              # multiply by collecting area
    print(planet_flux)
    # this should be received flux in \gamma s^-1 aperture^-1 Hz^-1

    if(starspec==''):
        bnu = 2.0*h*c*pow(uwave,3)/(np.exp(1.44*uwave/Teff)-1.0)
    else:
        fspec = open(starspec,'r')
        lines = fspec.readlines()
        flen = len(lines)
        swave = np.zeros(flen)
        sflux = np.zeros(flen)
        for i in range(0,flen):
            swave[i] = (float)(lines[i].split()[0])
            sflux[i] = (float)(lines[i].split()[1])
            
        sf = interpolate.interp1d(swave,sflux,kind='linear')
        bnu = sf(mod_wave)
        
    star_flux = bnu*omega_star
    star_flux = star_flux/(h*c*uwave)
    star_flux = star_flux*2.5e5

    fil = Filter.Filter()
    
    if mode<=8:
        xthruput = fil.NIRSpec(mode)
        oname = 'NIRSpec4'
    if mode==9:
        xthruput = fil.MIRILRS(True)  # slit spectroscopy
        oname = 'MIRISlit'
    if mode==10:
        xthruput = fil.MIRILRS(False) # slitless spectroscopy
        oname = 'MIRISlitless'
    if mode==11:
        xthruput = fil.NIRCamF322W2()
        oname = 'NIRCam'
    if mode==12:
        xthruput = fil.NIRISS()
        oname = 'NIRISS'
    if mode==13:
        xthruput = fil.MIRIMRS('A')
        oname = 'MIRIA'
    if mode==14:
        xthruput = fil.MIRIMRS('B')
        oname = 'MIRIB'
    if mode==15:
        xthruput = fil.MIRIMRS('C')
        oname = 'MIRIC'
        
    # begin observe.pro -- this does not change the units -- \gamma s^-1 aperture^-1 Hz^-1
    num = int(fil.resolve)
    num6 = int(6.0*fil.resolve)
    sigma = num/2.35
    kernel = np.zeros(num6)
    for i in range(0,num6):
        kernel[i] = np.exp(-0.5*(i-num6/2.)*(i-num6/2.)/sigma/sigma)
    kernel = kernel/np.sum(kernel)
    flux_planet = np.convolve(planet_flux,kernel,mode='same')
    flux_density = np.convolve(received_flux,kernel,mode='same')

    flux_star = np.convolve(star_flux,kernel,mode='same')

    if mode==-1: fil.calwave = mod_wave
    #fil.calwave = 10000./fil.calwave

    if mode!=-1:
        f = interpolate.interp1d(mod_wave,flux_star,kind='linear')
        flux_star = f(fil.calwave)
        f2 = interpolate.interp1d(mod_wave,flux_planet,kind='linear')
        flux_planet = f2(fil.calwave)
        f3 = interpolate.interp1d(mod_wave,flux_density,kind='linear')
        flux_density = f3(fil.calwave)
    
    # begin background.pro
    dum = len(fil.calwave)
    num = dum
    background = np.zeros(num)
    zodi = np.zeros(num)

    # implement standard method EULER to convert RA and DEC to ecliptic coordinates                  
    equinox = '(J2000)'
    psi = 0.0
    stheta = 0.39777715593
    ctheta = 0.91748206207
    phi = 0.0

    ao = RA/57.296 - phi
    bo = Dec/57.296
    sb = np.sin(bo)
    cb = np.cos(bo)
    cbsa = cb*np.sin(ao)
    bo = -stheta*cbsa + ctheta*sb
    elat = np.arcsin(np.modf(bo)[0])

    ao = np.arctan((ctheta*cbsa + stheta*sb) / cb*np.cos(ao)) # make sure output is in the right quadrant
    elong = np.fmod(ao+psi+4*pi,2*pi)
    # end EULER

    elat = np.fabs(elat)

    # localzodi constants
    tau_lz90 = 4.0e-8
    A_lz = 0.6
    alpha = -0.4
    T_lz = 265
    elatrad = elat/57.296

    for i in range(0,num):
        # implement localzodi
        l_cm = fil.calwave[i]*1.e-4
        B_nu = (2.0*h*c/pow(l_cm,3)) / (np.exp(h*c/kB/T_lz/l_cm)-1)
        B_lz = (tau_lz90 * B_nu) / np.sqrt( np.sin(elatrad)*np.sin(elatrad) + pow(A_lz*pow(fil.calwave[i]/11.0,alpha)*np.cos(elatrad),2))
        zodi[i] = B_lz*1.e17
        # end localzodi
        
        # implement JWST_back
        freq = 3.0e14/fil.calwave[i]
        hnk = 1.44e4/fil.calwave[i]
        const = 2.0*h*pow(freq,3)/9.0e20
        T1 = 38.3
        T2 = 44.5
        T3 = 92.0
        T4 = 155.0
        B1 = const/(np.exp(hnk/T1)-1.0)
        B2 = const/(np.exp(hnk/T2)-1.0)
        B3 = const/(np.exp(hnk/T3)-1.0)
        B4 = const/(np.exp(hnk/T4)-1.0)
        S = 0.48*B1 + 0.1*B2 + 3.0e-5*B3 + 9.9e-7*B4
        background[i] = 1.0e17*S
        # end JWST_back

    if(fil.slit<0.0):
        xback = np.sum(background)
        for i in range(0,len(background)): background[i] = xback
        zback = np.sum(zodi)
        for i in range(0,len(zodi)): zodi[i] = zback

    background = 1.0e-17*background
    zodi = 1.0e-17*zodi
    background = background/(h*c*(1.0e4/fil.calwave))
    zodi = zodi/(h*c/(fil.calwave/1.e4))
    # end background.pro

    if(fil.slit > 0.): omeg_back = (fil.slit*fil.scale/206265.)*(1.22*fil.calwave/6.0e6)
    if(fil.slit < 0.): omeg_back = (fil.scale/206265.)*(1.22*fil.calwave/6.0e6)

    if(np.min(fil.calwave) < 0.9): omeg_back = (fil.scale/206265.)*(22.0*fil.scale/206265.)

    background = background * omeg_back * 2.5e5
    zodi = zodi * omeg_back * 2.5e5

    freq1 = 1.e4/fil.calwave                      # cm^-1
    bandwidth = np.fabs(np.gradient(freq1)) * c   # s^-1

    if(np.max(fil.calwave) > 15.0):
        bandwidth[fil.npix1-1] = bandwidth[fil.npix1-2]
        bandwidth[fil.npix1] = bandwidth[fil.npix1+1]
        bandwidth[fil.npix1+fil.npix2-1] = bandwidth[fil.npix1+fil.npix2-2]
        bandwidth[fil.npix1+fil.npix2] = bandwidth[fil.npix1+fil.npix2+1]
        bandwidth[fil.npix1+fil.npix2+fil.npix3-1] = bandwidth[fil.npix1+fil.npix2+fil.npix3-2]
        bandwidth[fil.npix1+fil.npix2+fil.npix3] = bandwidth[fil.npix1+fil.npix2+fil.npix3+1]

    if mode==-1:   # Photometric filter
        fFilter = open(fname,'r')
        flines = fFilter.readlines()
        flen = len(flines)
        fwave = np.zeros(flen)
        fthru = np.zeros(flen)
        for i in range(0,flen):
            fwave[i] = 1.e8/float(flines[i].split()[0])
            fthru[i] = float(flines[i].split()[1])

        wmin = [i for i in range(0,len(mod_wave)) if (mod_wave[i]<fwave[0])]
        wmax = [i for i in range(0,len(mod_wave)) if (mod_wave[i]>fwave[-1])]
        shortwave = mod_wave[wmin[0]:wmax[-1]]
        
        ffil = interpolate.interp1d(fwave,fthru,kind='linear')
        # bandwidth = *width of filter*
        
        # Resolved/directly-imaged
        if obstype==0:
            flux_signal = flux_planet  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi
        # Secondary eclipse
        if obstype==1:
            flux_signal = flux_planet  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi + flux_star
        # Transit
        if obstype==2:
            flux_signal = flux_star  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi
            
        electron_signal = np.mean(flux_signal[wmin[0]:wmax[-1]] * ffil(shortwave))  # e- s^-1 Hz^-1
        electron_noise = np.mean(flux_noise[wmin[0]:wmax[-1]] * ffil(shortwave))

        electron_signal = electron_signal + np.mean(background) # thermal background from telescope, behind the filter
        electron_noise = electron_noise + np.mean(background)

        bandwidth = (fwave[0] - fwave[-1])*2.998e10
        electron_signal = bandwidth * electron_signal  # e- s^-1 bin^-1
        electron_noise = bandwidth * electron_noise
        flux_density = bandwidth * np.mean(flux_density)

        fil.effic = 1.

        signal_electrons = 3600.*ntran*duration*fil.effic*electron_signal  # e- exposure^-1 bin^-1
        noise_electrons = np.sqrt(3600.*ntran*duration*fil.effic*2.0*(electron_signal+electron_noise))
        
        real_signal = 3600.*ntran*duration*fil.effic*electron_signal
        real_noise = noise_electrons/real_signal
        fnoise = noise_electrons/signal_electrons

        return (fwave[0]+fwave[-1])/2.,flux_density,fnoise
    
    else:
        # Resolved/directly-imaged
        if obstype==0:
            flux_signal = flux_planet  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi
        # Secondary eclipse
        if obstype==1:
            flux_signal = flux_planet  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi + flux_star
        # Transit
        if obstype==2:
            flux_signal = flux_star  # \gamma s^-1 aperture^-1 Hz^-1
            flux_noise = zodi
        print(flux_signal)
        print(bandwidth)
        
        electron_signal = flux_signal * xthruput  # e- s^-1 Hz^-1
        electron_noise = flux_noise * xthruput
    
        electron_signal = electron_signal + background # thermal background from telescope, behind the filter
        electron_noise = electron_noise + background

        electron_signal = bandwidth * electron_signal  # e- s^-1 bin^-1
        electron_noise = bandwidth * electron_noise
        #flux_density = bandwidth * flux_density       # Was meant to return flux per bin, replaced by flux per unit wavelength.
        # end observe.pro

        signal_electrons = 3600.*ntran*duration*fil.effic*electron_signal  # e- exposure^-1 bin^-1
        #print(signal_electrons)
        
        noise_electrons = np.sqrt(3600.*ntran*duration*fil.effic*2.0*(electron_signal+electron_noise))
        #print(noise_electrons)

        print("Electrons\n")
        for i in range(0,len(signal_electrons)):
            if i%100==0: print(signal_electrons[i]/noise_electrons[i])
        
        real_signal = 3600.*ntran*duration*fil.effic*electron_signal
        real_noise = noise_electrons/real_signal
        fnoise = noise_electrons/signal_electrons

        return fil.calwave,flux_density,fnoise
