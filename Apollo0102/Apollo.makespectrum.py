import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from scipy.interpolate import interp1d
from distutils.util import strtobool
from src import wrapPlanet_auto
from src import wrapPlanet_layer
from src import AddNoise

# An attempt at adding a multi-nested sampling option.
#multinest = True
'''
try:
    import pymultinest
except:
    print 'PyMultiNest not available.'
    multinest = False
    pass
'''
'''
Required files:

Apollo.py
setup.py
AddNoise.py
Filter.py
FilterList.dat

Planet_layer.cpp
Planet_layer.h
Planet_auto.cpp
Planet_auto.h
Atmosphere.cpp
Atmosphere.h
constants.cpp
constants.h
specincld.h

wrapPlanet_layer.pyx
wrapPlanet_layer.pxd
wrapPlanet_auto.pyx
wrapPlanet_auto.pxd

Opacities directory
Filters directory

Various example files included
'''

def GetBins(specwave,obshi,obslo):
    binshi = np.zeros(len(obshi))
    binslo = np.zeros(len(obshi))
    if(obshi[0]>specwave[0] or obslo[-1]<specwave[-1]):
        print "Wavelength out of range."
        return 0.
    for i in range(0,len(obshi)):
        for j in range(0,len(specwave)):
            if(obshi[i]>specwave[j]):
                binshi[i] = float(j) + (obshi[i]-specwave[j-1])/(specwave[j]-specwave[j-1])
                break
        for j in range(0,len(specwave)):
            if(obslo[i]>=specwave[j]):
                binslo[i] = float(j) + (obslo[i]-specwave[j-1])/(specwave[j]-specwave[j-1])
                break
    return [binshi,binslo]

def GetBinnedSpec(mode,mod_wave,mod_flux,binshi,binslo):
    # Note: the bins are indices, not wavenumbers.
    answer = np.zeros(len(binshi))
    frach = (1.-np.modf(binshi)[0])-0.5
    fracl = np.modf(binslo)[0]-0.5
    width = binslo-binshi
    
    for i in range(0,len(binshi)):
        if binshi[i]>0 and binslo[i]>0:

            answer[i] = np.sum(mod_flux[(int)(np.ceil(binshi[i])):(int)(np.ceil(binslo[i]))])
            answer[i] = answer[i] + frach[i]*mod_flux[(int)(np.floor(binshi[i]))]
            if (int)(np.ceil(binslo[i]))==len(mod_wave):
                answer[i] = answer[i] + fracl[i]*mod_flux[(int)(np.floor(binslo[i]))]
            else:
                answer[i] = answer[i] + fracl[i]*mod_flux[(int)(np.ceil(binslo[i]))]
            answer[i] = answer[i]/width[i]

    return answer

def GetScaOpac(gases,abunds):
    h2 = 1. - np.sum(10**abunds)
    scaopac = 0.672e-27 * h2 # Was 0.662
    mmw = 2.28 * h2 # includes helium
    if 'h2o' in gases:
        scaopac = scaopac + 2.454e-27 * 10**(abunds[gases.index('h2o')-1]) # Was 2.45
        mmw = mmw + 18. * 10**(abunds[gases.index('h2o')-1])
    if 'ch4' in gases:
        scaopac = scaopac + 6.50e-27 * 10**(abunds[gases.index('ch4')-1]) # Was 7.49
        mmw = mmw + 16. * 10**(abunds[gases.index('ch4')-1])
    if 'co' in gases:
        scaopac = scaopac + 4.14e-27 * 10**(abunds[gases.index('co')-1]) # Was 4.14
        mmw = mmw + 28. * 10**(abunds[gases.index('co')-1])
    if 'co2' in gases:
        scaopac = scaopac + 6.82e-27 * 10**(abunds[gases.index('co2')-1]) # Was 8.24
        mmw = mmw + 44. * 10**(abunds[gases.index('co2')-1])
    if 'nh3' in gases:
        scaopac = scaopac + 4.80e-27 * 10**(abunds[gases.index('nh3')-1]) # Was 5.37
        mmw = mmw + 17. * 10**(abunds[gases.index('nh3')-1])
    if 'Burrows_alk' in gases:
        scaopac = scaopac + 718.9e-27 * 10**(abunds[gases.index('Burrows_alk')-1])
        mmw = mmw + 24.1 * 10**(abunds[gases.index('Burrows_alk')-1])
    if 'Lupu_alk' in gases:
        scaopac = scaopac + 718.9e-27 * 10**(abunds[gases.index('Lupu_alk')-1])
        mmw = mmw + 24.1 * 10**(abunds[gases.index('Lupu_alk')-1])
    # Note: these four are based on computed values, which can vary by up to 2 orders of magnitude.
    if 'crh' in gases:
        scaopac = scaopac + 84.0e-27 * 10**(abunds[gases.index('crh')-1]) # Assuming alpha = 8.8e-24 cm^3
        mmw = mmw + 53.0 * 10**(abunds[gases.index('crh')-1])
    if 'feh' in gases:
        scaopac = scaopac + 84.0e-27 * 10**(abunds[gases.index('tio')-1]) # Not available; set equal to CrH.
        mmw = mmw + 56.8 * 10**(abunds[gases.index('feh')-1])
    if 'tio' in gases:
        scaopac = scaopac + 183.3e-27 * 10**(abunds[gases.index('tio')-1]) # Assuming alpha = 13e-24 cm^3
        mmw = mmw + 63.9 * 10**(abunds[gases.index('tio')-1])
    if 'vo' in gases:
        scaopac = scaopac + 131.3e-27 * 10**(abunds[gases.index('vo')-1]) # Assuming alpha = 11e-24 cm^3
        mmw = mmw + 66.9 * 10**(abunds[gases.index('vo')-1])
    return mmw,scaopac

def GetMollist(gases):
    mollist = np.zeros(len(gases))
    gaslist = ["h2","h2only","he","h-","h2o","ch4","co","co2","nh3","h2s","Burrows_alk","Lupu_alk","crh","feh","tio","vo","hcn","n2","ph3"]
    for i in range(0,len(gases)):
        if gases[i] in gaslist: mollist[i] = gaslist.index(gases[i])
        else: mollist[i] = 0
    return mollist

# Start of main program

# Read in planet parameters

settings = 'examples/example.resolved.dat' # Bundled example file
if len(sys.argv)>1: settings = sys.argv[1]

try: fparams = open(settings,'r')
except:
    print "\nError: input file not found.\n"
    sys.exit()

lines1 = fparams.readlines()
pllen = -1
for i in range(0,len(lines1)):
    if len(lines1[i].split())>=6: pllen = pllen+1

fparams.close()
fparams = open(settings,'r')

# Default Settings
mode = 0             # Emission spectrum
modestr = 'Resolved' # Needed for output file name
name = 'example'     # Bundled example file
tstar = 5770.        # Solar temperature
rstar = 1.0          # Solar radius
sma = 1.0            # Semi-Major Axis
dist = 10.0          # Standardized distance, 10 pc
RA = 0.0             # Right ascension
dec = 0.0            # Declination
wavei = 10000./0.60  # Full NIR wavelength range
wavef = 10000./4.99
degrade = 1          # Always 1 for making the spectrum
wrange = 0
minP = 0.0           # Pressure range to integrate over in cgs
maxP = 9.0
vres = 71            # number of layers for radiative transfer
streams = 1          # Use 1-stream by default
hazetype = 0         # No Clouds
hazestr = 'None'     # No Clouds
cloudmod = 0         # No Clouds
natm = 0             # Placeholder in case T-P profile is omitted
verbatim = False     # Interpolate the T-P profile
gray = False         # Used to create a gray atmosphere for testing
tgray = 1500         # Temperature of gray atmosphere
norad = False        # Flags if no radius variable is in the input file
obsfile = ''         # Optional file of observations
opacdir = '../Opacities' # Default opacities directory
starspec = ''        # Optional stellar spectrum file
short = False        # Switch to create short output file names
hires = 'nir'        # Default set of opacity tables to compute the spectra
lores = 'lores'      # Default low-resolution tables to compute Teff
outmode = ''

hazelist = ['None','H2SO4','Polyacetylene','Tholin','Corundum','Enstatite','Forsterite','Iron','KCl','Na2S','NH3Ice','Soot','H2OCirrus','H2OIce','ZnS']

# Read in settings
nlines = 0
while(True):
    last_pos = fparams.tell()
    line = fparams.readline().split()
    if len(line) >= 6:     # ends the loop when the parameters start
        if line[0]=='Parameter':
            break
        else:
            fparams.seek(last_pos)
            break
        
    nlines = nlines+1
    if nlines>100: break # prevents getting stuck in an infinite loop

    if line[0]=='Mode':
        if len(line)>1: modestr = line[1]
        if modestr=='Resolved': mode = 0
        if modestr=='Eclipse': mode = 1
        if modestr=='Transit': mode = 2
    elif line[0]=='Object':
        if len(line)>1: name = line[1]
    elif line[0]=='Star':
        if len(line)>1: tstar = (float)(line[1])
        if len(line)>2: rstar = (float)(line[2])
        if len(line)>3: sma = (float)(line[3])
    elif line[0]=='Location':
        if len(line)>1: dist = (float)(line[1])
        if len(line)>2: RA = (float)(line[2])
        if len(line)>3: dec = (float)(line[3])
    elif line[0]=='Data':
        if len(line)>1: obsfile = line[1]
    elif line[0]=='Star_Spec':
        if len(line)>1: starspec = line[1]
    elif line[0]=='Spectrum':
        if len(line)>1:
            if line[1]=='NIR': wrange = 0
            if line[1]=='MIR': wrange = 1
    elif line[0]=='Pressure':
        if len(line)>1: minP = (float)(line[1]) + 6.0 # Convert from bars to cgs
        if len(line)>2: maxP = (float)(line[2]) + 6.0
    elif line[0]=='Vres':
        if len(line)>1: vres = (int)(line[1])
    elif line[0]=='Streams':
        if len(line)>1: streams = (int)(line[1])
    elif line[0]=='Gray':
        if len(line)>1: gray = strtobool(line[1])
        if len(line)>2: tgray = line[2]
    elif line[0]=='Output':
        if len(line)>1 and line[1]=='Short': short = True
        if len(line)>2 and line[2]=='Short': short = True
        if len(line)>3 and line[3]=='Short': short = True
    elif line[0]=='Opacities':
        if len(line)>1: opacdir = line[1]
    elif line[0]=='Tables':
        if len(line)>1: hires = line[1]
        if len(line)>2: lores = line[2]
    elif line[0]=='Output_Mode':
        if len(line)>1: outmode = line[1]
        
# End read in settings

# Set output file: Object name, type of observation, and # of parameters
if short:
    outfile = '/' + name + '.'
else:
    outfile = '/' + name + '.' + modestr + '.' + str(pllen) + 'params.'
    
plparams = np.zeros(pllen)     # parameter list
mu       = np.zeros(pllen)     # Gaussian means
sigma    = np.zeros(pllen)     # Standard errors
bounds   = np.zeros((pllen,2)) # Bounds

# Read in parameters
lines = fparams.readlines()

i=0
state = -1
basic = []
gases = []
clouds = []
atmtype = 'Layers' # Default layered atmosphere
smooth = False     # Default (needed to correctly count parameters)
igamma = -1        # Index of smoothing parameter if included

b1 = -1
bnum = 0
g1 = -1
gnum = 0
a1 = -1
anum = 0
c1 = -1
cnum = 0
e1 = -1
enum = 0

for j in range(0,len(lines)):
    if str(lines[j]) == 'Basic\n':
        state = 0
        b1 = i
    elif str(lines[j]) == 'Gases\n':
        state = 1
        g1 = i
        if len(lines[j].split())>1:
            gases.append(lines[j].split()[1])
        else:
            gases.append('h2only')  # Filler gas is H2-only, may change
    elif str(lines[j].split()[0]) == 'Atm':
        state = 2
        a1 = i
        atmtype = lines[j].split()[1]
        if len(lines[j].split())>2:
            if str(lines[j].split()[2])=='Verbatim': verbatim = True
    elif str(lines[j].split()[0]) == 'Haze' or str(lines[j].split()[0]) == 'Clouds':
        state = 3
        c1 = i
        cloudmod = int(lines[j].split()[1])
        if len(lines[j].split())>=3:
            hazestr = str(lines[j].split()[2])
            if hazestr in hazelist: hazetype = hazelist.index(hazestr)
    elif str(lines[j]) == 'End\n':
        state = 4
        e1 = i
    elif len(lines[j].split())<6:
        print 'Error: missing parameter values.'
        sys.exit()
    else:
        if state==0:
            basic.append(lines[j].split()[0])
            bnum = bnum+1
        if state==1:
            gases.append(lines[j].split()[0])
            gnum = gnum+1
        if state==2:
            anum = anum+1
        if state==3:
            clouds.append(lines[j].split()[0])
            cnum = cnum+1
        if state==4:
            enum = enum+1
        plparams[i] = (float)(lines[j].split()[1])
        mu[i]       = (float)(lines[j].split()[2])
        sigma[i]    = (float)(lines[j].split()[3])
        bounds[i,0] = (float)(lines[j].split()[4])
        bounds[i,1] = (float)(lines[j].split()[5])
        if lines[j].split()[0]=='gamma':
            smooth = True
            igamma = j
        i = i+1

if gray: gases = []

b2 = b1+bnum
g2 = g1+gnum
a2 = a1+anum
c2 = c1+cnum
e2 = e1+enum
if smooth:
    a2 = a2-1
    anum = anum-1
ilen = int(10 + c2-c1)
ngas = g2-g1+1

if 'RtoD2U' in basic:
    pos = basic.index('RtoD2U')
    plparams[b1+pos] = 10**plparams[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2
    
# Meant to make the temperatures uniform in log space
'''
if atmtype == 'Layers':
    print 'Layers\n'
    plparams[a1:a2] = np.log10(plparams[a1:a2])
    bounds[a1:a2,:] = np.log10(bounds[a1:a2,:])
'''

# End of read in parameters


'''
The C++ function calls are as follows:

planet = wrapPlanet_layer.PyPlanet()
Calls the wrapper class for the constructor

planet.MakePlanet(mode,specwave,hazetype,mollist)
mode: 0 for emission, 1 for transit
specwave: array of wavelengths at which to compute the forward model
hazetype: 0 for no haze, 1 for H2SO4, 2 for polyacetylene, 3 for tholin
mollist: array of integer codes for molecules in the model

planet.set_Params(params1,abund,tplong)
params1: array of basic system properties
abund: array of molecular abundances
tplong: detailed temperature profile dervied from the T-P parameterization

planet.get_Spectrum(streams)
streams: 1 for 1-stream radiative transfer, 2 for 2-stream
'''

# Set wavelength spectrum, keeping full range for spectrum generation

# Compute hires spectrum wavelengths
opacfile = opacdir + '/h2o.' + hires + '.dat'
fopac = open(opacfile,'r')
opacshape = fopac.readline().split()
fopac.close()
nwave = (int)(opacshape[6])
lmin = (float)(opacshape[7])
resolv = (float)(opacshape[9])

speclen = (int)(nwave/degrade)
specwave = np.zeros(speclen)
for i in range(0,speclen):
    specwave[i] = 10000./(lmin*np.exp(i*degrade/resolv))
    
# Compute lores spectrum wavelengths
opacfile = opacdir + '/h2o.' + lores + '.dat'
fopac = open(opacfile,'r')
opacshape = fopac.readline().split()
fopac.close()
nwavelo = (int)(opacshape[6])
lminlo = (float)(opacshape[7])
resolvlo = (float)(opacshape[9])

specwavelo = np.zeros(nwavelo)
for i in range(0,nwavelo):
    specwavelo[i] = 10000./(lminlo*np.exp(i/resolvlo))
    
mmw,rxsec = GetScaOpac(gases,plparams[g1:g2])
mollist = GetMollist(gases)

natm = a2-a1
profin = np.zeros(natm)
for i in range(0,natm): profin[i] = maxP + (minP-maxP)*i/(natm-1)

if atmtype == 'Parametric' and natm != 5:
    print 'Error: wrong parameterization of T-P profile.'
    sys.exit()
    
# Create Planet and read in opacity tables
if atmtype == 'Layers':
    planet = wrapPlanet_layer.PyPlanet()
if atmtype == 'Parametric':
    planet = wrapPlanet_auto.PyPlanet()

print 'Haze type:',hazestr
print 'Cloud model:',cloudmod
mode = int(mode)
cloudmod = int(cloudmod)
hazetype = int(hazetype)
switches = [mode,cloudmod,hazetype]

# Switch for testing without calling the C++ code.
#sys.exit()

planet.MakePlanet(switches,specwave,specwavelo,mollist,opacdir,hires,lores)
print 'Setup complete.'
# End of setup

# Based on likelihood function
params = plparams

if len(gases)==0:
    abund = np.zeros(1)
    abund[0] = 1.
    mmw = 2.28
    rxsec = 0.
else:
    abund = np.zeros(len(gases))
    for i in range(1,len(gases)): abund[i] = params[g1+i-1]
    abund[0] = 1.-np.sum(10**params[g1:g2])
    mmw,rxsec = GetScaOpac(gases,abund[1:])

theta_planet = 0.
params1 = np.zeros(ilen)

# Radius handling
if 'Rad' in basic:
    pos = basic.index('Rad')
    params1[0] = params[b1+pos]
    theta_planet = params[b1+pos]*6.371e8/dist/3.086e18
elif 'RtoD' in basic:
    pos = basic.index('RtoD')
    params1[0] = 10**params[b1+pos]*dist*4.838e9 # convert R/D to Earth radii
    theta_planet = 10**params[b1+pos]
elif 'RtoD2U' in basic:
    pos = basic.index('RtoD2U')
    params1[0] = np.sqrt(params[b1+pos])
    theta_planet = np.sqrt(params[b1+pos])*6.371e8/dist/3.086e18
else:
    norad = True
    params1[0] = 11.2
    theta_planet = 11.2*6.371e8/dist/3.086e18
    # Default radius = Jupiter

# Gravity handling
if 'Log(g)' in basic:
    pos = basic.index('Log(g)')
    params1[1] = params[b1+pos]
else:
    params1[1] = 4.1
    # Default gravity = Jupiter

# Cloud deck handling
if 'Cloud_Base' in clouds:
    pos = clouds.index('Cloud_Base')
    params1[2] = params[c1+pos]+6.
elif 'P_cl' in clouds:
    pos = clouds.index('P_cl')
    params1[2] = params[c1+pos]+6.
else:
    params1[2] = 2.5+6.
    # Default cloudless
if params1[2] < minP: params1[2] = minP+0.01
# Ensures the cloud deck is inside the model bounds.

params1[3] = tstar
params1[4] = rstar
params1[5] = mmw
params1[6] = rxsec
params1[7] = minP
params1[8] = maxP
params1[9] = sma

if hazetype!=0:
    if cloudmod==2 or cloudmod==3:
        for i in range(0,4): params1[i+10] = params[c1+i]
        params1[11] = params1[11]
        params1[12] = params1[12] + 6.
    if cloudmod==2: params1[13] = params1[12] + params1[13]
    
tpprof = np.zeros(natm)

# Gray atmosphere approximation
if gray:
    tplong = np.zeros(vres)
    for i in range(0,vres): tplong[i] = tgray
    planet.set_Params(params1,abund,tplong)
    specflux = planet.get_Spectrum(streams)

# Build atmosphere temperature profile and compute spectrum
else:
    for i in range(0,len(tpprof)): tpprof[i] = params[i+a1]
    if atmtype=='Parametric': tpprof[1] = 10**tpprof[1]

    if atmtype == 'Layers':
        if verbatim: tplong = tpprof
        elif natm==0:
            tplong = np.zeros(vres)
            for i in range(0,vres): tplong[i] = 1500.
        elif natm==1:
            tplong = np.zeros(vres)
            for i in range(0,vres):
                if tpprof[0]<75.: tplong[i]=75.
                elif tpprof[0]>4000.: tplong[i]=4000.
                else: tplong[i]=tpprof[0]                     
        else:
            # compute cubic spline T-P profile
            tplong = np.zeros(vres)
            if natm<=4: f = interp1d(profin,tpprof,kind='linear')
            else: f = interp1d(profin,tpprof,kind='cubic')
            for i in range(0,vres):
                tplong[i] = f(maxP + (minP-maxP)*i/(float)(vres-1))
                if(tplong[i]<75.): tplong[i]=75.
                if(tplong[i]>4000.): tplong[i]=4000.

        # Compute spectrum
        planet.set_Params(params1,abund,tplong)        
        specflux = planet.get_Spectrum(streams)

    if atmtype == 'Parametric':
        # Compute spectrum
        planet.set_Params(params1,abund,tpprof)
        specflux = planet.get_Spectrum(streams)
        
# Set noise parameters
noise_params = np.zeros(7)
noise_params[0] = params1[0]  # radius
noise_params[1] = 1.0         # r_star
noise_params[2] = dist
noise_params[3] = 0.0
noise_params[4] = 0.0
noise_params[5] = 6500.       # t_star
noise_params[6] = 3.0         # exposure time in hours

wavels = 10000./specwave

# multiply by solid angle and collecting area
newtdepth = np.zeros(len(specflux))
if mode<=1:
    for i in range(0,len(specflux)):
        newtdepth[i] = specflux[i] * theta_planet*theta_planet
        #if i==0: print "newtdepth: ",i,specwave[i],newtdepth[i]
        # erg/s/aperture/Hz
        # theta_planet is actually the radius/distance ratio
        # so its square converts flux at the surface to flux at the telescope
if mode==2:
    newtdepth = specflux

# End of likelihood function emulator

# Full resolution output
    
foutname = 'modelspectra' + outfile + 'fullres.dat'
ftest = open(foutname,'w')
for i in range(0,len(specwave)-1):
    ftest.write('{0:8.2f} {1:8.2f} {2:8.5e} 0.0 0.0 {2:8.5e}\n'.format(specwave[i],specwave[i+1],newtdepth[i]))

teff = planet.get_Teff()
print 'Teff = {0:8.2f} K'.format(teff)
    
# Continue if an observation file was named
if len(obsfile)>0:
    
    # Output binned to the observations
    fobs = open(obsfile,'r')

    rlines = fobs.readlines()
    rlen   = len(rlines)
    rcalhi = np.zeros(rlen)
    rcallo = np.zeros(rlen)
    rflux  = np.zeros(rlen)
    rnoise = np.zeros(rlen)

    for i in range(0,rlen):
        rcalhi[i] = (float)(rlines[i].split()[0])
        rcallo[i] = (float)(rlines[i].split()[1])
        rflux[i]  = (float)(rlines[i].split()[2])
        rnoise[i] = (float)(rlines[i].split()[3])

    #rcalmid = 10000./((rcalhi + rcallo)/2.) # Is that the wavelength offset?
    rcalmid = (rcalhi + rcallo)/2.
    rbandwidth = np.fabs(np.gradient(rcalmid)) * 2.99792458e10
    
    fr = interp1d(specwave,newtdepth,kind='linear')  # Cubic?
    bflux = fr(rcalmid)

    foutname = 'modelspectra' + outfile + 'binned.dat'
    fout = open(foutname,'w')
    for i in range(0,rlen):
        fout.write('{0:8.2f} {1:8.2f} {2:8.5e} 0.0 0.0 {2:8.5e}\n'.format(rcalhi[i],rcallo[i],bflux[i]))

    rcalmid = 10000./rcalmid
    specwave = 10000./specwave
    
    xmin = max(min(specwave),min(rcalmid))
    xmax = min(max(specwave),max(rcalmid))
    xmin = xmin - 0.05*(xmax-xmin)
    xmax = xmax + 0.05*(xmax-xmin)
    
    dlist = [i for i in range(0,len(specwave)) if xmin < specwave[i] < xmax]
    rlist = [i for i in range(0,len(rcalmid)) if xmin < rcalmid[i] < xmax]
    
    yref = max(max(newtdepth[dlist]),max(rflux[rlist]))
    
    ymin = -0.20 * yref
    ymax =  1.05 * yref

    # Plot with binned spectrum

    plotname = 'plots' + outfile + 'binned.fit.png'  # Need to put this in the plots folder.

    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))

    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)

    # May need a way to set the labels based on the input files
    
    #ax.plot(rcalmid,rflux,'-',linewidth=1,label='Observations',c='k')
    ax.plot(rcalmid,bflux,'-',linewidth=1,label='Retrieved Spectrum (binned)',c='b')
    ax.errorbar(rcalmid,rflux,rnoise,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')

    residuals = bflux-rflux
    ax.plot(rcalmid,residuals+ymin/2.,'-',linewidth=1,label='Residuals (binned and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')
    
    # Add average error bar symbol
    
    plt.legend(fontsize=12)
    plt.savefig(plotname) # Need to procedurally generate filename

    # Plot with full-res spectrum
    
    plotname = 'plots' + outfile + 'fullres.fit.png'  # Need to put this in the plots folder.
    
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))

    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)

    # May need a way to set the labels based on the input files
    
    #ax.plot(rcalmid,rflux,'-',linewidth=1,label='Observations',c='k')
    ax.plot(specwave,newtdepth,'-',linewidth=1,label='Retrieved Spectrum',c='b')
    ax.errorbar(rcalmid,rflux,rnoise,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')

    residuals = bflux-rflux
    ax.plot(rcalmid,residuals+ymin/2.,'-',linewidth=1,label='Residuals (binned and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')
    
    # Add average error bar symbol
    
    plt.legend(fontsize=12)
    plt.savefig(plotname) # Need to procedurally generate filename
    
# Create binned files for particular JWST modes.
    
if outmode == '':
    sys.exit()
    
elif outmode == 'NIRISS':

    calwave,flux_density,fnoise = AddNoise.addNoise(12,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.

    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    if mode<=1: noise = fnoise*obsdepth
    if mode==2: noise = fnoise
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.

    foutname = 'modelspectra' + outfile + 'NIRISS.dat'
    ftest = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        ftest.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRCam':

    calwave,flux_density,fnoise = AddNoise.addNoise(11,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.

    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRCam.dat'
    f1 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f1.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRSpec1H':

    calwave,flux_density,fnoise = AddNoise.addNoise(0,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.

    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRSpec1H.dat'
    f1 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f1.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRSpec2H':

    calwave,flux_density,fnoise = AddNoise.addNoise(1,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.
    
    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRSpec2H.dat'
    f2 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f2.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRSpec3H':

    calwave,flux_density,fnoise = AddNoise.addNoise(2,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.
    
    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRSpec3H.dat'
    f3 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f3.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRSpec4H':

    calwave,flux_density,fnoise = AddNoise.addNoise(3,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.
    
    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRSpec4H.dat'
    f4 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f4.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'NIRSpec4M':

    calwave,flux_density,fnoise = AddNoise.addNoise(7,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.

    obsdepth = np.zeros(len(calwave))
    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'NIRSpec4M.dat'
    f5 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f5.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

elif outmode == 'Prism':

    calwave,flux_density,fnoise = AddNoise.addNoise(8,mode,specwave,newtdepth,noise_params,starspec)
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2.
        callo[i] = (calwave[i]+calwave[i+1])/2.
    
    obsdepth = GetBinnedSpec(mode,specwave,newtdepth,10000./calhi,10000./callo)

    noise = fnoise*obsdepth
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise
    for i in range(0,len(calwave)):
        if(obs_flux[i]<0.): obs_flux[i] = 0.
   
    foutname = 'modelspectra' + outfile + 'Prism.dat'
    f6 = open(foutname,'w')
    for i in range(0,len(calwave)-1):
        f6.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(calhi[i],callo[i],obsdepth[i],noise[i],obs_flux[i]))

else:
    if os.path.exists(outmode):
        calwave,flux_density,fnoise = AddNoise.addNoise(-1,mode,specwave,newtdepth,noise_params,starspec,outmode)
        print 'Band Center: {0:8.2f} cm^-1'.format(calwave)
        print 'Flux:     {0:12.5e} cgs'.format(flux_density)
        print 'Noise:    {0:12.5e} cgs'.format(fnoise*flux_density)
    else:
        print 'Error: filter file not found.'