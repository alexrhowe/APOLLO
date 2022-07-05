from __future__ import print_function
from inspect import getmembers, isfunction
import os
import sys
import numpy as np
import pandas as pd
import pickle
import scipy.optimize as op
from scipy.interpolate import interp1d
from scipy.stats import invgamma
from distutils.util import strtobool
import emcee
import schwimmbad

import matplotlib
# Used on Discover because the GUI backend does not work there. Comment out to use 'Manual' mode.
#matplotlib.use('pdf')
import matplotlib.pyplot as plt

# Comment out corner if not supported on your platform.
#import corner
from src import wrapPlanet
from src import ApolloFunctions as af
from src import AddNoise
from src.defaults import *
from src.P_points import P_profiles
from src import TP_profiles

# An attempt at adding a multi-nested sampling option.
#multinest = True

REarth = 6.371e8   # R_Earth in cm
parsec = 3.086e18  # parsec in cm
RJup = 11.2        # R_Jupiter in R_Earth

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
ApolloFunctions.py
AddNoise.py
Filter.py
setup.py

Planet.cpp
Planet.h
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

Various example files included
'''

#----------------------------------------------------------------------------------------#
# Read in input file

settings = 'examples/example.resolved.dat'  # Bundled example file
if len(sys.argv)>1: settings = sys.argv[1]

task = 'Spectrum'    # Default: makes a spectrum plot
if len(sys.argv)>2:
    if sys.argv[2]=='Spectrum': task = 'Spectrum'
    elif sys.argv[2]=='Retrieval': task = 'Retrieval'
    elif sys.argv[2]=='Ensemble': task = 'Ensemble'
    else:
        print('Error: specify "Spectrum" or "Retrieval" or "Ensemble".')
        sys.exit()
    
override = False
manual = False
if len(sys.argv)>3 and sys.argv[3]=='Serial': override = True
if len(sys.argv)>3 and sys.argv[3]=='Manual': manual = True

try: fparams = open(settings,'r')
except:
    print("\nError: input file not found.\n")
    sys.exit()

lines1 = fparams.readlines()
pllen = -1
for i in range(0,len(lines1)):
    if len(lines1[i].split())>=6: pllen = pllen+1

fparams.close()
fparams = open(settings,'r')
#----------------------------------------------------------------------------------------#
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

    elif line[0]=='Object':
        if len(line)>1: name = line[1]
    if line[0]=='Mode':
        if len(line)>1: modestr = line[1]
        if modestr=='Resolved': mode = 0
        if modestr=='Eclipse': mode = 1
        if modestr=='Transit': mode = 2
    elif line[0]=='Parallel':
        if len(line)>1: parallel = strtobool(line[1])
    elif line[0]=='Plotting':
        if len(line)>1: plotting = strtobool(line[1])
    elif line[0]=='Data':
        if len(line)>1: datain = line[1]
        if len(line)>2:
            if line[2]=='Polyfit': polyfit = True
    elif line[0]=='Convolve':
        if len(line)>1: dataconv = (int)(line[1])
    elif line[0]=='Binning':
        if len(line)>1: databin = (int)(line[1])
    elif line[0]=='Degrade':
        if len(line)>1: degrade = (int)(line[1])  # Only works with low-res data; mainly used to speed up execution for testing
    elif line[0]=='Prior':
        if len(line)>1: prior = line[1]
    elif line[0]=='N_Walkers':
        if len(line)>1: nwalkers = (int)(line[1])
        if override: nwalkers = 2
    elif line[0]=='N_Steps':
        if len(line)>1: nsteps = (int)(line[1])
    elif line[0]=='Star':
        if len(line)>1: tstar = (float)(line[1])
        if len(line)>2: rstar = (float)(line[2])
        if len(line)>3: sma = (float)(line[3])
    elif line[0]=='Star_Spec':
        if len(line)>1: starspec = line[1]
    elif line[0]=='Location':
        if len(line)>1: dist = (float)(line[1])
        if len(line)>2: RA = (float)(line[2])
        if len(line)>3: dec = (float)(line[3])
    elif line[0]=='Mass_Limits':
        if len(line)>1: minmass = (float)(line[1])
        if len(line)>2: maxmass = (float)(line[2])
    elif line[0]=='Tables':
        if len(line)>1: hires = line[1]
        if len(line)>2: lores = line[2]
    elif line[0]=='Pressure':
        if len(line)>1: minP = (float)(line[1]) + 6.0  # Convert from bars to cgs
        if len(line)>2: maxP = (float)(line[2]) + 6.0
        if maxP <= minP: maxP = minP + 0.01
        if len(line)>3: P_profile = P_profiles[line[3]] + 6.0
        else: P_profile = None
    elif line[0]=='Gray':
        if len(line)>1: gray = strtobool(line[1])
        if len(line)>2: tgray = line[2]
    elif line[0]=='Vres':
        if len(line)>1: vres = (int)(line[1])
    elif line[0]=='Streams':
        if len(line)>1: streams = (int)(line[1])
    elif line[0]=='Output_Mode':
        if len(line)>1: outmode = line[1]
        if len(line)>2: exptime = (float)(line[2])
    elif line[0]=='Output':
        if len(line)>1: outdir = line[1]
        if len(line)>2:
            if line[2]=='Short': short = True
            if line[2]=='Full': printfull = True
        if len(line)>3:
            if line[3]=='Short': short = True
            if line[3]=='Full': printfull = True
    elif line[0]=='Opacities':
        if len(line)>1: opacdir = line[1]
        
# End read in settings

# Output file name: Object name, type of observation, # of parameters, and # of steps.
if short:
    outfile = '/' + name + '.'
else:
    outfile = '/' + name + '.' + modestr + '.' + str(pllen) + 'params' + str(int(nsteps/1000)) + 'k.'

#----------------------------------------------------------------------------------------#
# Read in observations

print('Reading in observations.')
fobs = open(datain,'r')

obslines = fobs.readlines()
obslength = len(obslines)

wavelo = np.zeros(obslength)
wavehi = np.zeros(obslength)
flux = np.zeros(obslength)
errlo = np.zeros(obslength)
errhi = np.zeros(obslength)

for i in range(0,obslength):
    wavelo[i] = obslines[i].split()[0]
    wavehi[i] = obslines[i].split()[1]
    flux[i] = obslines[i].split()[5]
    errlo[i] = obslines[i].split()[3]
    errhi[i] = obslines[i].split()[4]

# End of read in observations

#----------------------------------------------------------------------------------------#
# Process observations for retrieval

# Separate out individual bands
bandindex,bandlo,bandhi,bandflux,banderr = af.FindBands(wavelo,wavehi,flux,errhi)
nband = len(bandhi)

# Convolve the observations to account for effective resolving power or fit at lower resolving power
convflux,converr = af.ConvBands(bandflux,banderr,dataconv)

# Bin the observations to fit a lower sampling resolution
binlo,binhi,binflux,binerr = af.BinBands(bandlo,bandhi,convflux,converr,databin)
binlen = len(binflux)
binmid = np.zeros(len(binlo))
for i in range(0,len(binlo)): binmid[i] = (binlo[i]+binhi[i])/2.

totalflux = 0
for i in range(0,len(binflux)): totalflux = totalflux + binflux[i]*(binhi[i]-binlo[i])*1.e-4

#----------------------------------------------------------------------------------------#
# Read in model parameters

print('Reading in parameters.')

lines = fparams.readlines()

plparams = np.zeros(pllen)     # Parameter list
mu       = np.zeros(pllen)     # Gaussian means
sigma    = np.zeros(pllen)     # Standard errors
bounds   = np.zeros((pllen,2)) # Bounds
guess    = np.zeros(pllen)     # Used for initial conditions

i=0
state = -1
pnames = []
pvars = []
nvars = []
basic = []
gases = []
atm = []
clouds = []
end   = []
atmtype = 'Layers' # Default layered atmosphere
smooth = False     # Default no smoothing
igamma = -1        # Index of gamma if included

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

ensparams = []

for j in range(0,len(lines)):
    
    if str(lines[j]) == 'Basic\n':
        state = 0
        b1 = i
    elif str(lines[j].split()[0]) == 'Gases':
        state = 1
        g1 = i
        if len(lines[j].split())>1:
            gases.append(lines[j].split()[1])
        else:
            gases.append('h2he')  # Filler gas is h2+he, may get more reliable results from h2 only
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
        print('Error: missing parameter values.')
        sys.exit()
            
    else:
        line = lines[j].split()
        pnames.append(line[0])
        if state==0:
            basic.append(line[0])
            bnum = bnum+1
        if state==1:
            gases.append(line[0])
            gnum = gnum+1
        if state==2:
            atm.append(line[0])
            anum = anum+1
        if state==3:
            clouds.append(line[0])
            cnum = cnum+1
        if state==4:
            if line[0][0]=='S':
                snum = int(line[0][1:])
                if snum<0 or snum>=nband: pass
            end.append(line[0])
            enum = enum+1
        plparams[i] = (float)(line[1])
        guess[i]    = (float)(line[1])
        mu[i]       = (float)(line[2])
        sigma[i]    = (float)(line[3])
        bounds[i,0] = (float)(line[4])
        bounds[i,1] = (float)(line[5])

        if sigma[i] > 0.:
            pvars.append(plparams[i])
            nvars.append(i)
        
        if guess[i]==0: guess[i] = 0.1*bounds[i,1]        # Prevents errors from occuring where parameters are zero.

        # Convert pressure units from bars to cgs
        if pnames[i]=='Cloud_Base' or pnames[i]=='P_cl':
            plparams[i] = plparams[i] + 6.
            guess[i] = guess[i] + 6.
            mu[i] = mu[i] + 6.
            bounds[i,0] = bounds[i,0] + 6.
            bounds[i,1] = bounds[i,1] + 6.
            if bounds[i,0] < minP: bounds[i,0] = minP
            if bounds[i,1] > maxP: bounds[i,1] = maxP
            
        if line[0]=='gamma':
            smooth = True
            igamma = j
        if len(line)>6 and line[6]=='Ensemble':
            ensparams.append(i)
        i = i+1

ndim = len(nvars)

if gray: gases = []
elif gases==[]: gases = ['h2he']

if nwalkers==0: nwalkers = ndim*8            # Default number of walkers
if nwalkers<2*ndim: nwalkers = ndim*2 + 2    # Minimum number of walkers
if nwalkers%2==1: nwalkers = nwalkers + 1    # Number of walkers must be even

# Forces a serial execution for command line debugging.
if override:
    parallel = False
    printfull = True
    nsteps = 2
    nwalkers = ndim*2 + 2

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
    
# Meant to make the temperatures uniform in log space, not currently used
'''
if atmtype == 'Layers':
    print 'Layers\n'
    plparams[a1:a2] = np.log10(plparams[a1:a2])
    guess[a1:a2] = np.log10(guess[a1:a2])
    bounds[a1:a2,:] = np.log10(bounds[a1:a2,:])
'''

minDL = 0
maxDL = 0.

# Set statistical parameters
for n in range(0,enum):
    if end[n][0:4]=='logf':
        plparams[e1+n] = np.log(max(errhi**2) * min(errhi**2))/2.
        guess[e1+n]    = plparams[e1+n]
        mu[e1+n]       = plparams[e1+n]
        sigma[e1+n]    = np.log(10.)
        bounds[e1+n,0] = min(errhi**2) + bounds[e1+n,0]
        bounds[e1+n,1] = max(errhi**2) + bounds[e1+n,1]
    elif end[n][0:6]=='deltaL':
        check0 = bounds[e1+n,0]*0.001*1.1
        if minDL > check0: minDL = check0
        check1 = bounds[e1+n,1]*0.001*1.1
        if maxDL < check1: maxDL = check1

# End of read in model parameters
# End of read in input file

#----------------------------------------------------------------------------------------#
# Set the cross section tables if not already set.
# Note that the "default" assumes a particular set of tables.
if 'deltaL' in end:
    pos = end.index('deltaL')
    
wavei = max(wavelo) + minDL
wavef = min(wavehi) + maxDL
if hires=='':
    if wavei < 5.0 and wavef < 5.0:
        if wavei < 0.6: wavei = 0.6
        hires = 'nir'
    elif wavei < 5.0 and wavef > 5.0:
        if wavei < 0.6: wavei = 0.6
        if wavef > 30.0: wavef = 30.0
        hires = 'wide'
    elif wavei > 5.0 and wavef > 5.0:
        if wavef > 30.0: wavef = 30.0
        hires = 'mir'
        
#----------------------------------------------------------------------------------------#
# Set model spectrum wavelength range

# Compute hires spectrum wavelengths
opacfile = opacdir + '/gases/h2o.' + hires + '.dat'
fopac = open(opacfile,'r')
opacshape = fopac.readline().split()
fopac.close()
nwave = (int)(opacshape[6])
lmin = (float)(opacshape[7])
resolv = (float)(opacshape[9])

opaclen = (int)(np.floor(nwave/degrade)) # THIS SEEMED TO NEED A +1 IN CERTAIN CASES.
opacwave = np.zeros(opaclen)
for i in range(0,opaclen):
    opacwave[i] = lmin*np.exp(i*degrade/resolv)
'''
if wavehi[0]>opacwave[0] or wavelo[-1]<opacwave[-1]:
    trim = [i for i in range(0,len(wavehi)) if (wavehi[i]<opacwave[0] and wavelo[i]>opacwave[-1])]
    wavehi = wavehi[trim]
    wavelo = wavelo[trim]
    flux = flux[trim]
    errlo = errlo[trim]
    errhi = errhi[trim]
    obslength = len(trim)
'''
# Compute lores spectrum wavelengths
opacfile = opacdir + '/gases/h2o.' + lores + '.dat'
fopac = open(opacfile,'r')
opacshape = fopac.readline().split()
fopac.close()
nwavelo = (int)(opacshape[6])
lminlo = (float)(opacshape[7])
resolvlo = (float)(opacshape[9])

modwavelo = np.zeros(nwavelo)
for i in range(0,nwavelo):
    modwavelo[i] = lminlo*np.exp(i/resolvlo)
    
# Set up wavelength ranges
'''
imin = np.where(opacwave<np.max(wavehi))[0]-1
imax = np.where(opacwave<np.min(wavelo))[0]+2
istart = np.where(opacwave<wavehi[0])[0]-1
iend = np.where(opacwave<wavelo[-1])[0]+2

if len(imin)==0: imin = [0]
elif imin[0]<0: imin[0] = 0
if len(imax)==0: imax = [len(opacwave)-1]
elif imax[-1]>=len(opacwave): imax[-1] = len(opacwave)-1
if len(istart)==0: istart = [0]
elif istart[0]<0: istart[0]=0
if len(iend)==0: iend = [len(opacwave)-1]
elif iend[-1]>=len(opacwave): iend[-1] = len(opacwave)-1

# Truncated and band-limited spectrum that does not extend beyond the range of the observations
modwave = np.array(opacwave)[(int)(imin[0]):(int)(imax[0])]
lenmod = len(modwave)
'''
# End set up model spectrum wavelength range

#----------------------------------------------------------------------------------------#
# Handle bands and optional polynomial fitting
bindex, modindex, modwave = af.SliceModel(bandlo,bandhi,opacwave,minDL,maxDL)
polyindex = -1
for i in range(1,len(bindex)):
    if bindex[i][0] < bindex[i-1][0]:
        polyindex = i        
if polyindex==-1: polyfit = False

if polyfit:
    normlo = bandlo[0]
    normhi = bandhi[0]
    normflux = bandflux[0]
    normerr = banderr[0]
    for i in range(1,len(bandlo)):
        normlo = np.r_[normlo,bandlo[i]]
        normhi = np.r_[normhi,bandhi[i]]
        normflux = np.r_[normflux,bandflux[i]]
        normerr = np.r_[normerr,banderr[i]]
    normmid = (normlo+normhi)/2.

    slennorm = []
    elennorm = []
    for i in range(polyindex,len(bandindex)):
        slennorm.append(bandindex[i][0])
        elennorm.append(bandindex[i][1])

    masternorm = af.NormSpec(normmid,normflux,slennorm,elennorm)
    fluxspecial = np.concatenate((normerr[0:slennorm[0]],normflux[slennorm[0]:]),axis=0)
    mastererr = af.NormSpec(normmid,fluxspecial,slennorm,elennorm)

else:
    masternorm = binflux
    mastererr = binerr

if not task=='Retrieval':    
    modwave = opacwave
    
# Get indices of the edges of the observation bins in the model spectrum
bins = af.GetBins(modwave,binlo,binhi)
ibinlo = bins[0]
ibinhi = bins[1]

# Needed to calculate the spectrum with the wavelength offset later.
delmodwave = modwave + 0.001
delbins = af.GetBins(delmodwave,binlo,binhi)
delibinlo = delbins[0]-ibinlo
delibinhi = delbins[1]-ibinhi
    
# End of band handling

#----------------------------------------------------------------------------------------#
# Create Planet and read in opacity tables

mmw,rxsec = af.GetScaOpac(gases,plparams[g1:g2])
mollist = af.GetMollist(gases)

natm = a2-a1
if P_profile is not None:
    profin = P_profile
else:
    profin = maxP + (minP-maxP)*np.arange(natm)/(natm-1)

if atmtype == 'Parametric' and natm != 5:
    print('Error: wrong parameterization of T-P profile.')
    sys.exit()

planet = wrapPlanet.PyPlanet()
print('Haze type:',hazestr)
print('Cloud model:',cloudmod)
mode = int(mode)
cloudmod = int(cloudmod)
hazetype = int(hazetype)

atmmod = 0
TP_functions = dict(getmembers(TP_profiles, isfunction))
if atmtype=='Layers' or atmtype in TP_functions: atmmod = 0
if atmtype=='Parametric': atmmod = 1

switches = [mode,cloudmod,hazetype,streams,atmmod]

finalnames = []
for i in nvars: finalnames.append(pnames[i])
guess  = guess[nvars]
mu     = mu[nvars]
sigma  = sigma[nvars]
bounds = bounds[nvars]

#----------------------------------------------------------------------------------------#
# Ensemble setup

if task=='Ensemble':

    esize = 1
    enstable = []
    n=0
    for i in ensparams:
        epmin = bounds[i,0]
        epmax = bounds[i,1]
        estep = sigma[i]
        enum = int((epmax-epmin)/estep)+1
        esize = esize*enum
        enstable.append([])
        for j in range(0,enum):
            epoint = epmin + j*estep
            enstable[n].append(epoint)
        n=n+1
        
    if esize>1000:
        prompt = input('Warning: this ensemble is more than 1,000 entries and may take significant time. Continue? (y/n): ')
        while not (prompt=='y' or prompt=='Y'):
            if prompt=='n' or prompt=='N':
                sys.exit()
            else:
                prompt = input('Error. Please type \'y\' or \'n\': ')
    elif esize>10000:
        print('Error: this ensemble is more than 10,000 entries and may require several hours and several GB of memory. Please use a smaller ensemble.')
        print('If you wish to continue with this ensemble, comment out this warning in the source code.')
        sys.exit()
    print('Ensemble size: {0:d}'.format(esize))

    eplist = np.zeros((esize,len(plparams)))
    for i in range(0,esize):
        eplist[i] = plparams
        eindices = []
        n = i
        for j in range(0,len(enstable)):
            k = n%len(enstable[j])
            n = int(n/len(enstable[j]))
            eindices.append(k)
            
        for j in range(0,len(eindices)):
            eplist[i][ensparams[j]] = enstable[j][eindices[j]]
    
    foutnamee = 'modelspectra' + outfile + 'ensemble.dat'
    fout = open(foutnamee,'w')

    for i in range(0,len(ensparams)):
        fout.write('     {0:s}'.format(pnames[ensparams[i]]))
        for j in range(0,esize):
            fout.write(' {0:f}'.format(eplist[j][ensparams[i]]))
        fout.write('\n')

# Halt before the C++ code for testing purposes.
#sys.exit()

planet.MakePlanet(switches,modwave,modwavelo,mollist,opacdir.encode('utf-8'),hires.encode('utf-8'),lores.encode('utf-8'))
print('Setup complete.')
# End of setup

'''
C++ functions from wrapPlanet_layer and wrapPlanet_auto

MakePlanet(switches,modwave,modwavelo,mollist,opacdir.encode('utf-8'),hires.encode('utf-8'),lores.encode('utf-8'))
switches  = [mode, cloudmod, hazetype, streams], instruct the code which options to use
modwave   = wavelengths over which to compute the model
modwavelo = wavelengths to use for computing the bolometric luminosity and effective temperature
opacdir   = directory where the cross section tables are found
hires     = set of cross section tables to use for the model
lores     = set of cross section tables to use for the effective temperature

set_Params(params1,abund,tpprof)
params1 = array of parameters describing the planet model that aren't included in the other groups
abund   = table of molcular abundances, the filler gas being first
tpprof  = a temperature profile for wrapPlanet_layer or a set of temperature profile parameters for wrapPlanet_auto

get_Spectrum()
get_Teff
'''

#----------------------------------------------------------------------------------------#
# Function to compute forward model

def GetModel(x):
    if len(gases)==0:
        abund = np.zeros(1)
        abund[0] = 1.
        mmw = 2.28
        rxsec = 0.
    else:
        abund = np.zeros(g2-g1+1)
        asum = 0
        for i in range(g1,g2):
            abund[i-g1+1] = x[i]
            asum = asum + 10**x[i]
        abund[0] = 1. - asum
        mmw,rxsec = af.GetScaOpac(gases, abund[1:])
        
    params1 = np.zeros(ilen)

    # Dummy variables in case they cannot be calculated.
    mass = 0.
    ctoo = 0.
    fetoh = 0.
    teff = 0.
    
    # Radius handling
    if 'Rad' in basic:
        pos = basic.index('Rad')
        params1[0] = x[b1+pos]
        radius = REarth*x[b1+pos]
    else:
        global norad
        norad = True
        params1[0] = RJup
        radius = RJup*REarth
        # Default radius = Jupiter

    # Gravity handling
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        params1[1] = x[b1+pos]
        grav = 10**x[b1+pos]
    else:
        params1[1] = 4.1
        grav = 10**4.1

    # Cloud deck handling
    if 'Cloud_Base' in clouds:
        pos = clouds.index('Cloud_Base')
        params1[2] = x[c1+pos]
    elif 'P_cl' in clouds:
        pos = clouds.index('P_cl')
        params1[2] = x[c1+pos]
    elif 'Cloud_Base' in basic:
        pos = basic.index('Cloud_Base')
        params1[2] = x[b1+pos]
    elif 'P_cl' in basic:
        pos = basic.index('P_cl')
        params1[2] = x[b1+pos]
    else:
        params1[2] = 8.5
        # Default cloudless
    if params1[2] < minP: params1[2] = minP+0.01
    # Ensures the cloud deck is inside the model bounds.

    mass = grav*radius*radius/6.67e-8/1.898e30

    # Compute C/O and [Fe/H]
    carbon = 0.
    oxygen = 0.
    metals = 0.
    ccompounds = ['ch4','co','co2','hcn']
    cmult = [1.,1.,1.,1.]
    ocompounds = ['h2o','co','co2','tio','vo']
    omult = [1.,1.,2.,1.,1.]
    zcompounds = ['h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
    zmult = [16.,12.,28.,44.,14., 32.,24.,24.,52.,56., 64.,67.,26.,28.,31.]

    for i in range(0,len(ccompounds)):
        if ccompounds[i] in gases:
            j = gases.index(ccompounds[i])
            carbon = carbon + cmult[i]*(10**x[g1+j-1]) # -1 because of hydrogen
    for i in range(0,len(ocompounds)):
        if ocompounds[i] in gases:
            j = gases.index(ocompounds[i])
            oxygen = oxygen + omult[i]*(10**x[g1+j-1])
    for i in range(0,len(zcompounds)):
        if zcompounds[i] in gases:
            j = gases.index(zcompounds[i])
            metals = metals + zmult[i]*(10**x[g1+j-1])

    ctoo = carbon/oxygen
    fetoh = np.log10(metals/0.0196)

    # For fractional cloud coverage, generate the cloud-free
    # and cloudy spectra, then mix according to the fraction.
    if 'Cloud_Fraction' in clouds:
        pos = clouds.index('Cloud_Fraction')
        cloud_fraction = x[c1+pos]
    else:
        cloud_fraction = 1
    
    params1[3] = tstar
    params1[4] = rstar
    params1[5] = np.sum(mmw)
    params1[6] = np.sum(rxsec)
    params1[7] = minP
    params1[8] = maxP
    params1[9] = sma
    
    if hazetype!=0:
        if cloudmod==2 or cloudmod==3:
            for i in range(0,4): params1[i+10] = x[c1+i]
            params1[11] = params1[11]
            params1[12] = params1[12] + 6.
        if cloudmod==2: params1[13] = params1[12] + params1[13]

    if cloudmod==4:
        for i in range(0,5): params1[i+10] = x[c1+i]
        params1[11] = params1[11] + 6.
        params1[12] = params1[12] + 6.
        
    tpprof = np.zeros(natm)
    # Gray atmosphere approximation
    if gray:
        tplong = np.zeros(vres)
        for i in range(0,vres): tplong[i] = tgray
        planet.set_Params(params1,abund,rxsec,tplong)
        if cloud_fraction == 0:
            specflux = planet.get_ClearSpectrum()
        elif cloud_fraction == 1:
            specflux = planet.get_Spectrum()
        else:
            specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())

    # Build atmosphere temperature profile and compute spectrum
    else:
        for i in range(0,len(tpprof)): tpprof[i] = x[i+a1]
        if atmtype=='Parametric': tpprof[1] = 10**tpprof[1]

        if atmtype in TP_functions:
            tplong = (TP_functions[atmtype])(*tpprof,num_layers_final=vres,P_min=minP-6,P_max=maxP-6)
            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tplong)
            if cloud_fraction == 0:
                specflux = planet.get_ClearSpectrum()
            elif cloud_fraction == 1:
                specflux = planet.get_Spectrum()
            else:
                specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())
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

            # Reverse the T-P profile--required to convert to the current radiative transfer paradigm.                     
            tplong = tplong[::-1]
            np.save(outdir+outfile+"T-P_array_linear", tplong)

            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tplong)
            if cloud_fraction == 0:
                specflux = planet.get_ClearSpectrum()
            elif cloud_fraction == 1:
                specflux = planet.get_Spectrum()
            else:
                specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())
        if atmtype == 'Parametric':
            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tpprof)
            if cloud_fraction == 0:
                specflux = planet.get_ClearSpectrum()
            elif cloud_fraction == 1:
                specflux = planet.get_Spectrum()
            else:
                specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())

    for i in range(0,nband):
        sname = 'S' + str(i)
        if sname in end:
            pos = end.index(sname)
            bandscale = x[e1+pos]
            for j in range(modindex[i][0],modindex[i][1]):
                specflux[j] = specflux[j] * bandscale
            
    teff = planet.get_Teff()
    if task!='Ensemble':
        print('M/Mj: ',mass)
        print('C/O: ',ctoo)
        print('[Fe/H]: ',fetoh)
        print('Teff: ',teff)

    # Plot the results of get_Spectrum() directly for testing purposes.
    '''
    figtest = plt.figure()
    ax = figtest.add_subplot()
    ax.plot(modwave,specflux)
    plt.savefig('plots/test.png')
    sys.exit()
    '''
    
    return specflux, [mass,ctoo,fetoh,teff]

# End of GetModel function

#----------------------------------------------------------------------------------------#
# Likelihood function for "Retrieval" mode.

def lnlike(x,ibinlo,ibinhi):
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]
        
    modflux, derived = GetModel(params)
    
    theta_planet = 0.
    if 'Rad' in basic:
        pos = basic.index('Rad')
        theta_planet = params[b1+pos]*REarth/dist/parsec
    else:
        theta_planet = params[b1+pos]*REarth/dist/parsec
        # Default radius = Jupiter
        
    # Statistical parameters

    deltaL = np.zeros(nband)
    lnf = np.zeros(nband)-100.0
    
    for n in range(0,enum):
        if end[n][0:6]=='deltaL':
            if end[n]=='deltaL': pos = 0
            else: pos = int(end[n][6:])
            deltaL[pos] = params[e1+n]
        if end[n][0:4]=='logf':
            if end[n]=='logf': pos = 0
            else: pos = int(end[n][4:])
            lnf[pos] = params[e1+n]

    # Multiply by solid angle and collecting area
    fincident = np.zeros(len(modflux))
    if mode<=1:
        for i in range(0,len(modflux)):
            fincident[i] = modflux[i] * theta_planet*theta_planet
            #if i==0: print "newtdepth: ",i,specwave[i],fincident[i] # print output for debugging purposes
            # erg/s/aperture/Hz
            # theta_planet is actually the radius/distance ratio
            # so its square converts flux at the surface to flux at the telescope
    if mode==2:
        fincident = modflux

    # This test figure plots the forward model against the observations to the screen
    # and updates with every sample in Serial mode.
    '''
    li1.set_xdata(modwave)
    li1.set_ydata(fincident)
    li2.set_xdata(wavehi)
    li2.set_ydata(masternorm)
    ax.relim()
    ax.autoscale_view(True,True,True)
    figa.canvas.draw()
    plt.pause(0.05)
    '''

    # Outputs a forward model spectrum pre-processing. Part of the testing suite.
    '''
    fff = open('Sample.dat','w')
    for i in range(0,len(modwave)):
        fff.write('{0:e} {1:e}\n'.format(modwave[i],fincident[i]))
    '''
        
    # Adjust for wavelength calibration error
    newibinlo = np.zeros(len(ibinlo))
    newibinhi = np.zeros(len(ibinhi))

    di = 0
    for i in range(0,nband):
        for j in range(di,di+len(bandlo[i])):
            newibinlo[j] = ibinlo[j] + deltaL[i] * delibinlo[j]
            newibinhi[j] = ibinhi[j] + deltaL[i] * delibinhi[j]
        di = di + len(bandlo[i])
        
    # Bin and normalize spectrum
    if norm:
        normspec = af.NormSpec(modwave,fincident,snormtrunc,enormtrunc)
    else:
        normspec = fincident
        
    # Normalize if no radius was given
    if norad:
        normspec = normspec * totalflux/np.sum(normspec)
        
    # normspec is the final forward model spectrum
    binw = (newibinlo[1]-newibinlo[0])*(dataconv/databin)
    convmod = []
    for i in range(0,len(modindex)):
        convmod.append(af.ConvSpec(normspec[modindex[i][0]:modindex[i][1]],binw))

    convmodfl = [item for sublist in convmod for item in sublist]
    
    binmod = af.BinModel(convmodfl,newibinlo,newibinhi)
    
    iw = [i for i in range(0,len(binmod)) if ((i<polyindex or polyindex==-1) and binmod[i]!=0)]
    s2 = np.zeros(len(mastererr))
    di = 0
    for i in range(0,len(lnf)):
        for j in range(di,di+len(bandlo[i])):
            s2[j] = mastererr[j]*mastererr[j] + np.exp(lnf[i])
        di = di + len(bandlo[i])
        
    likelihood = 0
    
    for i in iw:
        likelihood = likelihood + (masternorm[i]-binmod[i])**2/s2[i] + np.log(2.*np.pi*s2[i])
    likelihood = likelihood * -0.5
    
    if(np.isnan(likelihood) or np.isinf(likelihood)):
        print("Error: ")
        print(params)

    # Uncomment to halt execution after the first sample for testing.
    #sys.exit()
        
    return likelihood, derived
    
# End of likelihood function

#----------------------------------------------------------------------------------------#
# Prior probability

def lnprior(x,derived):
    mass = derived[0]
    teff = derived[3]
    
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]

    priors = np.zeros(ndim)
    for i in range(0,ndim):
        if not bounds[i,0] <= x[i] <= bounds[i,1]:
            print('Out of Bound: {0:s} {1} {2} {3}'.format(pnames[nvars[i]],x[i],bounds[i,0],bounds[i,1]))
            return -np.inf
        if smooth and nvars[i]==a2:
            priors[i] = invgamma.pdf(x[i],1,scale=5.e-5) # gamma with inverse gamma function prior, alpha=1, beta=5.e-5
        else:
            if prior == 'Uniform':
                priors[i] = 1/(bounds[i,1]-bounds[i,0])
            if prior == 'Normal':
                priors[i] = 1/sigma[i]/2.5066 * np.exp(-(x[i]-mu[i])*(x[i]-mu[i])/2/sigma[i]/sigma[i])
    abundsum = np.sum(10**params[g1:g2])
    if abundsum>1.0:
        print('Sum of abundances > 1\n')
        return -np.inf
    
    grav = 0.
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        grav = 10**params[b1+pos]
    else: grav = 10**4.1
    
    if mass<minmass or mass>maxmass:
        print('Mass out of Bound. Rad={0} log(g)={1} Mass={2}'.format(radius/REarth/RJup,np.log10(grav),mass))
        return -np.inf

    if len(gases)==0:
        mmw = 2.28
    else:
        mmw,rxsec = af.GetScaOpac(gases,params[g1:g2])
    scale = 1.38e-16*teff/np.sum(mmw)/grav;
    if scale/radius > 0.05:
        print('Gravity too low for reliable convergence. teff={0}, mu={1}, log(g)={2}, scale={3}, rad={4}, ratio={5}'.format(teff,np.sum(mmw),np.log10(grav),scale/1.e5,radius/1.e5,scale/radius))
        return -np.inf
    if np.isnan(teff):
        print('Teff returned nan.')
        return -np.inf
    
    penalty = 0
    if smooth:
        gamma = params[a2]
        for i in range(a1+1,a2-1):
            penalty = penalty - 0.5/gamma*(params[i+1] - 2*params[i] + params[i-1])**2
        penalty = penalty - 0.5*np.log(2*np.pi*gamma)*(a2-a1)

    # ada: Add a check to make sure the Madhusudhan-Seager T-P parametrization has correct monotonicity.
    if atmtype == "power_linear":
        if (params[a2-3] < 0) or (params[a2-2] < 0) or (params[a2-1] < 0):
            print("Negative temperatures in one or more proposed parameters.")
            return -np.inf
        # ada: Removing the condition that the middle temperature node be >= the TOA temperature means we can have an inversion layer.
        elif params[a2-2] < params[a2-3]:
            print("Prior failed: profile set up for monotonically increasing temperature structure only. T_mid={}K < T_TOA={}K.".format(params[a2-2],params[a2-3]))
            return -np.inf
        elif params[a2-1] < params[a2-2]:
            print("Prior failed: profile set up for monotonically increasing temperature structure only. T_max={}K < T_mid={}K.".format(params[a2-1],params[a2-2]))
            return -np.inf

    # These lines weight the parameters based on the width of the prior if the boundaries cut off a significant amount of the normal distribution
    # However, it's not clear how important they are to comparing relative likelihoods
    #priors[1] = priors[1]/0.520
    #for i in range(9,14):
    #    priors[i] = priors[i]/0.9772
    
    return np.log(np.prod(priors))+penalty

# End of prior function

#----------------------------------------------------------------------------------------#
# Probability function

def lnprob(x,binslo,binshi,fluxrat,frathigh):
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]
        
    like, derived = lnlike(x,binslo,binshi)
    lp = lnprior(x,derived)
    
    # Check if any of the priors were out of bounds.
    if not np.isfinite(lp):
        return -np.inf, derived

    # Check if an error returned a non-result.
    prob = lp + like
    
    if np.isnan(prob):
        return -np.inf, derived
        
    return prob, derived

# End of probability function

#----------------------------------------------------------------------------------------#
# Sets up the runtime plot in lnlike(). Part of the testing suite.
'''
figa = plt.figure()
ax = figa.add_subplot(111)

x1=[0]
y1=[0]
x2=[0]
y2=[0]
li1, = ax.plot(x1,y1)
li2, = ax.plot(x2,y2)
figa.canvas.draw()
plt.show(block=False)
'''

#----------------------------------------------------------------------------------------#
# Set up the MCMC run

if task=='Ensemble' or task=='Spectral_Range':

    print(eplist)
    print(len(eplist))
    
    allspec = np.zeros((len(eplist),len(modwave)))
    for ii in range(0,len(eplist)):
        print('Model #{0:d}'.format(ii))
        radfinal = RJup

        pos = 0
        theta_planet = 0.
        if 'Rad' in basic:
            pos = basic.index('Rad')
            theta_planet = eplist[ii][b1+pos]*REarth/dist/parsec
            radfinal = eplist[ii][b1+pos]
        else:
            theta_planet = RJup*REarth/dist/parsec
            # Default radius = Jupiter
            
        # Probably need to handle deltaL here.
        print(eplist[ii])
        spectrum = GetModel(eplist[ii])
        
        # Multiply by solid angle and collecting area
        fincident = np.zeros(len(spectrum))
        if mode<=1:
            for j in range(0,len(spectrum)):
                fincident[j] = spectrum[j] * theta_planet*theta_planet
                #if i==0: print "newtdepth: ",j,specwave[i],fincident[i] # print output for debugging purposes
                # erg/s/aperture/Hz
                # theta_planet is actually the radius/distance ratio
                # so its square converts flux at the surface to flux at the telescope
        if mode==2:
            fincident = spectrum

        allspec[ii] = fincident

    for i in range(0,len(modwave)-1):
        fout.write('{0:8.5f} {1:8.5f}'.format(modwave[i],modwave[i+1]))
        for j in range(0,len(eplist)):
            fout.write(' {0:8.5e}'.format(allspec[j][i]))
        fout.write('\n')
        
    sys.exit()

# End of ensemble-specific section

#----------------------------------------------------------------------------------------#
# Set up retrieval

if task=='Retrieval':
    # Used to test the serial part of the code at the command line
    print('Reduced log-likelihood of input parameters: {0:f}'.format(lnprob(guess,ibinlo,ibinhi,binflux,binerr)[0]/len(ibinlo)))
    #sys.exit()

    eps = 0.01
    for i in range(0,len(guess)):
        if guess[i] < bounds[i,0]+eps*(mu[i]-bounds[i,0]): guess[i] = bounds[i,0]+eps*(mu[i]-bounds[i,0])
        if guess[i] > bounds[i,1]-eps*(bounds[i,1]-mu[i]): guess[i] = bounds[i,1]-eps*(bounds[i,1]-mu[i])
        
    if parallel:
        # MPI Setup
        from schwimmbad import MPIPool
        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
            
    pos = [guess + 0.1*eps*guess*np.random.randn(ndim) for i in range(nwalkers)]

    # Part of planned Multi-Nested Sampling capability
    #fprog = open('chain.dat','w')
    #fprog.close()
    #result = pymultinest.solve(LogLikelihood=lnlike, Prior=lnprior, n_dims=ndim, outputfiles_basename=pmntest, verbose=True)

    reader = emcee.backends.Backend()
    
    # MPI sampler
    if parallel:
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,backend=reader,pool=pool,args=(ibinlo,ibinhi,binflux,binerr))
    
    # Non-MPI sampler for testing purposes
    if not parallel:
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,backend=reader,args=(ibinlo,ibinhi,binflux,binerr))
    
    print('Walking...{0} {1} {2}...'.format(nwalkers,pllen,nsteps))

    if printfull: foutname = outdir + outfile + 'full.dat'
    else: foutname = outdir + outfile + 'dat'
    fchain = open(foutname,'w')

    if printfull: fchain.write('Samples {0:d} {1:d} {2:d} '.format(nwalkers,nsteps,ndim+5))
    else: fchain.write('Samples {0:d} {1:d} {2:d} '.format(nwalkers,(int)(nsteps/10),ndim+5))

    if atmtype=='Layers': fchain.write('Layered {0:d} {1:f} {2:f}\n'.format(a2-a1,minP-6.,maxP-6.))
    else: fchain.write('Parametric\n')

    for i in range(0,len(finalnames)):
        fchain.write('{0:s} '.format(finalnames[i]))
    fchain.write('Mass C/O [Fe/H] Teff Likelihood\n')

    #----------------------------------------------------------------------------------------#
    # Actual MCMC run, and Write samples to output file

    maxlikli = 0.
    medianparams = guess
    nsteps = 50
    i = 0
    for sample in sampler.sample(pos, iterations=nsteps):
        coords = reader.get_last_sample().coords
        blobs = reader.get_blobs()[i]
        likli = reader.get_log_prob()[i]
        if max(likli) > maxlikli or i==0:
            j = np.argmax(likli)
            maxlikli = likli[j]
            medianparams = coords[j]
        if i%100==0 or printfull: print('Step number {0:d}'.format(i+1))
        if printfull or i>=nsteps*0.9:
            for i2 in range(0,len(coords)):
                for i3 in range(0,len(coords[i2])):
                    fchain.write('{0:f} '.format(coords[i2][i3]))
                for i3 in range(0,len(blobs[i2])):
                    fchain.write('{0:f} '.format(blobs[i2][i3]))
                fchain.write('{0:f}\n'.format(likli[i2]/len(masternorm)))
        i = i+1
    fchain.close()

    finalparams = medianparams
    
    # End of MCMC run

    if parallel: pool.close()

    xplot = np.linspace(1,nsteps,nsteps)
    first = int(0.9*len(sampler.chain[0]))
    if printfull: first=0

    print("Retrieval Complete")

    #----------------------------------------------------------------------------------------#
    # Create waterfall plots of results

    gn1i = [i for i in nvars if g1<=i]
    if len(gn1i)>0: gn1 = nvars.index(gn1i[0])
    else: gn1 = -1
    gn2i = [i for i in nvars if g2<=i]
    if len(gn2i)>0: gn2 = nvars.index(gn2i[0])
    else: gn2 = -1
    gnames2 = []
    if gn2<gn1:
        gsamples = sampler.chain[:,first:,gn1:]
        for i in range(gn1,len(pnames)): gnames2.append(pnames[nvars[i]])
    else:
        gsamples = sampler.chain[:,first:,gn1:gn2]
        for i in range(gn1,gn2): gnames2.append(pnames[nvars[i]])

    bn1i = [i for i in nvars if b1<=i]
    if len(bn1i)>0: bn1 = nvars.index(bn1i[0])
    else: bn1 = -1
    bn2i = [i for i in nvars if b2<=i]
    if len(bn2i)>0: bn2 = nvars.index(bn2i[0])
    else: bn2 = -1    
    bnames2 = []
    if bn2<bn1:
        bsamples = sampler.chain[:,first:,bn1:]
        for i in range(bn1,len(pnames)): bnames2.append(pnames[nvars[i]])
    else:
        bsamples = sampler.chain[:,first:,bn1:bn2]
        for i in range(bn1,bn2): bnames2.append(pnames[nvars[i]])

    an1i = [i for i in nvars if a1<=i]
    if len(an1i)>0: an1 = nvars.index(an1i[0])
    else: an1 = -1
    an2i = [i for i in nvars if a2<=i]
    if len(an2i)>0: an2 = nvars.index(an2i[0])
    else: an2 = -1
    if an2<an1: tsamples = sampler.chain[:,first:,an1:]
    else: tsamples = sampler.chain[:,first:,an1:an2]
    
    dsamples = np.transpose(sampler.blobs[first:],axes=[1,0,2])

    gsamples2 = gsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),gn2-gn1),order='F')
    bsamples2 = bsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),bn2-bn1),order='F')
    tsamples2 = tsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),an2-an1),order='F')
    dsamples2 = dsamples.reshape((len(dsamples[0])*len(dsamples),4),order='F')

    bnames2.append('Mass')
    bnames2.append('C/O')
    bnames2.append('[Fe/H]')
    bnames2.append('Teff')

    baselist = ['Rad','Log(g)','Cloud_Base','P_cl','Mass','C/O','[Fe/H]','Teff']
    bnamelist = ['Radius (R$_J$)','log(g)','P$_{cloud}$ (bar)','Base Pressure (bar)','Mass (M$_J$)','C/O','[Fe/H]','T$_{eff}$ (K)']
    gaslist = ['h2he','h2','he','h-','h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
    gnamelist = ['H$_2$+He','H$_2$','He','[H-]','[H$_2$O]','[CH$_4$]','[CO]','[CO$_2$]','[NH$_3$]','[H$_2$S]','[Na,K]','[Na,K]','[CrH]','[FeH]','[TiO]','[VO]','[HCN]','[N2]','[PH3]']
    bnames = []
    gnames = []

    for i in range(0,len(bnames2)):
        if bnames2[i] in baselist:
            j = baselist.index(bnames2[i])
            bnames.append(bnamelist[j])
    
    for i in range(0,len(gnames2)):
        if gnames2[i] in gaslist:
            j = gaslist.index(gnames2[i])
            gnames.append(gnamelist[j])

    rpos = -1
    ppos = -1
    if 'P_cl' in basic:
        ppos = basic.index('P_cl')
    if 'Cloud_Base' in basic:
        ppos = basic.index('Cloud_Base')
    
    bsamples3 = np.zeros((len(bsamples2),bn2-bn1+4))
    lenbasic = len(bsamples2[0])

    for i in range(0,len(bsamples3)):
        for j in range(0,lenbasic):
            if j==rpos:
                bsamples3[i,j] = bsamples2[i,j]
            if j==ppos:
                bsamples3[i,j] = bsamples3[i,j]-6.
        for j in range(0,4):
            bsamples3[i,lenbasic+j] = dsamples2[i,j]
            
    igood = [i for i in range(0,len(bsamples3)) if not np.isinf(bsamples3[i,-1])]
    gsamples3 = gsamples2[igood]
    bsamples4 = bsamples3[igood]
    tsamples3 = tsamples2[igood]

    grange = np.zeros(len(gnames))
    for i in range(0,len(gnames)): grange[i]=0.99
    brange = np.zeros(len(bnames))
    for i in range(0,len(bnames)): brange[i]=0.99
    
    if plotting:
        fig = corner.corner(gsamples3,labels=gnames,range=grange,plot_datapoints=False,labelsize=24)
        fig.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
        fig1name = 'plots' + outfile + 'gases.png'
        fig.savefig(fig1name)
        
        fig2 = corner.corner(bsamples4,labels=bnames,range=brange,plot_datapoints=False,labelsize=24)
        fig2.subplots_adjust(left=0.10,bottom=0.10,wspace=0,hspace=0)
        fig2name = 'plots' + outfile + 'basic.png'
        fig2.savefig(fig2name)
        
        # Plot the T-P profile
        
        plist = np.zeros(anum)
        for i in range(0,anum): plist[i] = 10**(maxP + (minP-maxP)*i/(anum-1)) * 1.e-6
        tlist = np.percentile(tsamples3,[16,50,84],axis=0)
        
        fig3 = plt.figure(figsize=(8,6))
        ax = fig3.add_subplot(111)
        plt.axis((0,3000,10**(maxP-6.),10**(minP-6.)))
        ax.set_yscale('log')
        plt.xlabel('T (K)',fontsize=14)
        plt.ylabel('P (bar)', fontsize=14)
        
        ax.fill_betweenx(plist,tlist[0],tlist[2],facecolor='#ff8080')
        ax.plot(tlist[0],plist,c='r')
        ax.plot(tlist[1],plist,c='k')
        ax.plot(tlist[2],plist,c='r')
        
        fig3name = 'plots' + outfile + 'TP.png'
        fig3.savefig(fig3name)
    
    # End of retrieval plots

    #----------------------------------------------------------------------------------------#
    # Write parameter file of best fit model
    finallower, medianparams_all, finalupper = np.percentile(sampler.chain[:,first:,:],[16,50,84],axis=0)[:,0]

    finalparams2 = np.percentile(sampler.chain[:,first:,:],[50,84],axis=0)[:,0]

    finalsigma2 = np.zeros(ndim)
    for i in range(0,ndim):
        finalsigma2[i] = (float)(finalparams2[1][i]) - (float)(finalparams2[0][i])

    finalparams = np.zeros(pllen)
    finalsigma = np.zeros(pllen)
    finalbounds = np.zeros((pllen,2))
    for i in range(0,pllen):
        if i in nvars:
            j = nvars.index(i)
            finalparams[i] = finalparams2[0][j]
            finalsigma[i] = finalsigma2[j]
            finalbounds[i,0] = bounds[j,0]
            finalbounds[i,1] = bounds[j,1]
        else:
            finalparams[i] = plparams[i]
            finalsigma[i] = 0.
            finalbounds[i,0] = finalparams[i]
            finalbounds[i,1] = finalparams[i]
        
    outparams = '.' + outfile + 'retrieved.dat'
    ffout = open(outparams,'w')
    
    ffout.write('Mode         {0:s}\n'.format(modestr))
    ffout.write('Object       {0:s}\n'.format(name))
    ffout.write('Star         {0:5.0f} {1:5.2f} {2:8.3f}\n'.format(tstar,rstar,sma))
    if not starspec=='': ffout.write('Star_Spec   {0:s}\n'.format(starspec))
    ffout.write('Location     {0:6.2f} {1:6.2f} {2:6.2f}\n'.format(dist,RA,dec))
    ffout.write('Data         {0:s}\n'.format(datain))
    ffout.write('Convolve     {0:5.1f}\n'.format(dataconv))
    ffout.write('Binning      {0:5.1f}\n'.format(databin))
    ffout.write('Degrade      {0:5.1f}\n'.format(degrade))
    ffout.write('Pressure     {0:5.1f} {1:5.1f}\n'.format(minP-6.,maxP-6.))
    ffout.write('Streams      {0:d}\n'.format(streams))
    ffout.write('Vres         {0:d}\n'.format(vres))
    if gray: ffout.write('Gray        {0:5.0f}\n'.format(tgray))
    ffout.write('N_Steps         {0:d}\n'.format(nsteps))
    ffout.write('Parallel        {}\n'.format(parallel))
    ffout.write('Prior           {0:s}\n'.format(prior))
    ffout.write('Mass_Limits {0:5.2f} {1:5.2f}\n'.format(minmass,maxmass))
    ffout.write('Tables       {0:s} {1:s}\n'.format(hires,lores))
    ffout.write('Output       modelspectra    Short\n')
    ffout.write('Opacities    {0:s}\n'.format(opacdir))
    if not outmode=='': ffout.write('Output_Mode {0:s}\n'.format(outmode))
    
    ffout.write('Parameter    Initial    Mu    Sigma    Min    Max\n')
    if b1>=0:
        ffout.write('Basic\n')
        for i in range(b1,b2):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if g1>=0:
        ffout.write('Gases     {0:s}\n'.format(gases[0]))
        for i in range(g1,g2):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if a1>=0:
        ffout.write('Atm       {0:s}\n'.format(atmtype))
        if smooth:
            for i in range(a1,a2+1):
                ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))
        else:
            for i in range(a1,a2):
                ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))            
                
    if c1>=0:
        ffout.write('Clouds    {0:d}    {1:s}\n'.format(cloudmod,hazestr))
        for i in range(c1,c2):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if e1>=0:
        ffout.write('End\n')
        for i in range(e1,e2):
            ffout.write('{0:s}    {1:8.4f}    {2:8.4f}    {3:8.4f}    {4:8.4f}    {5:8.4f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

# End of retrieval-specific section
    
#----------------------------------------------------------------------------------------#
# Plot the best fit spectrum or spectrum requested at the command line
    
if task=='Spectrum':
    finalparams = plparams

radfinal = RJup

while True:
    pos = 0
    theta_planet = 0.
    if 'Rad' in basic:
        pos = basic.index('Rad')
        theta_planet = finalparams[b1+pos]*REarth/dist/parsec
        radfinal = finalparams[b1+pos]
    else:
        theta_planet = RJup*REarth/dist/parsec
        # Default radius = Jupiter
        
    spectrum, derived = GetModel(finalparams)
    
    # Multiply by solid angle and collecting area
    fincident = np.zeros(len(spectrum))
    if mode<=1:
        for i in range(0,len(spectrum)):
            fincident[i] = spectrum[i] * theta_planet*theta_planet
            #if i==0: print("newtdepth: ",i,spectrum[i],fincident[i]) # print output for debugging purposes
            # erg/s/aperture/Hz
            # theta_planet is actually the radius/distance ratio
            # so its square converts flux at the surface to flux at the telescope
    if mode==2:
        fincident = spectrum
        
    # Adjust for wavelength calibration error
    deltaL = np.zeros(nband)
    
    for n in range(0,enum):
        if end[n][0:6]=='deltaL':
            if end[n]=='deltaL': pos = 0
            else: pos = int(end[n][6:])
            deltaL[pos] = finalparams[e1+n]
            
    newibinlo = np.zeros(len(ibinlo))
    newibinhi = np.zeros(len(ibinhi))

    di = 0
    for i in range(0,nband):
        for j in range(di,di+len(bandlo[i])):
            newibinlo[j] = ibinlo[j] + deltaL[i] * delibinlo[j]
            newibinhi[j] = ibinhi[j] + deltaL[i] * delibinhi[j]
        di = di + len(bandlo[i])
    
    binw = (newibinlo[1]-newibinlo[0])*(dataconv/databin)
    convmod = af.ConvSpec(fincident,binw)
    binmod = af.BinModel(convmod,newibinlo,newibinhi)
    resid = binflux-binmod
    specout = binmid
    
    if plotting:
        xmin = min(specout)
        xmax = max(specout)
        xmin = xmin - 0.05*(xmax-xmin)
        xmax = xmax + 0.05*(xmax-xmin)
        
        yref = max(max(binmod),max(binflux))
        ymin = -0.20 * yref
        ymax =  1.05 * yref
        #ymin = 0.97*yref  # I think this is for transits.
        #ymax = 1.02*yref
        
        # Plot the BINNED model/retrieved spectrum against the observations.
        
        fig4 = plt.figure(figsize=(10,7))
        ax = fig4.add_subplot(111)
        plt.axis((xmin,xmax,ymin,ymax))
        
        plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
        plt.ylabel('Flux (cgs)',fontsize=14)
        plt.tick_params(axis='both',which='major',labelsize=12)
        
        ax.errorbar(specout,binflux,binerr,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k',zorder=-1)
        ax.plot(specout,binmod,'-',linewidth=1,label='Retrieved Spectrum',c='#8080ff',zorder=1)
        ax.plot(specout,resid+ymin/2.,'-',linewidth=1,label='Residuals (offset)',c='r')
        ax.plot([xmin,xmax],[0.,0.],'-',c='k')
        ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')
        
        plt.legend(fontsize=12)
        
        if manual:
            plt.show()
            
    if not manual:
        print('Computing final outputs.')
        break

    prompt = input('Do you want to compute a new spectrum (y/n)? ')
    while not (prompt in ['y', 'Y', 'n', 'N']):
        prompt = input('Error. Please type \'y\' or \'n\': ')
    if prompt=='n' or prompt=='N':
        print('Computing final outputs.')
        break

    while prompt!='Go':
        prompt = input('Enter a parameter name and new value or type \'Go\' to compute the new spectrum: ')
        psplit = prompt.split()
        if prompt!='Go' and (len(psplit)<2 or psplit[0] not in pnames):
            print('Error: invalid input.')
        else:
            pos = pnames.index(psplit[0])
            finalparams[pos] = float(psplit[1])

if task=='Spectrum': outfile = '/' + name + '.Spectrum.'

if plotting:
    fig4name = 'plots' + outfile + 'binned.png'
    fig4.savefig(fig4name)

#----------------------------------------------------------------------------------------#
# Create an output file of the BINNED model/retrieved spectrum.

foutnameb = 'modelspectra' + outfile + 'binned.dat'
fout = open(foutnameb,'w')
for i in range(0,len(specout)):
    fout.write('{0:8.5f} {1:8.5f} {2:12.9e} 0.0 0.0 {2:12.9e}\n'.format(binlo[i],binhi[i],binmod[i]))
fout.close()

# Plot the FULL-RES model/retrieved spectrum against the observations.

if plotting:
    fig5name = 'plots' + outfile + 'fullres.png'
    
    fig5 = plt.figure(figsize=(10,7))
    ax = fig5.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))
    
    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)

# Compute the residuals by binning to the observations without downsampling for resolving power.
wavemid = (wavelo+wavehi)/2.
binw = (wavelo[1]-wavelo[0])
convmod = af.ConvSpec(fincident,binw)
newbinslo,newbinshi = af.GetBins(modwave,wavelo,wavehi)
binmod = af.BinModel(convmod,newbinslo,newbinshi)

convflux2 = convflux[0]
converr2 = converr[0]
for i in range(1,len(convflux)):
    convflux2 = np.r_[convflux2,convflux[i]]
    converr2 = np.r_[converr2,converr[i]]
resid2 = convflux2-binmod
    
rechisq = np.sum(resid2**2/binerr**2)/(len(convflux2)-len(pvars))
rmserr = np.sqrt(np.sum(resid**2)/len(convflux2))/np.mean(binerr)
print('Reduced chi-square (no correction):\t\t{0:f}'.format(rechisq))
print('RMS Error (sigma, no correction):\t\t{0:f}'.format(rmserr))

lnf = np.zeros(nband)
for n in range(0,enum):
    if end[n][0:4]=='logf':
        if end[n]=='logf': pos = 0
        else: pos = int(end[n][4:])
        lnf[pos] = finalparams[e1+n]
        print(finalparams[e1+n])

# This plots the scaled error bars, but it's passed over to avoid confusion.
'''
di = 0
for i in range(0,nband):
    for j in range(di,di+len(bandlo[i])):
        binerr[j] = np.sqrt(binerr[j]*binerr[j] + np.exp(lnf[i]))
    di = di + len(bandlo[i])
''' 
rechisq = np.sum(resid2**2/binerr**2)/(len(convflux2)-len(pvars))
rmserr = np.sqrt(np.sum(resid**2)/len(convflux2))/np.mean(binerr)
print('Reduced chi-square (error bar correction):\t{0:f}'.format(rechisq))
print('RMS Error (sigma, error bar correction):\t{0:f}'.format(rmserr))

print(binerr)

if plotting:
    ax.errorbar(specout,binflux,binerr,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k',zorder=-1)
    # I don't know what this was supposed to be, but it doesn't work.
    #ax.errorbar(wavemid,convflux2,converr2,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k',zorder=-1)
    ax.plot(modwave,fincident,'-',linewidth=1,label='Retrieved Spectrum',c='#8080ff',zorder=1)
    ax.plot(wavemid,resid2+ymin/2.,'-',linewidth=1,label='Residuals (convolved and offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')
    
    plt.legend(fontsize=12)
    fig5.savefig(fig5name)

# Create an output file of the FULL-RES model/retrieved spectrum.

foutnamef = 'modelspectra' + outfile + 'fullres.dat'
fout = open(foutnamef,'w')
for i in range(0,len(modwave)-1):
    fout.write('{0:8.5f} {1:8.5f} {2:12.9e} 0.0 0.0 {2:12.9e}\n'.format(modwave[i],modwave[i+1],fincident[i]))

# End of plot spectrum

# Calling getContribution, which returns "taulayer" from the C++ side.
contribution = planet.getContribution()
cloudcont = planet.getCloudContribution()
gascont = planet.getGasContribution()
speciescont = planet.getSpeciesContribution()

nlayers = np.shape(contribution)[-1] + 1
logPs = (minP-6.0) + (maxP-minP)*np.arange(nlayers)/(nlayers-1)
mean_logPs = (logPs[:-1]+logPs[1:])/2.

def bin_contributions(contribution, wavelengths=modwave, datain=datain, logPs=mean_logPs):
    # Continue if an observation file was named
    if len(datain)>0:

        # Output binned to the observations
        fobs = open(datain,'r')

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

        interpolate = interp1d(wavelengths, contribution, kind="linear", axis=-2)
        binned_wavelengths = rcalmid
        binned_contribution = interpolate(rcalmid)

    else:
        binned_wavelengths = wavelengths
        binned_contribution = interpolate(rcalmid)

    return pd.DataFrame(data=binned_contribution,
                        index=binned_wavelengths,
                        columns=mean_logPs)

speciescont = {gas: bin_contributions(cont)
               for gas, cont in zip(gases, speciescont)}

otherconts = {"total": bin_contributions(contribution),
              "cloud": bin_contributions(cloudcont),
              "gas": bin_contributions(gascont)}

contributions = speciescont.copy()
contributions.update(otherconts)

with open(name+"_contributions.p", "wb") as pickle_file:
    pickle.dump(contributions,
                pickle_file,
                protocol=pickle.HIGHEST_PROTOCOL)

if task=='Retrieval': sys.exit()

#----------------------------------------------------------------------------------------#
# Create binned files for particular JWST modes.
# This is only allowed in Spectrum mode because Retrieval mode doesn't compute the full opacity table for efficiency.
# The partial opacity table often does not cover the JWST modes.
# This should not be a significant issue because this section is used for modeling observations rather than retrievals.

# Set noise parameters
noise_params = np.zeros(7)

noise_params[0] = radfinal      # radius
noise_params[1] = rstar         # R_star (R_Sun)
noise_params[2] = dist          # Distance (pc)
noise_params[3] = RA            # Right Ascension
noise_params[4] = dec           # Declination
noise_params[5] = 6500.         # T_star
noise_params[6] = exptime/3600. # exposure time in hours

noisemode = -1

modenames = ['G140H-F070LP','G140H-F100LP','G235H-F170LP','G395H-F290LP','G140M-F070LP','G140M-F100LP','G235M-F170LP','G395M-F290LP','Prism-Clear','LRS_Slit','LRS_Slitless','NIRCam','NIRISS','MRS_A','MRS_B','MRS_C']
modefiles = ['NIRSpec1H','NIRSpec2H','NIRSpec3H','NIRSpec4H','NIRSpec1M','NIRSpec2M','NIRSpec3M','NIRSpec4M','Prism','MIRISlit','MIRISlitless','MIRIA','MIRIB','MIRIC']

if outmode == '': sys.exit()
if outmode in modenames:
    noisemode = modenames.index(outmode)
else:
    if os.path.exists(outmode):
        calwave,flux_density,fnoise = AddNoise.addNoise(noisemode,mode,opacwave,spectrum,noise_params,starspec,outmode)
        print('Band Center: {0:8.2f} cm^-1'.format(calwave))
        print('Flux:     {0:12.5e} cgs'.format(flux_density))
        print('Noise:    {0:12.5e} cgs'.format(fnoise*flux_density))
    else:
        print('Error: filter file not found.')
        
if noisemode >= 0:
    if (noisemode < 8 or noisemode == 11 or noisemode == 12) and (modwave[0] < 5.01):
        print('Requested spectral mode does not match input wavelengths.')
        sys.exit()
    if (noisemode > 12 or noisemode == 9 or noisemode == 10) and (modwave[-1] > 4.99):
        print('Requested spectral mode does not match input wavelengths.')
        sys.exit()
        
    calwave,flux_density,fnoise = AddNoise.addNoise(noisemode,mode,opacwave,spectrum,noise_params,starspec)
        
    callo = np.zeros(len(calwave))
    calhi = np.zeros(len(calwave))
    callo[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    if callo[0] < modwave[0]: callo[0] = modwave[0] + deltaL/1000.
    calhi[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    if calhi[-1] > modwave[-1]: calhi[-1] = modwave[-1] + deltaL/1000.
    
    for i in range(0,len(calwave)-1):
        callo[i+1] = (calwave[i]+calwave[i+1])/2. + deltaL/1000.
        calhi[i]   = (calwave[i]+calwave[i+1])/2. + deltaL/1000.
        
    obsdepth = flux_density
    if mode<=1: noise = fnoise*obsdepth
    if mode==2: noise = fnoise
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise

    foutname = 'modelspectra' + outfile + outmode + '.dat'
    ftest = open(foutname,'w')
    
    for i in range(0,len(calwave)-1):
        ftest.write('{0:8.5f} {1:8.5f} {2:12.9e} {3:12.9e} {3:12.9e} {4:12.9e}\n'.format(callo[i]-deltaL/1000.,calhi[i]-deltaL/1000.,obsdepth[i],noise[i],obs_flux[i]))

    # Plot the JWST mode against the observations.
        
    plotname = 'plots' + outfile + outmode  + '.fit.png'  # Need to put this in the plots folder.
    
    fig2 = plt.figure(figsize=(10,7))
    ax = fig2.add_subplot(111)

    snip = [i for i in range(0,len(modwave)) if calwave[0] > modwave[i] < calwave[-1]]
    
    xmin = min(calwave)
    xmax = max(calwave)
    xmin = xmin - 0.05*(xmax-xmin)
    xmax = xmax + 0.05*(xmax-xmin)
    
    yref = max(max(obs_flux),max(fincident[snip]))
    ymin = -0.20 * yref
    ymax =  1.05 * yref
    
    plt.axis((xmin,xmax,ymin,ymax))

    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)
    
    ax.plot(modwave[snip],fincident[snip],'-',linewidth=1,label='Model Spectrum',c='k')
    ax.errorbar(calwave,obsdepth,noise,capsize=3,marker='o',linestyle='',linewidth=1,label=outmode,c='b')
    # Alternate plot as a line without errorbars.
    #ax.plot(calwave,obsdepth,'-',linewidth=1,label=outmode,c='b')
    
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    
    plt.legend(fontsize=12)
    plt.savefig(plotname) # Need to procedurally generate filename
