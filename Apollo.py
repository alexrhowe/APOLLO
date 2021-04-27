from __future__ import print_function
import os
import sys
import numpy as np
import scipy.optimize as op
from scipy.interpolate import interp1d
from scipy.stats import invgamma
from distutils.util import strtobool
import emcee
import schwimmbad

import matplotlib
# Used on Discover because the GUI does not work there. Comment out to use 'Manual' mode.
#matplotlib.use('pdf')
import matplotlib.pyplot as plt

import corner
from src import wrapPlanet
from src import ApolloFunctions as af
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
ApolloFunctions.py
AddNoise.py
Filter.py
setup.py

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

Various example files included
'''

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

# Default Settings

name = 'example'     # Bundled example file
mode = 0             # Emission spectrum
modestr = 'Resolved' # Needed for output file name
parallel = True      # Parallel operation
datain = 'examples/example.obs.dat' # Bundled example file
polyfit = False      # Switch to normalize the spectrum to a polynomial fit
norm = False         # Dummy variable if polyfit is false
dataconv = 1         # Factor to convolve the observations to blur the spectrum or account for actual resolving power
databin = 1          # Factor to bin down the observations for simpler fitting
degrade = 1          # Factor to degrade the model spectrum for faster calculation
prior = 'Uniform'    # Uniform priors
nwalkers = 0         # Placeholder for later adjustment
nsteps = 30000       # Tested minimum required number of steps
tstar = 5770.        # Solar temperature
rstar = 1.0          # Solar radius
sma = 1.0            # Semi-Major Axis
starspec = ''        # Optional stellar spectrum file
dist = 10.0          # Standardized distance, 10 pc
RA = 0.0             # Right ascension
dec = 0.0            # Declination
minmass = 0.5        # Minimum mass in Jupiter masses
maxmass = 80.0       # Maximum mass in Jupiter masses
hires = ''           # Default set of opacity tables to compute the spectra
lores = 'lores'      # Default low-resolution tables to compute Teff
minP = 0.0           # Pressure range to integrate over in cgs, default 1 mubar to 1 kbar.
maxP = 9.0
gray = False         # Used to create a gray atmosphere for testing
vres = 71            # Number of layers for radiative transfer
streams = 1          # Use 1-stream by default
wavei = 10000./0.60  # Full NIR wavelength range
wavef = 10000./5.00
outmode = ''         # JWST observing mode
exptime = 1000.      # Exposure time in seconds
outdir = 'samples'   # Default output directory
short = False        # Switch to create short output file names
printfull = False    # Switch to print the full sample array instead of the last 10%
opacdir = '../Opacities' # Default opacities directory

norad = False        # Flags if no radius variable is in the input file
natm = 0             # Placeholder in case T-P profile is omitted
verbatim = False     # Interpolate the T-P profile
tgray = 1500         # Temperature of gray atmosphere
hazetype = 0         # No Clouds
hazestr = 'None'     # No Clouds
cloudmod = 0         # No Clouds

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

    elif line[0]=='Object':
        if len(line)>1: name = line[1]
    if line[0]=='Mode':
        if len(line)>1: modestr = line[1]
        if modestr=='Resolved': mode = 0
        if modestr=='Eclipse': mode = 1
        if modestr=='Transit': mode = 2
    elif line[0]=='Parallel':
        if len(line)>1: parallel = strtobool(line[1])
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

# Special handling if area ratio is used instead of radius
if 'RtoD2U' in basic:
    pos = basic.index('RtoD2U')
    plparams[b1+pos] = 10**plparams[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2
    guess[b1+pos] = 10**guess[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2    
    sigma[b1+pos] = guess[b1+pos]*(10**sigma[b1+pos]-1.)
    mu[b1+pos] = 10**mu[b1+pos] * dist**2 * 4.838e9**2
    sigma[b1+pos] = sigma[b1+pos]*mu[b1+pos]
    bounds[b1+pos] = 10**bounds[b1+pos] * dist**2 * 4.838e9**2
    
# Meant to make the temperatures uniform in log space, not currently used
'''
if atmtype == 'Layers':
    print 'Layers\n'
    plparams[a1:a2] = np.log10(plparams[a1:a2])
    guess[a1:a2] = np.log10(guess[a1:a2])
    bounds[a1:a2,:] = np.log10(bounds[a1:a2,:])
'''

# End of read in model parameters
# End of read in input file

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

# Read in observations
# Note: header contains info about star needed for JWST pipeline

print('Reading in observations.')
fobs = open(datain,'r')

obslines = fobs.readlines()
obslength = len(obslines)

wavehi = np.zeros(obslength)
wavelo = np.zeros(obslength)
flux = np.zeros(obslength)
errlo = np.zeros(obslength)
errhi = np.zeros(obslength)

for i in range(0,obslength):
    wavehi[i] = obslines[i].split()[0]
    wavelo[i] = obslines[i].split()[1]
    flux[i] = obslines[i].split()[5]
    errlo[i] = obslines[i].split()[3]
    errhi[i] = obslines[i].split()[4]

# End of read in observations

# Set statistical parameters    
if 'logf' in end:
    pos = end.index('logf')
    plparams[e1+pos] = np.log(max(errhi**2))
    guess[e1+pos] = guess[e1+pos]
    mu[e1+pos] = np.log(max(errhi**2))
    sigma[e1+pos] = abs(mu[e1+pos])/10.
    bounds[e1+pos,0] = np.log(min(errhi**2) * bounds[e1+pos,0])
    bounds[e1+pos,1] = np.log(max(errhi**2) * bounds[e1+pos,1])

# Set the cross section tables if not already set.
wavei = max(wavehi) * 1.001
wavef = min(wavelo) * 0.999
if hires=='':
    if wavei > 10000/5.0 and wavef > 10000/5.0:
        if wavei > 10000/0.6: wavei = 10000./0.6
        hires = 'nir'
    elif wavei > 10000/5.0 and wavef < 10000/5.0:
        if wavei > 10000/0.6: wavei = 10000./0.6
        if wavef < 10000/30.0: wavef = 10000./30.0
        hires = 'wide'
    elif wavei < 10000/5.0 and wavef < 10000/5.0:
        if wavef < 10000/30.0: wavef = 10000./30.0
        hires = 'mir'

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
    opacwave[i] = 10000./(lmin*np.exp(i*degrade/resolv))

if wavehi[0]>opacwave[0] or wavelo[-1]<opacwave[-1]:
    trim = [i for i in range(0,len(wavehi)) if (wavehi[i]<opacwave[0] and wavelo[i]>opacwave[-1])]
    wavehi = wavehi[trim]
    wavelo = wavelo[trim]
    flux = flux[trim]
    errlo = errlo[trim]
    errhi = errhi[trim]
    obslength = len(trim)

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
    modwavelo[i] = 10000./(lminlo*np.exp(i/resolvlo))
    
# Set up wavelength ranges
imin = np.where(opacwave<np.max(wavehi))[0]-1
imax = np.where(opacwave<np.min(wavelo))[0]+2
istart = np.where(opacwave<wavehi[0])[0]-1
iend = np.where(opacwave<wavelo[-1])[0]+2

if len(imin)==0: imin = [0]
elif imin[0]<0: imin[0] = 0
if len(imax)==0: imax = [opacwave[-1]]
elif imax[-1]>=len(opacwave): imax[-1] = opacwave[-1]
if len(istart)==0: istart = [0]
elif istart[0]<0: istart[0]=0
if len(iend)==0: iend = [opacwave[-1]]
elif iend[-1]>=len(opacwave): iend[-1] = opacwave[-1]

if len(imin)==0: imin=[0]
if len(imax)==0: imax=[len(opacwave)-1]

# Truncated and band-limited spectrum that does not extend beyond the range of the observations
modwave = opacwave[(int)(imin[0]):(int)(imax[0])]
lenmod = len(modwave)

# End set up model spectrum wavelength range

# Process observations for retrieval

# Separate out individual bands
bandindex,bandhi,bandlo,bandflux,banderr = af.FindBands(wavehi,wavelo,flux,errhi)
nband = len(bandhi)

# Convolve the observations to account for effective resolving power or fit at lower resolving power
convflux,converr = af.ConvBands(bandflux,banderr,dataconv)

bandhifl = [item for sublist in bandhi for item in sublist]
bandlofl = [item for sublist in bandlo for item in sublist]
convfluxfl = [item for sublist in convflux for item in sublist]
converrfl = [item for sublist in converr for item in sublist]

# Bin the observations to fit a lower sampling resolution
binhi,binlo,binflux,binerr = af.BinBands(bandhifl,bandlofl,convfluxfl,converrfl,databin)
binlen = len(binflux)
binmid = np.zeros(len(binhi))
for i in range(0,len(binhi)): binmid[i] = (binhi[i]+binlo[i])/2.

totalflux = 0
for i in range(0,len(binflux)): totalflux = totalflux + binflux[i]*(1./binhi[i]-binlo[i])

# Handle bands and optional polynomial fitting
bindex, modindex, modwave = af.SliceModel(bandhi,bandlo,modwave)
polyindex = -1
for i in range(1,len(bindex)):
    if bindex[i][0] < bindex[i-1][0]:
        polyindex = i        
if polyindex==-1: polyfit = False

if polyfit:
    normhi = bandhi[0]
    normlo = bandlo[0]
    normflux = bandflux[0]
    normerr = banderr[0]
    for i in range(1,len(bandhi)):
        normhi = np.r_[normhi,bandhi[i]]
        normlo = np.r_[normlo,bandlo[i]]
        normflux = np.r_[normflux,bandflux[i]]
        normerr = np.r_[normerr,banderr[i]]
    normmid = (normhi+normlo)/2.

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
    
# End of band handling

if task=='Spectrum' or task=='Ensemble': modwave = opacwave

# Get indices of the edges of the observation bins in the model spectrum
bins = af.GetBins(modwave,binhi,binlo)
ibinhi = bins[0]
ibinlo = bins[1]

# Needed to calculate the spectrum with the wavelength offset later.
delmodwave = modwave+0.001
delbins = af.GetBins(delmodwave,binhi,binlo)
delibinhi = ibinhi-delbins[0]
delibinlo = ibinlo-delbins[1]

mmw,rxsec = af.GetScaOpac(gases,plparams[g1:g2])
mollist = af.GetMollist(gases)

natm = a2-a1
profin = np.zeros(natm)
for i in range(0,natm): profin[i] = maxP + (minP-maxP)*i/(natm-1)

if atmtype == 'Parametric' and natm != 5:
    print('Error: wrong parameterization of T-P profile.')
    sys.exit()

# Create Planet and read in opacity tables
planet = wrapPlanet.PyPlanet()
print('Haze type:',hazestr)
print('Cloud model:',cloudmod)
mode = int(mode)
cloudmod = int(cloudmod)
hazetype = int(hazetype)

atmmod = 0
if atmtype=='Layers': atmmod = 0
if atmtype=='Parametric': atmmod = 1

switches = [mode,cloudmod,hazetype,streams,atmmod]

guess  = guess[nvars]
mu     = mu[nvars]
sigma  = sigma[nvars]
bounds = bounds[nvars]

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
    
    # Radius handling
    if 'Rad' in basic:
        pos = basic.index('Rad')
        params1[0] = x[b1+pos]
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        params1[0] = 10**x[b1+pos]*dist*4.838e9 # convert R/D to Earth radii
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        params1[0] = np.sqrt(x[b1+pos])
    else:
        global norad
        norad = True
        params1[0] = 11.2
        # Default radius = Jupiter

    # Gravity handling
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        params1[1] = x[b1+pos]
    else:
        params1[1] = 4.1

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
    
    params1[3] = tstar
    params1[4] = rstar
    params1[5] = mmw
    params1[6] = rxsec
    params1[7] = minP
    params1[8] = maxP
    params1[9] = sma
    
    if hazetype!=0:
        if cloudmod==2 or cloudmod==3:
            for i in range(0,4): params1[i+10] = x[c1+i]
            params1[11] = params1[11]
            params1[12] = params1[12] + 6.
        if cloudmod==2: params1[13] = params1[12] + params1[13]
        
    tpprof = np.zeros(natm)
    # Gray atmosphere approximation
    if gray:
        tplong = np.zeros(vres)
        for i in range(0,vres): tplong[i] = tgray
        planet.set_Params(params1,abund,tplong)
        specflux = planet.get_Spectrum()

    # Build atmosphere temperature profile and compute spectrum
    else:
        for i in range(0,len(tpprof)): tpprof[i] = x[i+a1]
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
            specflux = planet.get_Spectrum()
        if atmtype == 'Parametric':
            # Compute spectrum
            planet.set_Params(params1,abund,tpprof)
            specflux = planet.get_Spectrum()
            
    teff = planet.get_Teff()
    if task!='Ensemble': print('Teff: ',teff)
    '''
    figtest = plt.figure()
    ax = figtest.add_subplot()
    ax.plot(modwave,specflux)
    plt.savefig('plots/test.png')
    sys.exit()
    '''
    return specflux

# End of GetModel function

# Likelihood function for "Retrieval" mode.

def lnlike(x,ibinhi,ibinlo,binflux,binerrhi):
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]
        
    modflux = GetModel(params)
    
    theta_planet = 0.
    if 'Rad' in basic:
        pos = basic.index('Rad')
        theta_planet = params[b1+pos]*6.371e8/dist/3.086e18
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        theta_planet = 10**params[b1+pos]
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        theta_planet = np.sqrt(params[b1+pos])*6.371e8/dist/3.086e18
    else:
        theta_planet = 11.2*6.371e8/dist/3.086e18
        # Default radius = Jupiter
    
    # Statistical parameters
    if 'deltaL' in end:
        pos = end.index('deltaL')
        deltaL = params[e1+pos]
    else:
        deltaL = 0.0
    if 'logf' in end:
        pos = end.index('logf')
        lnf = params[e1+pos]
    else:
        lnf = 1.0

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
        
    # Adjust for wavelength calibration error
    newibinhi = ibinhi + delibinhi*deltaL
    newibinlo = ibinlo + delibinlo*deltaL
    
    # Bin and normalize spectrum
    if norm:
        normspec = af.NormSpec(modwave,fincident,snormtrunc,enormtrunc)
    else:
        normspec = fincident
        
    # Normalize if no radius was given
    if norad:
        normspec = normspec * totalflux/np.sum(normspec)
        
    # normspec is the final forward model spectrum
    binw = (newibinhi[1]-newibinhi[0])*(dataconv/databin)
    convmod = []
    for i in range(0,len(modindex)):
        convmod.append(af.ConvSpec(normspec[modindex[i][0]:modindex[i][1]],binw))
        
    convmodfl = [item for sublist in convmod for item in sublist]
    binmod = af.BinModel(convmodfl,newibinhi,newibinlo)
    
    iw = [i for i in range(0,len(binmod)) if ((i<polyindex or polyindex==-1) and binmod[i]!=0)]
    s2 = np.zeros(len(mastererr))
    for i in range(0,len(mastererr)): s2[i] = mastererr[i]*mastererr[i] + np.exp(lnf)
    likelihood = 0
    for i in iw:
        likelihood = likelihood + (masternorm[i]-binmod[i])**2/s2[i] + np.log(2.*np.pi*s2[i])
    likelihood = likelihood * -0.5
    
    if(np.isnan(likelihood) or np.isinf(likelihood)):
        print("Error: ")
        print(params)
        print(params1)
    return likelihood
    
# End of likelihood function

# Prior probability

def lnprior(x,teff):
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]

    priors = np.zeros(ndim)
    for i in range(0,ndim):
        if not bounds[i,0] <= x[i] <= bounds[i,1]:
            if not e1 <= i < e2:
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
        print('Prior Failed\n')
        return -np.inf
    
    grav = 0.
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        grav = 10**params[b1+pos]
    if 'Rad' in basic:
        pos = basic.index('Rad')
        radius = 6.371e8*params[b1+pos]
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        radius = 10**params[b1+pos]*dist*4.838e9*6.371e8 # convert R/D to Earth radii
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        radius = np.sqrt(params[b1+pos])*6.371e8
    else: radius = 11.2*6.371e8
    mass = grav*radius*radius/6.67e-8/1.898e30
    if mass<minmass or mass>maxmass:
        print('Mass out of Bound. Rad={0} log(g)={1} Mass={2}'.format(radius/6.371e8/11.2,np.log10(grav),mass))
        return -np.inf

    if len(gases)==0:
        mmw = 2.28
    else:
        mmw,rxsec = af.GetScaOpac(gases,params[g1:g2])
    scale = 1.38e-16*teff/mmw/grav;
    if scale/radius > 0.05:
        print('Gravity too low for reliable convergence. teff={0}, mu={1}, log(g)={2}, scale={3}, rad={4}, ratio={5}'.format(teff,mmw,np.log10(grav),scale/1.e5,radius/1.e5,scale/radius))
        return -np.inf
    
    penalty = 0
    if smooth:
        gamma = params[a2]
        for i in range(a1+1,a2-1):
            penalty = penalty - 0.5/gamma*(params[i+1] - 2*params[i] + params[i-1])**2
        penalty = penalty - 0.5*np.log(2*np.pi*gamma)*(a2-a1)
    # These lines weight the parameters based on the width of the prior if the boundaries cut off a significant amount of the normal distribution
    # However, it's not clear how important they are to comparing relative likelihoods
    #priors[1] = priors[1]/0.520
    #for i in range(9,14):
    #    priors[i] = priors[i]/0.9772

    return np.log(np.prod(priors))+penalty

# End of prior function

# Probability function

def lnprob(x,binshi,binslo,fluxrat,frathigh):
    params = plparams
    for i in range(0,len(nvars)):
        params[nvars[i]] = x[i]
        
    # Dummy variables in case they cannot be calculated.
    mass = 0.
    ctoo = 0.
    fetoh = 0.
    teff = 0.
    
    # Compute mass
    grav = 0.
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        grav = 10**params[b1+pos]
    if 'Rad' in basic:
        pos = basic.index('Rad')
        radius = 6.371e8*params[b1+pos]
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        radius = 10**params[b1+pos]*dist*4.838e9*6.371e8 # convert R/D to Earth radii
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        radius = np.sqrt(params[b1+pos])*6.371e8
    else: radius = 11.2*6.371e8
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
            carbon = carbon + cmult[i]*(10**params[g1+j-1]) # -1 because of hydrogen
    for i in range(0,len(ocompounds)):
        if ocompounds[i] in gases:
            j = gases.index(ocompounds[i])
            oxygen = oxygen + omult[i]*(10**params[g1+j-1])
    for i in range(0,len(zcompounds)):
        if zcompounds[i] in gases:
            j = gases.index(zcompounds[i])
            metals = metals + zmult[i]*(10**params[g1+j-1])

    ctoo = carbon/oxygen
    fetoh = np.log10(metals/0.0196)
    
    teff = planet.get_Teff()
    lp = lnprior(x,teff)
    # Check if any of the priors were out of bounds.
    if not np.isfinite(lp):
        return -np.inf, [mass,ctoo,fetoh,teff]

    # Check if an error returned a non-result.
    prob = lp + lnlike(x,binshi,binslo,fluxrat,frathigh)
    if np.isnan(prob):
        return -np.inf, [mass,ctoo,fetoh,teff]
    
    return prob, [mass,ctoo,fetoh,teff]

# End of probability function

# Set up the MCMC run
#print('Likelihood of input parameters: {0:f}'.format(lnlike(guess,ibinhi,ibinlo,binflux,binerr)))

# I have no idea why, but the first time GetModel() runs, it spits out a blackbody spectrum.
# This is "burn-in" step to avoid that and should not be commented out or removed.
testspectrum = GetModel(plparams)

if task=='Ensemble':

    print(eplist)
    print(len(eplist))
    
    allspec = np.zeros((len(eplist),len(modwave)))
    for ii in range(0,len(eplist)):
        print('Model #{0:d}'.format(ii))
        radfinal = 11.2

        pos = 0
        theta_planet = 0.
        if 'Rad' in basic:
            pos = basic.index('Rad')
            theta_planet = eplist[i][b1+pos]*6.371e8/dist/3.086e18
            radfinal = eplist[i][b1+pos]
        elif 'RtoD' in basic:
            pos = basic.index('RtoD')
            theta_planet = 10**eplist[i][b1+pos]
            radfinal = 10**eplist[i][b1+pos]*dist*3.086e18/6.371e8
        elif 'RtoD2U' in basic:
            pos = basic.index('RtoD2U')
            theta_planet = np.sqrt(eplist[i][b1+pos])*6.371e8/dist/3.086e18
            radfinal = np.sqrt(eplist[i][b1+pos])
        else:
            theta_planet = 11.2*6.371e8/dist/3.086e18
            # Default radius = Jupiter

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
        fout.write('{0:8.2f} {1:8.2f}'.format(modwave[i],modwave[i+1]))
        for j in range(0,len(eplist)):
            fout.write(' {0:8.5e}'.format(allspec[j][i]))
        fout.write('\n')
    sys.exit()

# End of ensemble-specific section

if task=='Retrieval':
    # Used to test the serial part of the code at the command line
    print('Test')
    print('Likelihood of input parameters: {0:f}'.format(lnlike(guess,ibinhi,ibinlo,binflux,binerr)))
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
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,backend=reader,pool=pool,args=(ibinhi,ibinlo,binflux,binerr))
    
    # Non-MPI sampler for testing purposes
    if not parallel:
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,backend=reader,args=(ibinhi,ibinlo,binflux,binerr))
    
    print('Walking...{0} {1} {2}...'.format(nwalkers,pllen,nsteps))

    if printfull: foutname = outdir + outfile + 'full.dat'
    else: foutname = outdir + outfile + 'dat'
    fchain = open(foutname,'w')

    if printfull: fchain.write('{0:d} {1:d} {2:d}'.format(nwalkers,nsteps,ndim+4))
    else: fchain.write('{0:d} {1:d} {2:d}'.format(nwalkers,(int)(nsteps/10),ndim+4))

    if atmtype=='Layers': fchain.write(' {0:d} {1:f} {2:f}\n'.format(a2-a1,minP-6.,maxP-6.))
    else: fchain.write('\n')

    for i in range(0,len(pnames)):
        if pnames[i]=='RtoD' or pnames[i]=='RtoD2U':
            fchain.write('Rad ')
        else:
            fchain.write('{0:s} '.format(pnames[i]))
    fchain.write('Mass C/O [Fe/H] Teff\n')

    # Actual MCMC run, and Write samples to output file

    maxlikli = 0.
    medianparams = guess
    
    i = 0
    for sample in sampler.sample(pos, iterations=nsteps):
        print('Sample: {0:d}\n'.format(i))
        coords = reader.get_last_sample().coords
        blobs = reader.get_blobs()[i]
        likli = reader.get_log_prob()[i]
        if max(likli) > maxlikli:
            j = np.argmax(likli)
            maxlikli = likli[j]
            medianparams = coords[j]
        print(reader.get_log_prob())
        if i%100==0 or printfull: print('Step number {0:d}'.format(i+1))
        if printfull or i>=nsteps*0.9:
            for i2 in range(0,len(coords)):
                for i3 in range(0,len(coords[i2])):
                    if pnames[i3]=='RtoD2U':
                        fchain.write('{0:f} '.format(np.sqrt(coords[i2][i3])/11.2))
                    elif pnames[i3]=='Rad' or pnames[i3]=='RtoD':
                        fchain.write('{0:f} '.format(coords[i2][i3]/11.2))
                    else:
                        fchain.write('{0:f} '.format(coords[i2][i3]))
                for i3 in range(0,len(blobs[i2])):
                    fchain.write('{0:f} '.format(blobs[i2][i3]))
                fchain.write('\n')
        i = i+1
    fchain.close()

    finalparams = medianparams
    
    # End of MCMC run

    if parallel: pool.close()

    xplot = np.linspace(1,nsteps,nsteps)
    first = int(0.9*len(sampler.chain[0]))
    if printfull: first=0

    print("Retrieval Complete")

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

    baselist = ['Rad','RtoD','RtoD2U','Log(g)','Cloud_Base','P_cl','Mass','C/O','[Fe/H]','Teff']
    bnamelist = ['Radius (R$_J$)','Radius (R$_J$)','Radius (R$_J$)','log(g)','P$_{cloud}$ (bar)','Base Pressure (bar)','Mass (M$_J$)','C/O','[Fe/H]','T$_{eff}$ (K)']
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
    if 'RtoD2U' in basic:
        rpos = basic.index('RtoD2U')
        bnames2[rpos] = 'Rad'
    if 'RtoD' in basic:
        rtemp = basic.index('RtoD')
        bnames2[rtemp] = 'Rad'
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
                bsamples3[i,j] = np.sqrt(bsamples2[i,j])/11.2
            elif basic[j]=='Rad' or basic[j]=='RtoD':
                bsamples3[i,j] = bsamples2[i,j]/11.2
            else:
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

    # Write parameter file of best fit model

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
    ffout.write('Data         {0:s} {1:5.1f} {2:5.1f}\n'.format(datain,databin,dataconv))
    ffout.write('Streams      {0:d}\n'.format(streams))
    ffout.write('Pressure     {0:5.1f} {1:5.1f}\n'.format(minP-6.,maxP-6.))
    ffout.write('Vres         {0:d}\n'.format(vres))
    if gray: ffout.write('Gray        {0:5.0f}\n'.format(tgray))
    ffout.write('Minimum_Mass {0:5.2f}\n'.format(minmass))
    ffout.write('Output       modelspectra    Short\n')
    ffout.write('Opacities    {0:s}\n'.format(opacdir))
    ffout.write('Tables       {0:s} {1:s}\n'.format(hires,lores))
    if not outmode=='': ffout.write('Output_Mode {0:s}\n'.format(outmode))
    
    ffout.write('Parameter    Initial    Mu    Sigma    Min    Max\n')
    if b1>=0:
        ffout.write('Basic\n')
        for i in range(b1,b2):
            if pnames[i]=='Rad' and 'RtoD2U' in basic:
                finalsigma[i] = np.sqrt((float)(finalparams[i])) - np.sqrt((float)(finalparams[i]))
                ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],np.sqrt((float)(finalparams[i])),np.sqrt((float)(finalparams[i])),finalsigma[i],np.sqrt(finalbounds[i,0]),np.sqrt(finalbounds[i,1])))
            else:
                ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if g1>=0:
        ffout.write('Gases     {0:s}\n'.format(gases[0]))
        for i in range(g1,g2):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if a1>=0:
        ffout.write('Atm       {0:s}\n'.format(atmtype))
        for i in range(a1,a2+1):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if c1>=0:
        ffout.write('Clouds    {0:d}    {1:s}\n'.format(cloudmod,hazestr))
        for i in range(c1,c2):
            ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[i]),(float)(finalparams[i]),finalsigma[i],finalbounds[i,0],finalbounds[i,1]))

    if e1>=0:
        ffout.write('End\ndeltaL	       0.      0.     1.0  -10.0    10.0\nlogf	       1.      1.     0.1    0.01  100.')

# End of retrieval-specific section
    
# Plot the best fit spectrum or spectrum requested at the command line
    
if task=='Spectrum':
    finalparams = plparams

radfinal = 11.2

while True:
    pos = 0
    theta_planet = 0.
    if 'Rad' in basic:
        pos = basic.index('Rad')
        theta_planet = finalparams[b1+pos]*6.371e8/dist/3.086e18
        radfinal = finalparams[b1+pos]
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        theta_planet = 10**finalparams[b1+pos]
        radfinal = 10**finalparams[b1+pos]*dist*3.086e18/6.371e8
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        theta_planet = np.sqrt(finalparams[b1+pos])*6.371e8/dist/3.086e18
        radfinal = np.sqrt(finalparams[b1+pos])
    else:
        theta_planet = 11.2*6.371e8/dist/3.086e18
        # Default radius = Jupiter
    
    spectrum = GetModel(finalparams)

    # Multiply by solid angle and collecting area
    fincident = np.zeros(len(spectrum))
    if mode<=1:
        for i in range(0,len(spectrum)):
            fincident[i] = spectrum[i] * theta_planet*theta_planet
            #if i==0: print "newtdepth: ",i,specwave[i],fincident[i] # print output for debugging purposes
            # erg/s/aperture/Hz
            # theta_planet is actually the radius/distance ratio
            # so its square converts flux at the surface to flux at the telescope
    if mode==2:
        fincident = spectrum
        
    if 'deltaL' in end:
        pos = end.index('deltaL')
        deltaL = finalparams[e1+pos]
    else:
        deltaL = 0.0
        
    # Adjust for wavelength calibration error
    newibinhi = ibinhi + delibinhi*deltaL
    newibinlo = ibinlo + delibinlo*deltaL        
    binw = (newibinhi[1]-newibinhi[0])*(dataconv/databin)
    convmod = af.ConvSpec(fincident,binw)
    binmod = af.BinModel(convmod,newibinhi,newibinlo)
    resid = binflux-binmod
    
    specout = 10000./binmid
    xmin = min(specout)
    xmax = max(specout)
    xmin = xmin - 0.05*(xmax-xmin)
    xmax = xmax + 0.05*(xmax-xmin)
    
    yref = max(max(binmod),max(binflux))
    ymin = -0.20 * yref
    ymax =  1.05 * yref

    # Plot the BINNED model/retrieved spectrum against the observations.
    
    fig4 = plt.figure(figsize=(10,7))
    ax = fig4.add_subplot(111)
    plt.axis((xmin,xmax,ymin,ymax))
    
    plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
    plt.ylabel('Flux (cgs)',fontsize=14)
    plt.tick_params(axis='both',which='major',labelsize=12)
    
    ax.errorbar(specout,binflux,binerr,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')
    ax.plot(specout,binmod,'-',linewidth=1,label='Retrieved Spectrum',c='b')
    ax.plot(specout,resid+ymin/2.,'-',linewidth=1,label='Residuals (offset)',c='r')
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')
    
    plt.legend(fontsize=12)
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
fig4name = 'plots' + outfile + 'binned.png'
fig4.savefig(fig4name)

# Create an output file of the BINNED model/retrieved spectrum.

foutnameb = 'modelspectra' + outfile + 'binned.dat'
fout = open(foutnameb,'w')
for i in range(0,len(specout)):
    fout.write('{0:8.2f} {1:8.2f} {2:8.5e} 0.0 0.0 {2:8.5e}\n'.format(binhi[i],binlo[i],binmod[i]))
fout.close()

# Plot the FULL-RES model/retrieved spectrum against the observations.
    
fig5name = 'plots' + outfile + 'fullres.png'

fig5 = plt.figure(figsize=(10,7))
ax = fig5.add_subplot(111)
plt.axis((xmin,xmax,ymin,ymax))

plt.xlabel('$\lambda$ ($\mu$m)',fontsize=14)
plt.ylabel('Flux (cgs)',fontsize=14)
plt.tick_params(axis='both',which='major',labelsize=12)

# Compute the residuals by binning to the observations without downsampling for resolving power.
wavemid = (wavehi+wavelo)/2.
binw = (wavehi[0]-wavehi[1])
convmod = af.ConvSpec(fincident,binw)
newbinshi,newbinslo = af.GetBins(modwave,wavehi,wavelo)
binmod = af.BinModel(convmod,newbinshi,newbinslo)

convflux2 = convflux[0]
converr2 = converr[0]
for i in range(1,len(convflux)):
    convflux2 = np.r_[convflux2,convflux[i]]
    converr2 = np.r_[converr2,converr[i]]
resid2 = convflux2-binmod

ax.errorbar(10000./wavemid,convflux2,converr2,capsize=3,marker='o',linestyle='',linewidth=1,label='Observations',c='k')
ax.plot(10000./modwave,fincident,'-',linewidth=1,label='Retrieved Spectrum',c='b')
ax.plot(10000./wavemid,resid2+ymin/2.,'-',linewidth=1,label='Residuals (convovled and offset)',c='r')
ax.plot([xmin,xmax],[0.,0.],'-',c='k')
ax.plot([xmin,xmax],[ymin/2.,ymin/2.],'--',c='k')

plt.legend(fontsize=12)
fig5.savefig(fig5name)

# Create an output file of the FULL-RES model/retrieved spectrum.

foutnamef = 'modelspectra' + outfile + 'fullres.dat'
fout = open(foutnamef,'w')
for i in range(0,len(modwave)-1):
    fout.write('{0:8.2f} {1:8.2f} {2:8.5e} 0.0 0.0 {2:8.5e}\n'.format(modwave[i],modwave[i+1],fincident[i]))

# End of plot spectrum

if task=='Retrieval': sys.exit()

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
    if (noisemode < 8 or noisemode == 11 or noisemode == 12) and (10000./modwave[0] > 5.01):
        print('Requested spectral mode does not match input wavelengths.')
        sys.exit()
    if (noisemode > 12 or noisemode == 9 or noisemode == 10) and (10000./modwave[-1] < 4.99):
        print('Requested spectral mode does not match input wavelengths.')
        sys.exit()

    print(starspec)
    calwave,flux_density,fnoise = AddNoise.addNoise(noisemode,mode,opacwave,spectrum,noise_params,starspec)
        
    calhi = np.zeros(len(calwave))
    callo = np.zeros(len(calwave))
    calhi[0] = calwave[0] + (calwave[0]-calwave[1])/2.
    if calhi[0] < 10000./modwave[0]: calhi[0] = 10000./modwave[0] + deltaL/1000.
    callo[-1] = calwave[-1] + (calwave[-2]-calwave[-1])/2.
    if callo[-1] > 10000./modwave[-1]: callo[-1] = 10000./modwave[-1] + deltaL/1000.
    
    for i in range(0,len(calwave)-1):
        calhi[i+1] = (calwave[i]+calwave[i+1])/2. + deltaL/1000.
        callo[i] = (calwave[i]+calwave[i+1])/2. + deltaL/1000.
        
    obsdepth = flux_density
    if mode<=1: noise = fnoise*obsdepth
    if mode==2: noise = fnoise
    obs_flux = obsdepth + np.random.normal(size=len(calwave))*noise

    foutname = 'modelspectra' + outfile + outmode + '.dat'
    ftest = open(foutname,'w')
    
    for i in range(0,len(calwave)-1):
        ftest.write('{0:8.2f} {1:8.2f} {2:8.5e} {3:8.5e} {3:8.5e} {4:8.5e}\n'.format(10000./(calhi[i]-deltaL/1000.),10000./(callo[i]-deltaL/1000.),obsdepth[i],noise[i],obs_flux[i]))

    # Plot the JWST mode against the observations.
        
    plotname = 'plots' + outfile + outmode  + '.fit.png'  # Need to put this in the plots folder.
    
    fig2 = plt.figure(figsize=(10,7))
    ax = fig2.add_subplot(111)

    snip = [i for i in range(0,len(modwave)) if calwave[0] < 10000./modwave[i] < calwave[-1]]
    
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
    
    ax.plot(10000./modwave[snip],fincident[snip],'-',linewidth=1,label='Model Spectrum',c='k')
    ax.errorbar(calwave,obsdepth,noise,capsize=3,marker='o',linestyle='',linewidth=1,label=outmode,c='b')
    # Alternate plot as a line without errorbars.
    #ax.plot(calwave,obsdepth,'-',linewidth=1,label=outmode,c='b')
    
    ax.plot([xmin,xmax],[0.,0.],'-',c='k')
    
    plt.legend(fontsize=12)
    plt.savefig(plotname) # Need to procedurally generate filename
