import sys
import numpy as np
import scipy.optimize as op
from scipy.interpolate import interp1d
from scipy.stats import invgamma
from distutils.util import strtobool
import emcee
from src import wrapPlanet_auto
from src import wrapPlanet_layer
import corner
import matplotlib.pyplot as plt

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

def GetBins(specwave,obshi,obslo):
    binshi = np.zeros(len(obshi))
    binslo = np.zeros(len(obshi))
    '''
    if(obshi[0]>specwave[0] or obslo[-1]<specwave[-1]):
        print "Wavelength out of range."
        return 0.
    '''
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

def GetNorm(specwave,fincident,startsnorm,endsnorm):
    normwave = np.zeros(len(startsnorm))
    normpoints = np.zeros(len(startsnorm))

    for i in range(0,len(normpoints)):
        stemp = (int)(startsnorm[i])
        etemp = (int)(endsnorm[i])
        normwave[i] = (specwave[stemp] + specwave[etemp-1])/2.
        normpoints[i] = np.mean(fincident[stemp:etemp])
    
    fit = np.polyfit(normwave,normpoints,len(normwave))
    poly = np.poly1d(fit)
    return fincident/poly(specwave)

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

override = False
if len(sys.argv)>2 and sys.argv[2]=='Serial': override = True

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
parallel = True      # Parallel operation
name = 'example'     # Bundled example file
tstar = 5770.        # Solar temperature
rstar = 1.0          # Solar radius
sma = 1.0            # Semi-Major Axis
dist = 10.0          # Standardized distance, 10 pc
RA = 0.0
dec = 0.0
datain = 'example.obs.dat' # Bundled example file
databin = 1          # Factor to bin down the observations
wavei = 10000./0.60  # Full NIR wavelength range
wavef = 10000./5.00
degrade = 1          # Undegraded spectrum
nwalkers = 0         # Placeholder for later adjustment
nsteps = 30000       # Tested minimum required number of steps
minP = 0.0           # Pressure range to integrate over in cgs
maxP = 9.0
vres = 71            # Number of layers for radiative transfer
streams = 1          # Use 1-stream by default
hazetype = 0         # No Clouds
hazestr = 'None'     # No Clouds
cloudmod = 0         # No Clouds
natm = 0             # Placeholder in case T-P profile is omitted
prior = 'Uniform'    # Uniform priors
verbatim = False     # Interpolate the T-P profile
gray = False         # Used to create a gray atmosphere for testing
tgray = 1500         # Temperature of gray atmosphere
norad = False        # Flags if no radius variable is in the input file
opacdir = '../Opacities' # Default opacities directory
outdir = 'samples'   # Default output directory
short = False        # Switch to create short output file names
polyfit = False      # Switch to normalize the spectrum to a polynomial fit
norm = False         # Dummy variable is polyfit is false
printfull = False    # Switch to print the full sample array instead of the last 10%
hires = 'nir'        # Default set of opacity tables to compute the spectra
lores = 'lores'      # Default low-resolution tables to compute Teff

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
    elif line[0]=='Parallel':
        if len(line)>1: parallel = strtobool(line[1])
    elif line[0]=='Object':
        if len(line)>1: name = line[1]
    elif line[0]=='Star':
        if len(line)>1: tstar = (float)(line[1])
        if len(line)>2: rstar = (float)(line[2])
        if len(line)>3: sma = (float)(line[3])
    elif line[0]=='Location':
        if len(line)>1: dist = (float)(line[1])
        if len(line)>2: ra = (float)(line[2])
        if len(line)>3: dec = (float)(line[3])
    elif line[0]=='Data':
        if len(line)>1: datain = line[1]
        if len(line)>2: databin = (int)(line[2])
        if len(line)>3 and line[3]=='Polyfit': polyfit = True
    elif line[0]=='Degrade':
        if len(line)>1: degrade = (int)(line[1])
    elif line[0]=='N_Walkers':
        if len(line)>1: nwalkers = (int)(line[1])
        if override: nwalkers = 2
    elif line[0]=='N_Steps':
        if len(line)>1: nsteps = (int)(line[1])
        if override: nsteps = 2
    elif line[0]=='Pressure':
        if len(line)>1: minP = (float)(line[1]) + 6.0 # Convert from bars to cgs
        if len(line)>2: maxP = (float)(line[2]) + 6.0
        if maxP <= minP: maxP = minP + 0.01
    elif line[0]=='Vres':
        if len(line)>1: vres = (int)(line[1])
    elif line[0]=='Streams':
        if len(line)>1: streams = (int)(line[1])
    elif line[0]=='Prior':
        if len(line)>1: prior = line[1]
    elif line[0]=='Gray':
        if len(line)>1: gray = strtobool(line[1])
        if len(line)>2: tgray = line[2]
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
    elif line[0]=='Tables':
        if len(line)>1: hires = line[1]
        if len(line)>2: lores = line[2]
        
# End read in settings

if override:
    parallel = False
    printfull = True
if nwalkers==0: nwalkers = pllen*8           # Default number of walkers
if nwalkers<2*pllen: nwalkers = pllen*2 + 2  # Minimum number of walkers
if nwalkers%2==1: nwalkers = nwalkers + 1    # Number of walkers must be even

# Output file name: Object name, type of observation, # of parameters, and # of steps.
if short:
    outfile = '/' + name + '.'
else:
    outfile = '/' + name + '.' + modestr + '.' + str(pllen) + 'params' + str(nsteps/1000) + 'k.'

plparams = np.zeros(pllen)     # parameter list
mu       = np.zeros(pllen)     # Gaussian means
sigma    = np.zeros(pllen)     # Standard errors
bounds   = np.zeros((pllen,2)) # Bounds
guess    = np.zeros(pllen)

# Read in parameters
print 'Reading in parameters.'
lines = fparams.readlines()

i=0
state = -1
pnames = []
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
        pnames.append(lines[j].split()[0])
        if state==0:
            basic.append(lines[j].split()[0])
            bnum = bnum+1
        if state==1:
            gases.append(lines[j].split()[0])
            gnum = gnum+1
        if state==2:
            atm.append(lines[j].split()[0])
            anum = anum+1
        if state==3:
            clouds.append(lines[j].split()[0])
            cnum = cnum+1
        if state==4:
            end.append(lines[j].split()[0])
            enum = enum+1
        plparams[i] = (float)(lines[j].split()[1])
        guess[i]    = (float)(lines[j].split()[1])
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
    guess[b1+pos] = 10**guess[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2    
    sigma[b1+pos] = guess[b1+pos]*(10**sigma[b1+pos]-1.)
    mu[b1+pos] = 10**mu[b1+pos] * dist**2 * 4.838e9**2
    sigma[b1+pos] = sigma[b1+pos]*mu[b1+pos]
    bounds[b1+pos] = 10**bounds[b1+pos] * dist**2 * 4.838e9**2

if 'Cloud_Base' in clouds:
    pos = clouds.index('Cloud_Base')
    bounds[c1+pos,0] = minP
    bounds[c1+pos,1] = maxP
if 'P_cl' in clouds:
    pos = clouds.index('P_cl')
    bounds[c1+pos,0] = minP
    bounds[c1+pos,1] = maxP
    
# Meant to make the temperatures uniform in log space
'''
if atmtype == 'Layers':
    print 'Layers\n'
    plparams[a1:a2] = np.log10(plparams[a1:a2])
    guess[a1:a2] = np.log10(guess[a1:a2])
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

# Read in observations

# Header contains info about star needed for JWST pipeline
print 'Reading in observations.'
fobs = open(datain,'r')

obslines = fobs.readlines()
obslength2 = len(obslines)

obshi2 = np.zeros(obslength2)
obslo2 = np.zeros(obslength2)
fluxrat2 = np.zeros(obslength2)
fratlow2 = np.zeros(obslength2)
frathigh2 = np.zeros(obslength2)

for i in range(0,obslength2):
    obshi2[i] = obslines[i].split()[0]
    obslo2[i] = obslines[i].split()[1]
    fluxrat2[i] = obslines[i].split()[5]
    fratlow2[i] = obslines[i].split()[3]
    frathigh2[i] = obslines[i].split()[4]

# Bin down the observations, if applicable
    
obslength = int(obslength2/databin)

obshi = np.zeros(obslength)
obslo = np.zeros(obslength)
fluxrat = np.zeros(obslength)
fratlow = np.zeros(obslength)
frathigh = np.zeros(obslength)

if databin>1:
    for i in range(0,obslength):
        obshi[i] = obshi2[i*databin]
        obslo[i] = obslo2[(i+1)*databin-1]
        fluxrat[i] = np.mean(fluxrat2[i*databin:(i+1)*databin])
        fratlow[i] = np.mean(fratlow2[i*databin:(i+1)*databin])/min(1,np.sqrt(databin-1))
        frathigh[i] = np.mean(frathigh2[i*databin:(i+1)*databin])/min(1,np.sqrt(databin-1))
else:
    obshi = obshi2
    obslo = obslo2
    fluxrat = fluxrat2
    fratlow = fratlow2
    frathigh = frathigh2
    
fitobs = 0
if polyfit:
    for i in range(0,len(obshi)):
        if i>0 and fitobs==0 and obshi[i] > obshi[i-1]: fitobs = i
if fitobs==0: fitobs = len(obshi)

totalflux = np.sum(fluxrat*(1./obslo - 1./obshi))

if 'logf' in end:
    pos = end.index('logf')
    plparams[e1+pos] = np.log(max(frathigh**2))
    mu[e1+pos] = np.log(max(frathigh**2))
    sigma[e1+pos] = abs(mu[e1+pos])/10.
    bounds[e1+pos,0] = np.log(min(frathigh**2) * bounds[e1+pos,0])
    bounds[e1+pos,1] = np.log(max(frathigh**2) * bounds[e1+pos,1])

obsmid = (obshi+obslo)/2.

# Set wavelength range of computed spectrum just wider than the input spectrum.
wavei = max(obshi) * 1.01
wavef = min(obslo) * 0.99

if wavei > 10000/5.0 and wavef > 10000/5.0:
    if wavei > 10000/0.6: wavei = 10000./0.6
    hires = 'nir'
elif wavei > 10000/5.0 and wavef < 10000/5.0:
    if wavei > 10000/0.6: wavei = 10000./0.6
    if wavef < 10000/30.0: wavef = 10000./30.0
    if 'Lupu_alk' in gases: gases.remove('Lupu_alk') # Temporary fix until I get Lupu_alk working.
    hires = 'wide'
elif wavei < 10000/5.0 and wavef < 10000/5.0:
    if wavef < 10000/30.0: wavef = 10000./30.0
    hires = 'mir'
    if 'Lupu_alk' in gases: gases.remove('Lupu_alk')
    
# End of read in observations

# Set wavelength spectrum, full range, to be reduced to selected regions

# Compute hires spectrum wavelengths
opacfile = opacdir + '/h2o.' + hires + '.dat'
fopac = open(opacfile,'r')
opacshape = fopac.readline().split()
fopac.close()
nwave = (int)(opacshape[6])
lmin = (float)(opacshape[7])
resolv = (float)(opacshape[9])

speclen = (int)(nwave/degrade)
specwave2 = np.zeros(speclen)
for i in range(0,speclen):
    specwave2[i] = 10000./(lmin*np.exp(i*degrade/resolv))
    
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
    
# Set up wavelength ranges and normalization
imin = np.where(specwave2<np.max(obshi))[0]-1
imax = np.where(specwave2<np.min(obslo))[0]+2
istart = np.where(specwave2<obshi[0])[0]-1
iend = np.where(specwave2<obslo[-1])[0]+2

if len(imin)==0: imin=[0]
if len(imax)==0: imax=[len(specwave2)-1]

# Truncated spectrum that does not extend beyond the range of the observations
specwave = specwave2[imin[0]:imax[0]]

starts = [istart[0]]
ends = []
startsnorm = []
endsnorm = []

# Optional polynomial fitting
if polyfit:
    # This part finds the indices at the boundaries of the bands in both the full spectrum and the observations
    norm = False
    slennorm = [fitobs]
    elennorm = []
    for j in range(1,obslength):
        if norm:
            if obshi[j] < obslo[j-1]:
                js = np.where(specwave2<obshi[j])[0]-1
                je = np.where(specwave2<obslo[j-1])[0]+1
                startsnorm.append(js[0])
                endsnorm.append(je[0])
                slennorm.append(j)
                elennorm.append(j)
                
        if not norm:
            if obshi[j] < obslo[j-1]:
                js = np.where(specwave2<obshi[j])[0]-1
                je = np.where(specwave2<obslo[j-1])[0]+1
                starts.append(js[0])
                ends.append(je[0])
                
            if obshi[j] > obshi[j-1]:
                js = np.where(specwave2<obshi[j])[0]-1
                je = np.where(specwave2<obslo[j-1])[0]+1
                ends.append(je[0])
                startsnorm.append(js[0])
                norm = True

    if not norm:
        ends.append(iend[0])
    else:
        endsnorm.append(iend[0])
        elennorm.append(len(obsmid))

    # This part finds the indices at the boundaries of the bands in the truncated spectrum
    wstart = (int)(np.log(10000./specwave2[imin[0]]/lmin)*resolv/degrade)
    strunc = np.zeros(len(starts))
    etrunc = np.zeros(len(ends))
    
    for i in range(0,len(starts)):
        strunc[i] = starts[i] - wstart
        etrunc[i] = ends[i] - wstart
    
    snormtrunc = np.zeros(len(startsnorm))
    enormtrunc = np.zeros(len(endsnorm))        
    
    for i in range(0,len(startsnorm)):
        snormtrunc[i] = startsnorm[i] - wstart
        enormtrunc[i] = endsnorm[i] - wstart

    # This part gets the master normalization of the observations
    if norm:
        masternorm = GetNorm(obsmid,fluxrat,slennorm,elennorm)
        fluxspecial = np.concatenate((frathigh[0:fitobs],fluxrat[fitobs:]),axis=0)
        mastererr = GetNorm(obsmid,fluxspecial,slennorm,elennorm)
    else:
        masternorm = fluxrat
        mastererr = frathigh
    
else:
    masternorm = fluxrat
    mastererr = frathigh
    lenspec = len(specwave)
# End polynomial fitting

bins = GetBins(specwave,obshi,obslo)

binshi = bins[0]
binslo = bins[1]

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

# Likelihood function
def lnlike(x,binshi,binslo,fluxrat,frathigh):
    params = x

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
        global norad
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

    # Statistical parameters
    if 'deltaL' in end:
        pos = end.index('deltaL')
        deltaL = plparams[e1+pos]
    else:
        deltaL = 0.0
    if 'logf' in end:
        pos = end.index('logf')
        lnf = plparams[e1+pos]
    else:
        lnf = 1.0
    
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

    # multiply by solid angle and collecting area
    fincident = np.zeros(len(specflux))
    if mode<=1:
        for i in range(0,len(specflux)):
            fincident[i] = specflux[i] * theta_planet*theta_planet
            #if i==0: print "newtdepth: ",i,specwave[i],fincident[i]
            # erg/s/aperture/Hz
            # theta_planet is actually the radius/distance ratio
            # so its square converts flux at the surface to flux at the telescope
    if mode==2:
        fincident = specflux
        
    # Adjust for wavelength calibration error
    obsmid2 = np.zeros(len(obsmid))
    for i in range(0,len(obsmid)):
        obsmid2[i] = 10000./(10000./obsmid[i] + 0.001*deltaL)
        
    # Bin and normalize spectrum
    if norm:
        normspec = GetNorm(specwave,fincident,snormtrunc,enormtrunc)
    else:
        normspec = fincident
        
    #obsflux = GetBinnedSpec(mode,specwave,normspec,obsmid2[0:fitobs])
    obsflux = GetBinnedSpec(mode,specwave,normspec,binshi,binslo)
    iw = [i for i in range(0,len(obsflux)) if (i<fitobs and obsflux[i]!=0)]
    
    s2 = mastererr**2 + np.exp(lnf)
    likelihood = -0.5 * np.sum( (masternorm[iw]-obsflux[iw])**2/s2[iw] + np.log(2.*np.pi*s2[iw]) )
    
    if(np.isnan(likelihood) or np.isinf(likelihood)):
        print "Error: ",params
    
    return likelihood

# End of likelihood function

# Used to test the serial part of the code at the command line
print 'Likelihood of input parameters: {0:f}'.format(lnlike(guess,binshi,binslo,fluxrat,frathigh))
#sys.exit()

eps = 0.01
for i in range(0,len(guess)):
    if guess[i] < bounds[i,0]+eps*(mu[i]-bounds[i,0]): guess[i] = bounds[i,0]+eps*(mu[i]-bounds[i,0])
    if guess[i] > bounds[i,1]-eps*(bounds[i,1]-mu[i]): guess[i] = bounds[i,1]-eps*(bounds[i,1]-mu[i])
    
# Prior probability
def lnprior(x):
    params = x
    priors = np.zeros(pllen)
    for i in range(0,pllen):
        if not bounds[i,0] <= params[i] <= bounds[i,1]:
            if not e1 <= i < e2:
                print 'Out of Bound {0:d} {1} {2} {3}'.format(i,params[i],bounds[i,0],bounds[i,1])
            return -np.inf
        if smooth and i==a2:
            priors[i] = invgamma.pdf(params[i],1,scale=5.e-5) # gamma with inverse gamma function prior, alpha=1, beta=5.e-5
        else:
            if prior == 'Uniform':
                priors[i] = 1/(bounds[i,1]-bounds[i,0])
            if prior == 'Normal':
                priors[i] = 1/sigma[i]/2.5066 * np.exp(-(params[i]-mu[i])*(params[i]-mu[i])/2/sigma[i]/sigma[i])
    abundsum = np.sum(10**params[g1:g2])
    if abundsum>1.0:
        print 'Prior Failed\n'
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
    if mass<1. or mass>80.:
        print 'Mass out of Bound {0} {1} {2}\n'.format(radius,grav,mass)
        return -np.inf
    
    penalty = 0
    if smooth:
        gamma = params[a2]
        for i in range(a1+1,a2-1):
            penalty = penalty - 0.5/gamma*(params[i+1] - 2*params[i] + params[i-1])**2
        penalty = penalty - 0.5*np.log(2*np.pi*gamma)

    # These lines weight the parameters based on the width of the prior if the boundaries cut off a significant amount of the normal distribution
    # However, it's not clear how important they are to comparing relative likelihoods
    #priors[1] = priors[1]/0.520
    #for i in range(9,14):
    #    priors[i] = priors[i]/0.9772

    return np.log(np.prod(priors))+penalty

print 'Prior probability of input parameters: ', lnprior(guess)

# Probability function
def lnprob(x,binshi,binslo,fluxrat,frathigh):
    lp = lnprior(x)
    params = x

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
            carbon = carbon + cmult[i]*10**params[g1+j-1] # -1 because of hydrogen
    for i in range(0,len(ocompounds)):
        if ocompounds[i] in gases:
            j = gases.index(ocompounds[i])
            oxygen = oxygen + omult[i]*10**params[g1+j-1]
    for i in range(0,len(zcompounds)):
        if zcompounds[i] in gases:
            j = gases.index(zcompounds[i])
            metals = metals + zmult[i]*10**params[g1+j-1]

    ctoo = carbon/oxygen
    fetoh = np.log10(metals/0.0196)
    
    teff = planet.get_Teff()

    # Check if any of the priors were out of bounds.
    if not np.isfinite(lp):
        return -np.inf, [mass,ctoo,fetoh,teff]

    # Check if an error returned a non-result.
    prob = lp + lnlike(x,binshi,binslo,fluxrat,frathigh)
    if np.isnan(prob):
        return -np.inf, [mass,ctoo,fetoh,teff]
    
    return prob, [mass,ctoo,fetoh,teff]

pos1 = guess + 0.1*eps*np.random.randn(pllen)
testprior = lnprior(pos1)

if parallel:
    # MPI Setup
    from emcee.utils import MPIPool
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

# Walking
ndim = pllen
pos = [guess + 0.1*eps*guess*np.random.randn(ndim) for i in range(nwalkers)]

#fprog = open('chain.dat','w')
#fprog.close()
#result = pymultinest.solve(LogLikelihood=lnlike, Prior=lnprior, n_dims=ndim, outputfiles_basename=pmntest, verbose=True)

# MPI sampler
if parallel:
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool,args=(binshi,binslo,fluxrat,frathigh))
    
# Non-MPI sampler for testing purposes
if not parallel:
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(binshi,binslo,fluxrat,frathigh))

print 'Walking...{0} {1} {2}...'.format(nwalkers,pllen,nsteps)

rpos = -1
if 'RtoD2U' in pnames:
    rpos = pnames.index('RtoD2U')
    pnames[rpos] = 'Rad'
if 'RtoD' in pnames:
    rtemp = pnames.index('RtoD')
    pnames[rtemp] = 'Rad'

# Start samples output file
if printfull: foutname = outdir + outfile + 'full.dat'
else: foutname = outdir + outfile + 'dat'
fchain = open(foutname,'w')

if printfull: fchain.write('{0:d} {1:d} {2:d}'.format(nwalkers,nsteps,ndim+4))
else: fchain.write('{0:d} {1:d} {2:d}'.format(nwalkers,nsteps/10,ndim+4))
if atmtype=='Layers': fchain.write(' {0:d} {1:f} {2:f}\n'.format(a2-a1,minP-6.,maxP-6.))
else: fchain.write('\n')

for i in range(0,len(pnames)):
    if i==rpos:
        fchain.write('Rad ')
    else:
        fchain.write('{0:s} '.format(pnames[i]))
fchain.write('Mass C/O [Fe/H] Teff\n')

i = 0
for sample in sampler.sample(pos, iterations=nsteps):
    if i%100==0 or printfull: print 'Step number {0:d}'.format(i+1)
    if printfull or i>=nsteps*0.9:
        for i2 in range(0,len(sample[0])):
            for i3 in range(0,len(sample[0][i2])):
                if i3==rpos:
                    fchain.write('{0:f} '.format(np.sqrt(sample[0][i2][i3])/11.2))
                elif pnames[i3]=='Rad' or pnames[i3]=='RtoD':
                    fchain.write('{0:f} '.format(sample[0][i2][i3]/11.2))                    
                else:
                    fchain.write('{0:f} '.format(sample[0][i2][i3]))
            for i3 in range(0,len(sample[3][i2])):
                fchain.write('{0:f} '.format(sample[3][i2][i3]))   # sample[3] is the blob
            fchain.write('\n')
    i = i+1
#fchain.close()

if parallel: pool.close()

xplot = np.linspace(1,nsteps,nsteps)
first = int(0.9*len(sampler.chain[0]))

print "Retrieval Complete"

# Create waterfall plots of results

gsamples = sampler.chain[:,first:,g1:g2]
bsamples = sampler.chain[:,first:,b1:b2]
tsamples = sampler.chain[:,first:,a1:a2]
blobs = sampler.blobs

gsamples2 = gsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),g2-g1),order='F')
bsamples2 = bsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),b2-b1),order='F')
tsamples2 = tsamples.reshape(((len(sampler.chain[0])-first)*len(sampler.chain),a2-a1),order='F')
        
basic2 = basic
basic2.append('Mass')
basic2.append('C/O')
basic2.append('[Fe/H]')
basic2.append('Teff')

baselist = ['Rad','RtoD','RtoD2U','Log(g)','Cloud_Base','P_cl','Mass','C/O','[Fe/H]','Teff']
bnamelist = ['Radius (R$_J$)','Radius (R$_J$)','Radius (R$_J$)','log(g)','Base Pressure (bar)','Base Pressure (bar)','Mass (M$_J$)','C/O','[Fe/H]','T$_{eff}$ (K)']
gaslist = ['h2','h2only','he','h-','h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
gnamelist = ['H$_2$+He','H$_2$','He','[H-]','[H$_2$O]','[CH$_4$]','[CO]','[CO$_2$]','[NH$_3$]','[H$_2$S]','[Na,K]','[Na,K]','[CrH]','[FeH]','[TiO]','[VO]','[HCN]','[N2]','[PH3]']
bnames = []
gnames = []

for i in range(0,len(basic2)):
    if basic2[i] in baselist:
        j = baselist.index(basic2[i])
        bnames.append(bnamelist[j])

for i in range(0,len(pnames)):
    if pnames[i] in gaslist:
        j = gaslist.index(pnames[i])
        gnames.append(gnamelist[j])

rpos = -1
if 'RtoD2U' in basic:
    rpos = basic.index('RtoD2U')
    basic2[rpos] = 'Rad'
if 'RtoD' in basic:
    rtemp = basic.index('RtoD')
    basic2[rtemp] = 'Rad'
    
bsamples3 = np.zeros((len(bsamples2),b2-b1+4))
lenbasic = len(bsamples2[0])

for i in range(0,len(bsamples3)):
    for j in range(0,lenbasic):
        if j==rpos:
            bsamples3[i,j] = np.sqrt(bsamples2[i,j])/11.2
        elif basic[j]=='Rad' or basic[j]=='RtoD':
            bsamples3[i,j] = bsamples2[i,j]/11.2
        else:
            bsamples3[i,j] = bsamples2[i,j]
    for j in range(0,len(blobs[0][0])):
        bsamples3[i,lenbasic+j] = blobs[0][i][j]

#bsamples4 = [x for x in bsamples3 if not np.isinf(x[-1])]
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

plist = np.zeros(anum)
for i in range(0,anum): plist[i] = 10**(maxP + (minP-maxP)*i/(anum-1)) * 1.e-6
#tlist = np.percentile(sampler.chain[:,first:,a1:a2],[16,50,84],axis=0)
tlist = np.percentile(tsamples3,[16,50,84],axis=0)

# Plot the T-P profile

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

finalparams = np.percentile(sampler.chain[:,first:,:],[50,84],axis=0)[:,0]

finalsigma = np.zeros(pllen)
for i in range(0,pllen):
    finalsigma[i] = (float)(finalparams[1][i]) - (float)(finalparams[0][i])
outparams = '.' + outfile + 'retrieved.dat'
ffout = open(outparams,'w')

ffout.write('Mode        {0:s}\nObject      {1:s}\nStar        {2:5.0f} {3:5.2f} {4:8.3f}\nLocation    {5:6.2f} {6:6.2f} {7:6.2f}\nData        {8:s}    1\nDegrade     1\nN_Steps     2\nStreams     1\nPressure    {9:5.1f} {10:5.1f}\nPrior       Uniform\nOutput      modelspectra    Short\nOpacities {11:s}\nParameter    Initial    Mu    Sigma    Min    Max\n'.format(modestr,name,tstar,rstar,sma,dist,ra,dec,datain,minP-6.,maxP-6.,opacdir))
ffout.write('Basic\n')
for i in range(b1,b2):
    if pnames[i]=='Rad' and 'RtoD2U' in basic:
        finalsigma[i] = np.sqrt((float)(finalparams[1][i])) - np.sqrt((float)(finalparams[0][i]))
        ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],np.sqrt((float)(finalparams[0][i])),np.sqrt((float)(finalparams[0][i])),finalsigma[i],np.sqrt(bounds[i,0]),np.sqrt(bounds[i,1])))
    else:
        ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[0][i]),(float)(finalparams[0][i]),finalsigma[i],bounds[i,0],bounds[i,1]))
ffout.write('Gases\n')
for i in range(g1,g2):
    ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[0][i]),(float)(finalparams[0][i]),finalsigma[i],bounds[i,0],bounds[i,1]))
ffout.write('Atm       {0:s}\n'.format(atmtype))
for i in range(a1,a2+1):
    ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[0][i]),(float)(finalparams[0][i]),finalsigma[i],bounds[i,0],bounds[i,1]))
ffout.write('Clouds    {0:d}    {1:s}\n'.format(cloudmod,hazestr))
for i in range(c1,c2):
    ffout.write('{0:s}    {1:8.2f}    {2:8.2f}    {3:8.2f}    {4:8.2f}    {5:8.2f}\n'.format(pnames[i],(float)(finalparams[0][i]),(float)(finalparams[0][i]),finalsigma[i],bounds[i,0],bounds[i,1]))
ffout.write('End\ndeltaL	       0.      0.     1.0  -10.0    10.0\nlogf	       1.      1.     0.1    0.01  100.')
