import sys
import numpy as np
import scipy.optimize as op
from scipy.interpolate import interp1d
from scipy.stats import invgamma
from distutils.util import strtobool
import emcee
import wrapPlanet_auto
import wrapPlanet_layer

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

def GetBinnedSpec(mode,mod_wave,mod_flux,obsmid):
    f = interp1d(mod_wave,mod_flux,kind='linear')
    int_flux = f(obsmid)
    if mode==0: answer = int_flux
    if mode==1: answer = int_flux
    if mode==2: answer = int_flux
    return answer

def GetNorm(specwave,fincident,lenspec,lennorm):
    normwave = np.zeros(len(lennorm)-1)
    normpoints = np.zeros(len(lennorm)-1)
    
    for i in range(0,len(normpoints)):
        first = int(lennorm[i])
        last = int(lennorm[i+1])
        normwave[i] = (specwave[first] + specwave[last-1]) / 2.
        normpoints[i] = np.mean(fincident[first:last])
    
    # Polynomial normalization
    fit = np.polyfit(normwave,normpoints,len(normwave))
    poly = np.poly1d(fit)
    return fincident[0:lenspec]/poly(specwave[0:lenspec])

def GetFluxPerBin(mode,mod_wave,mod_flux,obsmid):
    f = interp1d(mod_wave,mod_flux,kind='linear')
    int_flux = f(obsmid)
    bandwidth = np.fabs(np.gradient(obsmid)) * 2.99792458e10
    if mode==0: answer = int_flux*bandwidth
    if mode==1: answer = int_flux*bandwidth
    if mode==2: answer = int_flux
    return answer

def GetScaOpac(gases,abunds):
    h2 = 1. - np.sum(10**abunds)
    scaopac = 0.662e-27 * h2
    mmw = 2.28 * h2 # includes helium
    if 'h2o' in gases:
        scaopac = scaopac + 2.45e-27 * 10**(abunds[gases.index('h2o')-1])
        mmw = mmw + 18. * 10**(abunds[gases.index('h2o')-1])
    if 'ch4' in gases:
        scaopac = scaopac + 7.49e-27 * 10**(abunds[gases.index('ch4')-1])
        mmw = mmw + 16. * 10**(abunds[gases.index('ch4')-1])
    if 'co' in gases:
        scaopac = scaopac + 4.34e-27 * 10**(abunds[gases.index('co')-1])
        mmw = mmw + 28. * 10**(abunds[gases.index('co')-1])
    if 'co2' in gases:
        scaopac = scaopac + 8.24e-27 * 10**(abunds[gases.index('co2')-1])
        mmw = mmw + 44. * 10**(abunds[gases.index('co2')-1])
    if 'nh3' in gases:
        scaopac = scaopac + 5.37e-27 * 10**(abunds[gases.index('nh3')-1])
        mmw = mmw + 17. * 10**(abunds[gases.index('nh3')-1])
    if 'Burrows_alk' in gases:
        mmw = mmw + 24.1 * 10**(abunds[gases.index('Burrows_alk')-1])
    if 'Lupu_alk' in gases:
        mmw = mmw + 24.1 * 10**(abunds[gases.index('Lupu_alk')-1])
    if 'crh' in gases:
        mmw = mmw + 53.0 * 10**(abunds[gases.index('crh')-1])
    if 'feh' in gases:
        mmw = mmw + 56.8 * 10**(abunds[gases.index('feh')-1])
    if 'tio' in gases:
        mmw = mmw + 63.9 * 10**(abunds[gases.index('tio')-1])
    if 'vo' in gases:
        mmw = mmw + 66.9 * 10**(abunds[gases.index('vo')-1])
    return mmw,scaopac

def GetMollist(gases):
    mollist = np.zeros(len(gases))
    gaslist = ["h2","h2o","ch4","co","co2","nh3","h2s","Burrows_alk","Lupu_alk","crh","feh","tio","vo"]
    for i in range(0,len(gases)):
        if gases[i] in gaslist: mollist[i] = gaslist.index(gases[i])
        else: mollist[i] = 0
    return mollist

# Start of main program

# Read in planet parameters

settings = 'example.resolved.dat' # Bundled example file
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
pllen = pllen+1 # add 1 to include teff

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
datain = 'example.obs.dat' # Bundled example file
datatype = 'Observations'  # Bundled file is observations
wavei = 10000./0.60  # Full NIR wavelength range
wavef = 10000./4.99
degrade = 1          # Undegraded spectrum
nwalkers = 0         # Placeholder for later adjustment
nsteps = 30000       # Tested minimum required number of steps
minP = 0.0           # Pressure range to integrate over in cgs
maxP = 9.0
vres = 71            # number of layers for radiative transfer
streams = 1          # Use 1-stream by default
hazetype = 0         # No Clouds
hazestr = 'None'     # No Clouds
cloudmod = 0         # No Clouds
teff0 = 0.0          # Dummy variable for handling T_eff
natm = 0             # Placeholder in case T-P profile is omitted
prior = 'Uniform'    # Uniform priors
verbatim = False     # Interpolate the T-P profile
gray = False         # Used to create a gray atmosphere for testing
tgray = 1500         # Temperature of gray atmosphere
norad = False        # Flags if no radius variable is in the input file
output = '/nfs/turbo/lsa-arhowe/' # Default output on Flux

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
    elif line[0]=='Distance':
        if len(line)>1: dist = (float)(line[1])
    elif line[0]=='Data':
        if len(line)>1: datain = line[1]
        if len(line)>2: datatype = line[2]
    elif line[0]=='Degrade':
        if len(line)>1: degrade = (int)(line[1])
    elif line[0]=='N_Walkers':
        if len(line)>1: nwalkers = (int)(line[1])
        if override: nwalkers = 2
    elif line[0]=='N_Steps':
        if len(line)>1: nsteps = (int)(line[1])
        if override: nsteps = 10
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
# End read in settings

if override: parallel = False
if nwalkers==0: nwalkers = pllen*8           # Default number of walkers
if nwalkers<=2*pllen: nwalkers = pllen*2 + 2 # Minimum number of walkers
if nwalkers%2==1: nwalkers = nwalkers + 1    # Number of walkers must be even

# Set output file
output = output + name + '.' + modestr + '.' + str(pllen) + 'params' + str(nsteps/1000) + 'k.dat'

plparams = np.zeros(pllen)     # parameter list
mu       = np.zeros(pllen)     # Gaussian means
sigma    = np.zeros(pllen)     # Standard errors
bounds   = np.zeros((pllen,2)) # Bounds
guess    = np.zeros(pllen)

# Read in parameters
lines = fparams.readlines()

i=0
state = -1
basic = []
gases = ['h2']     # Default filler is hydrogen
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

plparams[-1] = teff0 # placeholder for Teff
guess[-1] = teff0
bounds[-1,0] = -1.e9
bounds[-1,1] = 1.e9

if gray: gases = []

b2 = b1+bnum
g2 = g1+gnum
a2 = a1+anum
c2 = c1+cnum
e2 = e1+enum
if smooth: a2 = a2-1
ilen = int(10 + c2-c1)
ngas = g2-g1+1

if 'RtoD2U' in basic:
    pos = basic.index('RtoD2U')
    plparams[b1+pos] = 10**plparams[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2
    guess[b1+pos] = 10**guess[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2    
    sigma[b1+pos] = 10**(sigma[b1+pos]-mu[b1+pos])
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
print "Reading in observations."
fobs = open(datain,'r')

obslines = fobs.readlines()
obslength = len(obslines)

obshi = np.zeros(obslength)
obslo = np.zeros(obslength)
fluxrat = np.zeros(obslength)
fratlow = np.zeros(obslength)
frathigh = np.zeros(obslength)

fitobs = 0

for i in range(0,obslength):
    obshi[i] = obslines[i].split()[0]
    obslo[i] = obslines[i].split()[1]
    fluxrat[i] = obslines[i].split()[5]
    fratlow[i] = obslines[i].split()[3]
    frathigh[i] = obslines[i].split()[4]
    if i>0 and fitobs==0 and obshi[i] > obshi[i-1]: fitobs = i

if fitobs==0: fitobs = obslength

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
if wavei > 10000/0.60: wavei = 10000./0.60
if wavef < 10000/4.99: wavef = 10000./4.99

# End of read in observations

# Set wavelength spectrum, full range, to be reduced to selected regions
speclen = (int)(21205./degrade)
specwave2 = np.zeros(speclen)
for i in range(0,speclen):
    specwave2[i] = 10000./(0.6*np.exp(i*degrade/10000.))

# Set up wavelength ranges and normalization
imin = np.where(specwave2<obshi[0])[0]-1
imax = np.where(specwave2<obslo[-1])[0]+1

starts = [imin[0]]
ends = []
startsnorm = []
endsnorm = []

obslennorm = [fitobs]

norm = False
for j in range(1,obslength):
    if norm:
        if obshi[j] < obslo[j-1]:
            js = np.where(specwave2<obshi[j])[0]-1
            je = np.where(specwave2<obslo[j-1])[0]+1
            startsnorm.append(js[0])
            endsnorm.append(je[0])
            obslennorm.append(j)
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
    ends.append(imax[0])
else:
    endsnorm.append(imax[0])
    obslennorm.append(len(obsmid))

specwave = []

for i in range(0,len(starts)):
    for j in range(starts[i],ends[i]):
        specwave.append(specwave2[j])
lenspec = len(specwave)
        
for i in range(0,len(startsnorm)):
    for j in range(startsnorm[i],endsnorm[i]):
        specwave.append(specwave2[j])

if norm:
    lennorm = np.zeros(len(startsnorm)+1)
    lennorm[0] = int(lenspec)
    for i in range(1,len(lennorm)):
        lennorm[i] = lennorm[i-1] + endsnorm[i-1] - startsnorm[i-1]

if norm:
    masternorm = GetNorm(obsmid,fluxrat,fitobs,obslennorm)
    fluxspecial = np.concatenate((frathigh[0:fitobs],fluxrat[fitobs:]),axis=0)
    mastererr = GetNorm(obsmid,fluxspecial,fitobs,obslennorm)
else:
    masternorm = fluxrat
    mastererr = frathigh
# End set up normalization

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

planet.MakePlanet(mode,specwave,cloudmod,hazetype,mollist)
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

    teff = planet.get_Teff()
    x[-1] = teff

    # multiply by solid angle and collecting area
    fincident = np.zeros(len(specflux))
    if mode<=1:
        for i in range(0,len(specflux)):
            fincident[i] = specflux[i] * theta_planet*theta_planet
            if i%999==0: print "newtdepth: ",i,specwave[i],fincident[i]
            # erg/s/aperture/Hz
            # theta_planet is actually the radius/distance ratio
            # so its square converts flux at the surface to flux at the telescope
    if mode==2:
        fincident = specflux

    # Adjust for wavelength calibration error
    obsmid2 = np.zeros(len(obsmid))
    for i in range(0,len(obsmid)):
        obsmid2[i] = 10000./(10000./obsmid[i] + 0.001*deltaL)
    obsflux = GetBinnedSpec(mode,specwave,fincident,obsmid2)

    # Bin and normalize spectrum
    if norm:
        normspec = GetNorm(specwave,fincident,lenspec,lennorm)
    else:
        normspec = fincident
    obsflux1 = GetBinnedSpec(mode,specwave[0:lenspec],normspec[0:lenspec],obsmid2[0:fitobs])

    obsflux = obsflux1

    s2 = mastererr**2 + np.exp(lnf)
    likelihood = -0.5 * np.sum( (masternorm[0:fitobs]-obsflux)**2/s2 + np.log(2.*np.pi*s2) )
    
    if(np.isnan(likelihood) or np.isinf(likelihood)):
        print "Error: ",params

    return likelihood

# End of likelihood function

# Used to test the serial part of the code at the command line
print lnlike(guess,binshi,binslo,fluxrat,frathigh)
print 'Likelihood of input parameters.'

eps = 0.01
for i in range(0,len(guess)):
    if guess[i] < bounds[i,0]+eps: guess[i] = bounds[i,0]+eps
    if guess[i] > bounds[i,1]-eps: guess[i] = bounds[i,1]-eps

# Prior probability
def lnprior(x):
    params = x
    priors = np.zeros(pllen)
    for i in range(0,pllen):
        if not bounds[i,0] <= params[i] <= bounds[i,1]:
            print 'Out of Bound {0:d} {1} {2} {3}\n'.format(i,params[i],bounds[i,0],bounds[i,1])
            return -np.inf
        if smooth and i==a2-1:
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
        gamma = params[a2-1]
        for i in range(a1+1,a2-2):
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
    if not np.isfinite(lp):
        return -np.inf
    prob = lp + lnlike(x,binshi,binslo,fluxrat,frathigh)
    if np.isnan(prob):
        return -np.inf
    return prob

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
pos = [guess + 0.1*eps*np.random.randn(ndim) for i in range(nwalkers)]

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
# An attempt at getting something useful out of an incomplete run.
'''
for result in sampler.sample(pos, iterations=3, storechain=False):
    position = result[0]
    fprog = open('chain.dat','a')
    for i2 in range(position.shape[0]):
        fprog.write('{0:4d}\n'.format(i2))
    fprog.close()
'''
posb,probb,state = sampler.run_mcmc(pos,nsteps)

if parallel: pool.close()

xplot = np.linspace(1,nsteps,nsteps)

# Alternate output code to output the full chain.
'''
foutname2 = output + '.full'
fout2 = open(foutname2,'w')
fout2.write("{0:d} {1:d} {2:d}\n".format(nwalkers,nsteps,ndim))
for i in range(0,len(sampler.chain)):
    for j in range(0,len(sampler.chain[0])):
        for k in range(0,len(sampler.chain[0,0])):
            fout2.write("{0:10.5f} ".format(sampler.chain[i,j,k]))
        fout2.write("\n")
'''

# Reduced output code to output only the last 10% of the chain.
foutname = output
fout = open(foutname,'w')
fout.write("{0:d} {1:d} {2:d}\n".format(nwalkers,nsteps,ndim))
first = int(0.9*len(sampler.chain[0]))
for i in range(0,len(sampler.chain)):
    for j in range(first,len(sampler.chain[0])):
        for k in range(0,len(sampler.chain[0,0])):
            fout.write("{0:10.5f} ".format(sampler.chain[i,j,k]))
        fout.write("\n")

print "Retrieval Complete"
