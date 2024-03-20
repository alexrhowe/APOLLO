#from enum import auto, StrEnum
from src.ApolloFunctions import (
    GetScaOpac,
    NormSpec,
    ConvSpec,
    BinModel
    )
import src.ApolloFunctions as af
from collections.abc import Callable
from distutils.util import strtobool
from functools import partial
import numpy as np
from scipy.interpolate import interp1d
from src.wrapPlanet import PyPlanet as PlanetCCore
import sys

from planet import (
    CPlanet_constructor,
    MakeModel_constructor,
    ObserveModel_constructor,
    Planet
    )

from user.P_points import P_profiles
from user.priors import priors, evaluate_default_priors
from user.TP_models import TP_models
from user.cloud_models import cloud_models

REarth_in_cm = 6.371e8
parsec_in_cm = 3.086e18
RJup_in_REarth = 11.2

#class ObservationMode(StrEnum):
#    RESOLVED = auto()
#    ECLIPSE = auto()
 #   TRANSIT = auto()

#class Task(StrEnum):
#    SPECTRUM = auto()
#    SPECTRAL_RANGE = auto()
#    ENSEMBLE = auto()
#    RETRIEVAL = auto()

default_command_line_arguments = dict(
    settings="examples/example.resolved.dat",
    task="Spectrum",
    name_suffix="",
    override=False,
    manual=False
    )

available_tasks = ["Spectrum", "Retrieval", "Ensemble", "Spectral_Range"]


def ProcessCLI(command_line_argv, default=default_command_line_arguments):
    num_clargs = len(command_line_argv)
    settings = command_line_argv[1] if num_clargs>1 else default["settings"]
    task = command_line_argv[2] if num_clargs>2 else default["task"]
    if task not in available_tasks:
        print('Error: specify "Spectrum" or "Retrieval", "Ensemble", or "Spectral_Range".')
        sys.exit()
    name_suffix = command_line_argv[3] + "." if num_clargs>3 else default["name_suffix"]
    override = "Serial" in command_line_argv
    manual = "Manual" in command_line_argv

    return dict(
        settings=settings,
        task=task,
        name_suffix=name_suffix,
        override=override,
        manual=manual
        )


# Everything that is read in line by line from APOLLO-style text files.
def ReadInputsfromFile(
        settings="examples/example.resolved.dat",
        task="Spectrum",
        name_suffix="",
        override=True,
        manual=True
        ):

    # Default Settings
    ######################################################################################
    name = 'example'     # Bundled example file
    mode = 0             # Emission spectrum
    modestr = 'Resolved' # Needed for output file name
    parallel = True      # Parallel operation
    datain = 'examples/example.obs.dat' # Bundled example file
    polyfit = False      # Switch to normalize the spectrum to a polynomial fit
    norm = False         # Dummy variable if polyfit is false
    checkpoint_file = None
    sampler = None
    samples_file = None
    num_samples = 0
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
    outdir = '/nfs/arthurad/samples'   # Default output directory
    short = False        # Switch to create short output file names
    printfull = False    # Switch to print the full sample array instead of the last 10%
    opacdir = '/nfs/arthurad/Opacities_0v10' # Default opacities directory
    minT = 75            # Minimum temperature for opacity tables
    maxT = 4000          # Maximum temperature for opacity tables

    norad = False        # Flags if no radius variable is in the input file
    natm = 0             # Placeholder in case T-P profile is omitted
    verbatim = False     # Interpolate the T-P profile
    tgray = 1500         # Temperature of gray atmosphere
    hazetype = 0         # No Clouds
    hazestr = 'None'     # No Clouds
    cloudmod = 0         # No Clouds

    hazelist = ['None','H2SO4','Polyacetylene','Tholin','Corundum','Enstatite','Forsterite','Iron','KCl','Na2S','NH3Ice','Soot','H2OCirrus','H2OIce','ZnS']
    ######################################################################################

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
        elif line[0]=='Sampler':
            if len(line)>1: sampler = line[1]
        elif line[0]=='Samples':
            if len(line)>1: samples_file = line[1]
            if len(line)>2: num_samples = (int)(line[2])
        elif line[0]=='Checkpoint':
            if len(line)>1: checkpoint_file = line[1]
        elif line[0]=='Convolve':
            if len(line)>1: dataconv = (int)(line[1])
        elif line[0]=='Binning':
            if len(line)>1: databin = (int)(line[1])
        elif line[0]=='Degrade':
            if len(line)>1: degrade = (int)(line[1])  # Only works with low-res data; mainly used to speed up execution for testing
        elif line[0]=='Prior':
            if len(line)>1: prior_type = line[1]
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
            else:
                short = False
                printfull = False
            if len(line)>3:
                if line[3]=='Short': short = True
                if line[3]=='Full': printfull = True
        elif line[0]=='Opacities':
            if len(line)>1: opacdir = line[1]
            
    # End read in settings

    # Read in model parameters

    print('Reading in parameters.')

    lines = fparams.readlines()

    plparams         = np.zeros(pllen)     # Parameter list
    mu               = np.zeros(pllen)     # Gaussian means
    sigma            = np.zeros(pllen)     # Standard errors
    bounds           = np.zeros((pllen,2)) # Bounds
    guess            = np.zeros(pllen)     # Used for initial conditions
    prior_types  = [prior_type] * pllen

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
    ilogg = -1

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
                hazetype = 0
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
            if len(line) >= 9:
                prior_functions[i] = line[8]

            if sigma[i] > 0.:
                pvars.append(plparams[i])
                nvars.append(i)
            
            # if guess[i]==0: guess[i] = 0.1*bounds[i,1]        # Prevents errors from occuring where parameters are zero.

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
            # ada: We want to impose a normal prior on log g,
            # while keeping uniform priors on everything else.
            elif line[0]=='Log(g)':
                ilogg = i
            if len(line)>6 and line[6]=='Ensemble':
                ensparams.append(i)
            i = i+1

    return dict(
        task=task,
        pnames=pnames,
        nvars=nvars,
        override=override,
        pllen=pllen,
        # checkpoint_file=checkpoint_file, # Not in ProcessInputs
        name=name,
        mode=mode,
        modestr=modestr,
        parallel=parallel,
        datain=datain,
        # polyfit=polyfit, # Not in ProcessInputs
        sampler=sampler,
        # samples_file=samples_file, # Not in ProcessInputs
        # num_samples=num_samples, # Not in ProcessInputs
        dataconv=dataconv,
        databin=databin,
        degrade=degrade,
        # prior_type=prior_type, # Not in ProcessInputs
        nwalkers=nwalkers,
        nsteps=nsteps,
        smooth=smooth,
        # igamma=igamma, # Not in ProcessInputs
        # ilogg=ilogg, # Not in ProcessInputs
        atmtype=atmtype,
        verbatim=verbatim,
        cloudmod=cloudmod,
        hazestr=hazestr,
        hazetype=hazetype,
        plparams=plparams,
        norad=norad,
        guess=guess,
        mu=mu,
        sigma=sigma,
        bounds=bounds,
        a1=a1,
        anum=anum,
        b1=b1,
        bnum=bnum,
        basic=basic,
        c1=c1,
        cnum=cnum,
        clouds=clouds,
        e1=e1,
        enum=enum,
        end=end,
        g1=g1,
        gnum=gnum,
        gases=gases,
        tstar=tstar,
        rstar=rstar,
        sma=sma,
        # starspec=starspec, # Not in ProcessInputs
        dist=dist,
        # RA=RA, # Not in ProcessInputs
        # dec=dec, # Not in ProcessInputs
        # minmass=minmass, # Not in ProcessInputs
        # maxmass=maxmass, # Not in ProcessInputs
        hires=hires,
        lores=lores,
        minP=minP,
        maxP=maxP,
        P_profile=P_profile,
        gray=gray,
        tgray=tgray,
        vres=vres,
        streams=streams,
        ensparams=ensparams,
        # outmode=outmode, # Not in ProcessInputs
        # exptime=exptime, # Not in ProcessInputs
        # outdir=outdir, # Not in ProcessInputs
        short=short,
        printfull=printfull,
        opacdir=opacdir
        )


def ProcessInputs(
        task,
        parallel,
        short,
        modestr,
        mode,
        name,
        pllen,
        pnames,
        vres,
        tstar,
        rstar,
        sma,
        nsteps,
        nvars,
        gray,
        tgray,
        atmtype,
        verbatim,
        sampler,
        nwalkers,
        override,
        printfull,
        b1,
        bnum,
        basic,
        g1,
        gnum,
        gases,
        a1,
        anum,
        c1,
        cnum,
        clouds,
        cloudmod,
        hazestr,
        hazetype,
        e1,
        enum,
        end,
        smooth,
        dist,
        plparams,
        norad,
        guess,
        mu,
        sigma,
        bounds,
        datain,
        dataconv,
        databin,
        degrade,
        lores,
        hires,
        opacdir,
        P_profile,
        minP,
        maxP,
        ensparams,
        streams
        ):
    cluster_mode = (task == 'Retrieval' and parallel)

    # Output file name: Object name, type of observation, # of parameters, and # of steps.
    if short:
        outfile = '/' + name + '.'
    else:
        outfile = '/' + name + '.' + modestr + '.' + str(pllen) + 'params' + str(int(nsteps/1000)) + 'k.'

    ndim = len(nvars)

    if gray:
        TP_model = TP_models["gray"]
        gases = []
    elif gases==[]:
        gases = ['h2he']

    if atmtype in TP_models:
        TP_model = TP_models[atmtype]

    if sampler == "emcee":
        if nwalkers==0: nwalkers = ndim*8            # Default number of walkers
        if nwalkers<2*ndim: nwalkers = ndim*2 + 2    # Minimum number of walkers
        if nwalkers%2==1: nwalkers = nwalkers + 1    # Number of walkers must be even
    elif sampler == "dynesty":
        nwalkers = 1

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

    # Special handling if area ratio is used instead of radius
    if 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        plparams[b1+pos] = 10**plparams[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2
        guess[b1+pos] = 10**guess[b1+pos] * dist**2 * 4.838e9**2 # convert (R/D)^2 to Earth radii^2    
        sigma[b1+pos] = guess[b1+pos]*(10**sigma[b1+pos]-1.)
        mu[b1+pos] = 10**mu[b1+pos] * dist**2 * 4.838e9**2
        sigma[b1+pos] = sigma[b1+pos]*mu[b1+pos]
        bounds[b1+pos] = 10**bounds[b1+pos] * dist**2 * 4.838e9**2

        radius_parameter = plparams[b1+pos]
        radius_case = "RtoD2U"

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

    wavelo = np.round(wavelo, 5)
    wavehi = np.round(wavehi, 5)
    wavemid = (wavehi+wavelo)/2.
    # End of read in observations

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

    # Specific to HD 106906 b. Mark regions of known high telluric contamination,
    # then use a separate error correction factor for those regions.
    telluric_regions = [
        (binmid_band<1.16112) |
        ((binmid_band>1.32801) & (binmid_band<1.49731)) |
        ((binmid_band>1.76505) & (binmid_band<2.07946)) |
        (binmid_band>2.41511)
        for binmid_band in binmid
        ]
    # print(telluric_regions)

    totalflux = 0
    for i in range(0,len(binflux)): totalflux = totalflux + binflux[i]*(binhi[i]-binlo[i])*1.e-4

    # Set statistical parameters
    #if 'logf' in end:
    #    pos = end.index('logf')
    #    plparams[e1+pos] = np.log(max(errhi**2))
    #    guess[e1+pos] = plparams[e1+pos]
    #    mu[e1+pos] = np.log(max(errhi**2))
    #    sigma[e1+pos] = abs(mu[e1+pos])/10.
    #    bounds[e1+pos,0] = np.log(min(errhi**2) * bounds[e1+pos,0])
    #    bounds[e1+pos,1] = np.log(max(errhi**2) * bounds[e1+pos,1])

    # Set the cross section tables if not already set.
    # Note that the "default" assumes a particular set of tables.
    minDL = 0
    maxDL = 0.
    if 'deltaL' in end:
        pos = end.index('deltaL')
        minDL = bounds[e1+pos,0]*0.001
        maxDL = bounds[e1+pos,1]*0.001

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
    elif istart[0]<0: istart[0]=0
    if len(iend)==0: iend = [len(opacwave)-1]
    elif iend[-1]>=len(opacwave): iend[-1] = len(opacwave)-1

    # Truncated and band-limited spectrum that does not extend beyond the range of the observations
    modwave = np.array(opacwave)[(int)(imin[0]):(int)(imax[0])]
    lenmod = len(modwave)
    '''
    # End set up model spectrum wavelength range

    # Handle bands and optional polynomial fitting
    bindex, modindex, modwave = af.SliceModel(bandlo,bandhi,opacwave,minDL,maxDL)

    # print(f"modwave is {modwave}")
    #if np.any(np.isnan(modwave)):
    #    print("There are at least some nans in the modwave.")

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
        
    # End of band handling

    # if task=='Spectrum' or task=='Ensemble': modwave = opacwave

    # Get indices of the edges of the observation bins in the model spectrum
    bins = af.GetBins(modwave,binlo,binhi)
    ibinlo = bins[0]
    ibinhi = bins[1]
    # print(f"ibinlo: {ibinlo}, ibinhi: {ibinhi}")

    # Needed to calculate the spectrum with the wavelength offset later.
    delmodwave = modwave + 0.001
    delbins = af.GetBins(delmodwave,binlo,binhi)
    delibinlo = delbins[0]-ibinlo
    delibinhi = delbins[1]-ibinhi
    # print(f"delibinlo: {delibinlo}, delibinhi: {delibinhi}")

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

    # Create Planet and read in opacity tables
    # planet = PlanetCCore()
    print('Haze type:',hazestr)
    print('Cloud model:',cloudmod)
    mode = int(mode)
    cloudmod = int(cloudmod)
    hazetype = int(hazetype)

    atmmod = 0
    if atmtype=='Layers' or atmtype in TP_models: atmmod = 0
    if atmtype=='Parametric': atmmod = 1

    if cloudmod==4:
        cloud_model = cloud_models["verbatim"]

    # switches = [mode,cloudmod,hazetype,streams,atmmod]

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

    MakePlanet_kwargs = CPlanet_constructor(
        observation_mode_index=mode,
        cloud_model_index=cloudmod,
        hazetype_index=hazetype,
        number_of_streams=streams,
        TP_model_switch=atmmod,
        model_wavelengths=modwave,
        Teff_calculation_wavelengths=modwavelo,
        gas_species_indices=mollist,
        gas_opacity_directory=opacdir,
        opacity_catalog_name=hires,
        Teff_opacity_catalog_name=lores
        )

    MakeModel_initialization_kwargs = MakeModel_constructor(
        model_wavelengths=modwave,
        gas_species=gases,
        minimum_model_pressure=minP,
        maximum_model_pressure=maxP,
        model_pressures=profin,
        number_of_atmospheric_layers=vres,
        stellar_effective_temperature=tstar,
        stellar_radius=rstar,
        semimajor_axis=sma,
        distance_to_system=dist,
        gas_start_index=g1,
        gas_end_index=g2,
        TP_model_type=atmtype,
        TP_model_function=TP_model,
        is_TP_profile_gray=gray,
        isothermal_temperature=tgray,
        is_TP_profile_verbatim=verbatim,
        TP_start_index=a1,
        TP_end_index=a2,
        number_of_TP_arguments=natm,
        number_of_params1_arguments=ilen,
        basic_parameter_names=basic,
        basic_start_index=b1,
        cloud_model_index=cloudmod,
        hazetype_index=hazetype,
        cloud_parameter_names=clouds,
        cloud_start_index=c1,
        calibration_parameter_names=end,
        calibration_start_index=e1,
        APOLLO_use_mode=task
        )

    ModelObservable_initialization_kwargs = ObserveModel_constructor(
        data_wavelengths=wavemid,
        model_wavelengths=modwave,
        model_indices_by_opacity_bins=modindex,
        distance_to_system=dist,
        bin_low_wavelength_indices_by_opacity_bins=ibinlo,
        bin_high_wavelength_indices_by_opacity_bins=ibinhi,
        bin_low_wavelength_offset_by_index=delibinlo,
        bin_high_wavelength_offset_by_index=delibinhi,
        binned_data_flux=totalflux,
        data_downbinning_factor=databin,
        data_convolution_factor=dataconv,
        data_band_bounding_indices=bandindex,
        observation_mode_index=mode,
        normalize_spectrum=polyfit,
        use_default_radius=norad,
        data_errors_binned_to_model=mastererr
        )
    '''
        modindex = model_indices_by_opacity_bins
        modwave = model_wavelengths
        modflux = emission_flux
        deltaL = wavelength_offset
        dist = distance_to_system
        ibinlo = bin_low_wavelength_indices_by_opacity_bins
        ibinhi = bin_high_wavelength_indices_by_opacity_bins
        delibinlo = bin_low_wavelength_offset_by_index
        delibinhi = bin_high_wavelength_offset_by_index
        totalflux = binned_data_flux
        databin = data_downbinning_factor
        dataconv = data_convolution_factor
        bandindex = data_band_bounding_indices
        mode = observation_mode_index
        norm = normalize_spectrum
        mastererr = data_errors_binned_to_model
    '''
    return dict(
        parameters=plparams,
        MakePlanet_kwargs=MakePlanet_kwargs,
        MakeModel_initialization_kwargs=MakeModel_initialization_kwargs,
        ModelObservable_initialization_kwargs=ModelObservable_initialization_kwargs,
        )


def MakePlanet(PlanetCCore_kwargs: CPlanet_constructor) -> Planet:
    return Planet("b", PlanetCCore_kwargs)



def MakeObservation(observation_constructor: ObserveModel_constructor):
    def ObservationModel(
            emission_flux:                               list,
            radius_parameter:                            float,
            radius_case:                                 str,
            wavelength_offset:                           float,
            # begin ModelObservable kwargs
            data_wavelengths:                            list,
            model_wavelengths:                           list,
            model_indices_by_opacity_bins:               list,
            distance_to_system:                          float,
            bin_low_wavelength_indices_by_opacity_bins:  list,
            bin_high_wavelength_indices_by_opacity_bins: list,
            bin_low_wavelength_offset_by_index:          list,
            bin_high_wavelength_offset_by_index:         list,
            binned_data_flux:                            list,
            data_downbinning_factor:                     float,
            data_convolution_factor:                     float,
            data_band_bounding_indices:                  list, #list[list[int]]
            observation_mode_index:                      int,
            normalize_spectrum:                          bool,
            use_default_radius:                          bool,
            data_errors_binned_to_model:                 list
            ):
        ##########################################
        # Start by relabeling all the arguments. #
        ##########################################

        modflux = emission_flux
        modindex = model_indices_by_opacity_bins
        wavemid = data_wavelengths
        modwave = model_wavelengths
        deltaL = wavelength_offset
        dist = distance_to_system
        ibinlo = bin_low_wavelength_indices_by_opacity_bins
        ibinhi = bin_high_wavelength_indices_by_opacity_bins
        delibinlo = bin_low_wavelength_offset_by_index
        delibinhi = bin_high_wavelength_offset_by_index
        totalflux = binned_data_flux
        databin = data_downbinning_factor
        dataconv = data_convolution_factor
        bandindex = data_band_bounding_indices
        mode = observation_mode_index
        norm = normalize_spectrum
        norad = use_default_radius
        mastererr = data_errors_binned_to_model
        
        ##########################################
        
        # params = plparams.copy()
        #for i in range(0,len(nvars)):
        #    params[nvars[i]] = x[i]
        
        # modflux, _ = GetModel(params)
        #sampled_derived_parameter_chain.write(
        #    "{0:f} {1:f} {2:f} {3:f}\n".format(*derived_parameters)
        #    )
        
        if radius_case == "Rad":
            theta_planet = (radius_parameter*REarth_in_cm)/(dist*parsec_in_cm)
        elif radius_case == "RtoD":
            theta_planet = 10**radius_parameter
        elif radius_case=="RtoD2U":
            theta_planet = np.sqrt(radius_parameter)*REarth_in_cm/(dist*parsec_in_cm)
        else:
            theta_planet = (RJup_in_REarth*REarth_in_cm)/(dist*parsec_in_cm)
        '''
        theta_planet = 0.
        if 'Rad' in basic:
            pos = basic.index('Rad')
            theta_planet = radius_parameter*REarth_in_cm/dist/parsec_in_cm
        elif 'RtoD' in basic:
            pos = basic.index('RtoD')
            theta_planet = 10**radius_parameter
        elif 'RtoD2U' in basic:
            pos = basic.index('RtoD2U')
            theta_planet = np.sqrt(radius_parameter)*REarth_in_cm/dist/parsec_in_cm
        else:
            theta_planet = RJup_in_REarth*REarth_in_cm/dist/parsec_in_cm
            # Default radius = Jupiter
        '''
        '''
        # Statistical parameters
        if 'deltaL' in end:
            pos = end.index('deltaL')
            deltaL = params[e1+pos]
        else:
            deltaL = 0.0
        '''
        '''
        if 'logf' in end:
            pos = end.index('logf')
            lnf = params[e1+pos]
        else:
            lnf = -100.0
        if 'logf_telluric' in end:
            pos = end.index('logf_telluric')
            lnf_telluric = params[e1+pos]
        else:
            lnf_telluric = -100.0
        '''
        # Multiply by solid angle and collecting area
        fincident = modflux if mode == 2 else modflux*(theta_planet**2)
        print(f"{len(modflux)=}")
        ''' # ADA: proposed pythonic rewrite
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
        '''
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
        newibinlo = ibinlo + delibinlo*deltaL
        newibinhi = ibinhi + delibinhi*deltaL
        
        # Bin and normalize spectrum
        snormtrunc = None # Arthur: added b/c snormtrunc is undefined
        enormtrunc = None # Arthur: added b/c enormtrunc is undefined
        normspec = NormSpec(modwave,fincident,snormtrunc,enormtrunc) if norm else fincident
        '''
        if norm:
            snormtrunc = None # Arthur: added b/c snormtrunc is undefined
            enormtrunc = None # Arthur: added b/c enormtrunc is undefined
            normspec = NormSpec(modwave,fincident,snormtrunc,enormtrunc)
        else:
            normspec = fincident
        '''
        # Normalize if no radius was given
        normspec = normspec * totalflux/np.sum(normspec) if norad else normspec
        '''
        if norad:
            normspec = normspec * totalflux/np.sum(normspec)
        '''

        # normspec is the final forward model spectrum
        binw = (newibinlo[1]-newibinlo[0])*(dataconv/databin)
        print(f"{binw=}")
        convmod = []
        for i in range(0,len(modindex)):
            print(f"{modindex[i]=}")
            convmod.append(ConvSpec(normspec[modindex[i][0]:modindex[i][1]],binw))
        convmod = [item for sublist in convmod for item in sublist]
        print(f"{len(convmod)=}")
        # convmod = af.ConvSpec(fincident,binw)
        # convmod = fincident
        binmod_list = []
        for i in range(0,len(modindex)):
            binmod_piece = BinModel(convmod,newibinlo[bandindex[i][0]:(bandindex[i][1]+1)],newibinhi[bandindex[i][0]:(bandindex[i][1]+1)])
            binmod_list.append(binmod_piece)
        
        binmod = [item for sublist in binmod_list for item in sublist]
        # s2 = mastererr**2 #+ np.exp(lnf) # Can add the logf part in a subsequent routine
        
        binned_observation = dict(
            wavelengths=wavemid,
            data=binmod
        )

        model_observation = dict(
            wavelengths=modwave,
            data=fincident
        )

        return binned_observation, model_observation

    return partial(ObservationModel, **observation_constructor._asdict())
