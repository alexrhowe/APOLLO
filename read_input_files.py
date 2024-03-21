from distutils.util import strtobool
import numpy as np
import sys

from user.P_points import P_profiles


# Everything that is read in line by line from APOLLO-style text files.
def ReadInputsfromFile(
    settings="examples/example.resolved.dat",
    task="Spectrum",
    name_suffix="",
    override=True,
    manual=True,
):

    # Default Settings
    ######################################################################################
    name = "example"  # Bundled example file
    mode = 0  # Emission spectrum
    modestr = "Resolved"  # Needed for output file name
    parallel = True  # Parallel operation
    datain = "examples/example.obs.dat"  # Bundled example file
    polyfit = False  # Switch to normalize the spectrum to a polynomial fit
    norm = False  # Dummy variable if polyfit is false
    checkpoint_file = None
    sampler = None
    samples_file = None
    num_samples = 0
    dataconv = 1  # Factor to convolve the observations to blur the spectrum or account for actual resolving power
    databin = 1  # Factor to bin down the observations for simpler fitting
    degrade = 1  # Factor to degrade the model spectrum for faster calculation
    prior = "Uniform"  # Uniform priors
    nwalkers = 0  # Placeholder for later adjustment
    nsteps = 30000  # Tested minimum required number of steps
    tstar = 5770.0  # Solar temperature
    rstar = 1.0  # Solar radius
    sma = 1.0  # Semi-Major Axis
    starspec = ""  # Optional stellar spectrum file
    dist = 10.0  # Standardized distance, 10 pc
    RA = 0.0  # Right ascension
    dec = 0.0  # Declination
    minmass = 0.5  # Minimum mass in Jupiter masses
    maxmass = 80.0  # Maximum mass in Jupiter masses
    hires = ""  # Default set of opacity tables to compute the spectra
    lores = "lores"  # Default low-resolution tables to compute Teff
    minP = 0.0  # Pressure range to integrate over in cgs, default 1 mubar to 1 kbar.
    maxP = 9.0
    gray = False  # Used to create a gray atmosphere for testing
    vres = 71  # Number of layers for radiative transfer
    streams = 1  # Use 1-stream by default
    wavei = 10000.0 / 0.60  # Full NIR wavelength range
    wavef = 10000.0 / 5.00
    outmode = ""  # JWST observing mode
    exptime = 1000.0  # Exposure time in seconds
    outdir = "/nfs/arthurad/samples"  # Default output directory
    short = False  # Switch to create short output file names
    printfull = False  # Switch to print the full sample array instead of the last 10%
    opacdir = "/nfs/arthurad/Opacities_0v10"  # Default opacities directory
    minT = 75  # Minimum temperature for opacity tables
    maxT = 4000  # Maximum temperature for opacity tables

    norad = False  # Flags if no radius variable is in the input file
    natm = 0  # Placeholder in case T-P profile is omitted
    verbatim = False  # Interpolate the T-P profile
    tgray = 1500  # Temperature of gray atmosphere
    hazetype = 0  # No Clouds
    hazestr = "None"  # No Clouds
    cloudmod = 0  # No Clouds

    hazelist = [
        "None",
        "H2SO4",
        "Polyacetylene",
        "Tholin",
        "Corundum",
        "Enstatite",
        "Forsterite",
        "Iron",
        "KCl",
        "Na2S",
        "NH3Ice",
        "Soot",
        "H2OCirrus",
        "H2OIce",
        "ZnS",
    ]
    ######################################################################################

    try:
        fparams = open(settings, "r")
    except:
        print("\nError: input file not found.\n")
        sys.exit()

    lines1 = fparams.readlines()
    pllen = -1
    for i in range(0, len(lines1)):
        if len(lines1[i].split()) >= 6:
            pllen = pllen + 1

    fparams.close()
    fparams = open(settings, "r")

    # Read in settings

    nlines = 0
    while True:
        last_pos = fparams.tell()
        line = fparams.readline().split()
        if len(line) >= 6:  # ends the loop when the parameters start
            if line[0] == "Parameter":
                break
            else:
                fparams.seek(last_pos)
                break

        nlines = nlines + 1
        if nlines > 100:
            break  # prevents getting stuck in an infinite loop

        elif line[0] == "Object":
            if len(line) > 1:
                name = line[1]
        if line[0] == "Mode":
            if len(line) > 1:
                modestr = line[1]
            if modestr == "Resolved":
                mode = 0
            if modestr == "Eclipse":
                mode = 1
            if modestr == "Transit":
                mode = 2
        elif line[0] == "Parallel":
            if len(line) > 1:
                parallel = strtobool(line[1])
        elif line[0] == "Data":
            if len(line) > 1:
                datain = line[1]
            if len(line) > 2:
                if line[2] == "Polyfit":
                    polyfit = True
        elif line[0] == "Sampler":
            if len(line) > 1:
                sampler = line[1]
        elif line[0] == "Samples":
            if len(line) > 1:
                samples_file = line[1]
            if len(line) > 2:
                num_samples = (int)(line[2])
        elif line[0] == "Checkpoint":
            if len(line) > 1:
                checkpoint_file = line[1]
        elif line[0] == "Convolve":
            if len(line) > 1:
                dataconv = (int)(line[1])
        elif line[0] == "Binning":
            if len(line) > 1:
                databin = (int)(line[1])
        elif line[0] == "Degrade":
            if len(line) > 1:
                degrade = (int)(
                    line[1]
                )  # Only works with low-res data; mainly used to speed up execution for testing
        elif line[0] == "Prior":
            if len(line) > 1:
                prior_type = line[1]
        elif line[0] == "N_Walkers":
            if len(line) > 1:
                nwalkers = (int)(line[1])
            if override:
                nwalkers = 2
        elif line[0] == "N_Steps":
            if len(line) > 1:
                nsteps = (int)(line[1])
        elif line[0] == "Star":
            if len(line) > 1:
                tstar = (float)(line[1])
            if len(line) > 2:
                rstar = (float)(line[2])
            if len(line) > 3:
                sma = (float)(line[3])
        elif line[0] == "Star_Spec":
            if len(line) > 1:
                starspec = line[1]
        elif line[0] == "Location":
            if len(line) > 1:
                dist = (float)(line[1])
            if len(line) > 2:
                RA = (float)(line[2])
            if len(line) > 3:
                dec = (float)(line[3])
        elif line[0] == "Mass_Limits":
            if len(line) > 1:
                minmass = (float)(line[1])
            if len(line) > 2:
                maxmass = (float)(line[2])
        elif line[0] == "Tables":
            if len(line) > 1:
                hires = line[1]
            if len(line) > 2:
                lores = line[2]
        elif line[0] == "Pressure":
            if len(line) > 1:
                minP = (float)(line[1]) + 6.0  # Convert from bars to cgs
            if len(line) > 2:
                maxP = (float)(line[2]) + 6.0
            if maxP <= minP:
                maxP = minP + 0.01
            if len(line) > 3:
                P_profile = P_profiles[line[3]] + 6.0
            else:
                P_profile = None
        elif line[0] == "Gray":
            if len(line) > 1:
                gray = strtobool(line[1])
            if len(line) > 2:
                tgray = line[2]
        elif line[0] == "Vres":
            if len(line) > 1:
                vres = (int)(line[1])
        elif line[0] == "Streams":
            if len(line) > 1:
                streams = (int)(line[1])
        elif line[0] == "Output_Mode":
            if len(line) > 1:
                outmode = line[1]
            if len(line) > 2:
                exptime = (float)(line[2])
        elif line[0] == "Output":
            if len(line) > 1:
                outdir = line[1]
            if len(line) > 2:
                if line[2] == "Short":
                    short = True
                if line[2] == "Full":
                    printfull = True
            else:
                short = False
                printfull = False
            if len(line) > 3:
                if line[3] == "Short":
                    short = True
                if line[3] == "Full":
                    printfull = True
        elif line[0] == "Opacities":
            if len(line) > 1:
                opacdir = line[1]

    # End read in settings

    # Read in model parameters

    print("Reading in parameters.")

    lines = fparams.readlines()

    plparams = np.zeros(pllen)  # Parameter list
    mu = np.zeros(pllen)  # Gaussian means
    sigma = np.zeros(pllen)  # Standard errors
    bounds = np.zeros((pllen, 2))  # Bounds
    guess = np.zeros(pllen)  # Used for initial conditions
    prior_types = [prior_type] * pllen

    i = 0
    state = -1
    pnames = []
    pvars = []
    nvars = []
    basic = []
    gases = []
    atm = []
    clouds = []
    end = []
    atmtype = "Layers"  # Default layered atmosphere
    smooth = False  # Default no smoothing
    igamma = -1  # Index of gamma if included
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

    for j in range(0, len(lines)):

        if str(lines[j]) == "Basic\n":
            state = 0
            b1 = i
        elif str(lines[j].split()[0]) == "Gases":
            state = 1
            g1 = i
            if len(lines[j].split()) > 1:
                gases.append(lines[j].split()[1])
            else:
                gases.append(
                    "h2he"
                )  # Filler gas is h2+he, may get more reliable results from h2 only
        elif str(lines[j].split()[0]) == "Atm":
            state = 2
            a1 = i
            atmtype = lines[j].split()[1]
            if len(lines[j].split()) > 2:
                if str(lines[j].split()[2]) == "Verbatim":
                    verbatim = True
        elif str(lines[j].split()[0]) == "Haze" or str(lines[j].split()[0]) == "Clouds":
            state = 3
            c1 = i
            cloudmod = int(lines[j].split()[1])
            if len(lines[j].split()) >= 3:
                hazetype = 0
                hazestr = str(lines[j].split()[2])
                if hazestr in hazelist:
                    hazetype = hazelist.index(hazestr)
        elif str(lines[j]) == "End\n":
            state = 4
            e1 = i

        elif len(lines[j].split()) < 6:
            print("Error: missing parameter values.")
            sys.exit()

        else:
            line = lines[j].split()
            pnames.append(line[0])
            if state == 0:
                basic.append(line[0])
                bnum = bnum + 1
            if state == 1:
                gases.append(line[0])
                gnum = gnum + 1
            if state == 2:
                atm.append(line[0])
                anum = anum + 1
            if state == 3:
                clouds.append(line[0])
                cnum = cnum + 1
            if state == 4:
                end.append(line[0])
                enum = enum + 1
            plparams[i] = (float)(line[1])
            guess[i] = (float)(line[1])
            mu[i] = (float)(line[2])
            sigma[i] = (float)(line[3])
            bounds[i, 0] = (float)(line[4])
            bounds[i, 1] = (float)(line[5])
            if len(line) >= 9:
                prior_functions[i] = line[8]

            if sigma[i] > 0.0:
                pvars.append(plparams[i])
                nvars.append(i)

            # if guess[i]==0: guess[i] = 0.1*bounds[i,1]        # Prevents errors from occuring where parameters are zero.

            # Convert pressure units from bars to cgs
            if pnames[i] == "Cloud_Base" or pnames[i] == "P_cl":
                plparams[i] = plparams[i] + 6.0
                guess[i] = guess[i] + 6.0
                mu[i] = mu[i] + 6.0
                bounds[i, 0] = bounds[i, 0] + 6.0
                bounds[i, 1] = bounds[i, 1] + 6.0
                if bounds[i, 0] < minP:
                    bounds[i, 0] = minP
                if bounds[i, 1] > maxP:
                    bounds[i, 1] = maxP

            if line[0] == "gamma":
                smooth = True
                igamma = j
            # ada: We want to impose a normal prior on log g,
            # while keeping uniform priors on everything else.
            elif line[0] == "Log(g)":
                ilogg = i
            if len(line) > 6 and line[6] == "Ensemble":
                ensparams.append(i)
            i = i + 1

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
        opacdir=opacdir,
    )
