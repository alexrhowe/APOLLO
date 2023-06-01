from ApolloFunctions import GetScaOpac
from collections.abc import Callable
import numpy as np
from scipy.interpolate import interp1d
from src.wrapPlanet import PyPlanet as Planet
from user.TP_models import TP_models

REarth_in_cm = 6.371e8
parsec_in_cm = 3.086e18
RJup_in_REarth = 11.2


def MakePlanet(
        observation_type:             str,
        TP_model_type:                str,
        cloud_model_type:             int,
        haze_type:                    str,
        gas_species:          list[str],
        number_of_streams:            int,
        gas_opacity_directory:        str,
        opacity_catalog_name:         str,
        model_wavelengths:            list[list[float]],
        Teff_opacity_catalog_name:    str,
        Teff_calculation_wavelengths: list[list[float]],
        ) -> Planet:

    switches = [
        observation_type,
        cloud_model_type,
        haze_type,
        number_of_streams
        ]

    #opacities = get_opacities_from_catalog(
    #    gas_opacity_directory,
    #    opacity_catalog_name,
    #    list_of_gas_species
    #    )

    #list_of_molecular_species = list(opacities.index)
    #default_opacity_entry = opacities.loc[opacities.index[0]]

    planet = Planet()
    planet.MakePlanet(
        switches,
        model_wavelengths,
        Teff_calculation_wavelengths,
        gas_species,
        gas_opacity_directory,
        opacity_catalog_name.encode('utf-8'),
        Teff_opacity_catalog_name.encode('utf-8')
        )

    return planet


def GetModel(
        model_parameters:              list[float],
        planet:                        Planet,
        planet_parameters:             list[...],
        gas_species:                   list[str],
        minimum_model_pressure:        float,
        maximum_model_pressure:        float,
        model_pressures:               list[float],
        number_of_atmospheric_layers:  int,
        stellar_effective_temperature: float,
        stellar_radius:                float,
        semimajor_axis:                float,
        distance_to_system:            float,
        gas_start_index:               int,
        gas_end_index:                 int,
        TP_model_type:                 str,
        TP_model_function:             Callable[[...], list[float]],
        is_TP_profile_gray:            bool,
        isothermal_temperature:        float,
        is_TP_profile_verbatim:        bool
        TP_start_index:                int,
        TP_end_index:                  int,
        number_of_TP_arguments:        int,
        number_of_params1_arguments:   int,
        basic_parameter_names:         list[str],
        basic_start_index:             int,
        cloud_model_type:              int,
        haze_type:                     str,
        cloud_parameter_names:         list[str],
        cloud_start_index:             int,
        APOLLO_use_mode:               str,
        ):
    ##########################################
    # Start by relabeling all the arguments. #
    ##########################################
    x = model_parameters
    plparams = planet_parameters
    gases = gas_species
    minP = minimum_model_pressure
    maxP = maximum_model_pressure
    profin = model_pressures
    vres = number_of_atmospheric_layers
    tstar = stellar_effective_temperature
    rstar = stellar_radius
    sma = semimajor_axis
    dist = distance_to_system
    g1 = gas_start_index
    g2 = gas_end_index
    atmtype = TP_model_type
    TP_model = TP_model_function
    gray = is_TP_profile_gray
    tgray = isothermal_temperature
    verbatim = is_TP_profile_verbatim
    a1 = TP_start_index
    a2 = TP_end_index
    natm = number_of_TP_arguments
    ilen = number_of_params1_arguments
    basic = basic_parameter_names
    b1 = basic_start_index
    cloudmod = cloud_model_type
    hazetype = haze_type
    clouds = cloud_parameter_names
    c1 = cloud_start_index
    task = APOLLO_use_mode

    ##########################################

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
        mmw,rxsec = GetScaOpac(gases, abund[1:])

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
        radius = REarth_in_cm*x[b1+pos]
    elif 'RtoD' in basic:
        pos = basic.index('RtoD')
        params1[0] = 10**x[b1+pos]*dist*4.838e9 # convert R/D to Earth radii
        radius = 10**x[b1+pos]*dist*4.838e9*REarth_in_cm
    elif 'RtoD2U' in basic:
        pos = basic.index('RtoD2U')
        params1[0] = np.sqrt(x[b1+pos])
        radius = np.sqrt(x[b1+pos])*REarth_in_cm
    else:
        global norad
        norad = True
        params1[0] = RJup_in_REarth
        radius = RJup_in_REarth*REarth_in_cm
        # Default radius = Jupiter

    # Gravity handling
    if 'Log(g)' in basic:
        pos = basic.index('Log(g)')
        params1[1] = x[b1+pos]
        grav = 10**x[b1+pos]
    else:
        params1[1] = 4.1
        grav = 4.1

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

    # ada: For fractional cloud coverage, we generate the cloud-free
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
        
    tpprof = np.zeros(natm)
    # Gray atmosphere approximation
    if gray:
        tplong = TP_model.evaluate_model(tgray,
                                         num_layers_final=vres,
                                         P_min=minP,
                                         P_max=maxP)

        planet.set_Params(params1,abund,rxsec,tplong)
        if cloud_fraction == 1:
            specflux = planet.get_Spectrum()
        else:
            specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())

    # Build atmosphere temperature profile and compute spectrum
    else:
        for i in range(0,len(tpprof)): tpprof[i] = x[i+a1]
        if atmtype=='Parametric': tpprof[1] = 10**tpprof[1]
        if atmtype in TP_models:
            tplong = TP_model.evaluate_model(*tpprof,
                                             num_layers_final=vres,
                                             P_min=minP-6,
                                             P_max=maxP-6)
            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tplong)
            if cloud_fraction == 1:
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
                # ada: having a cubic interpolation with too many points can
                # introduce potentially unwarranted wiggles in the profile,
                # which can greatly affect the resulting spectrum.
                # if natm<=4: f = interp1d(profin,tpprof,kind='linear')
                if natm<=4 or natm>=13: f = interp1d(profin,tpprof,kind='linear')
                else: f = interp1d(profin,tpprof,kind='cubic')
                for i in range(0,vres):
                    tplong[i] = f(maxP + (minP-maxP)*i/(float)(vres-1))
                    if(tplong[i]<75.): tplong[i]=75.
                    if(tplong[i]>4000.): tplong[i]=4000.
            # The layer-by-layer profile is input from largest to smallest
            # pressures, so we need to reverse the profile before passing
            # it to the C++ side.
            tplong = tplong[::-1]
            # np.save(outdir+outfile+"T-P_array_linear", tplong) # vestigial

            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tplong)
            if cloud_fraction == 1:
                specflux = planet.get_Spectrum()
            else:
                specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())
        if atmtype == 'Parametric':
            # Compute spectrum
            planet.set_Params(params1,abund,rxsec,tpprof)
            if cloud_fraction == 1:
                specflux = planet.get_Spectrum()
            else:
                specflux = cloud_fraction*np.asarray(planet.get_Spectrum()) + (1-cloud_fraction)*np.asarray(planet.get_ClearSpectrum())

    if atmtype == "Piette":
        monotonic = (np.diff(plparams[a1:a2]) >= 0)
        if not np.all(monotonic):
            print("Profile isn't purely monotonic. The array: {}".format(plparams[a1:a2]))
            print("Truth array: {}".format(monotonic))

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

    return specflux, [mass, ctoo, fetoh, teff]