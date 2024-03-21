from __future__ import annotations

from astropy.constants import R_earth, R_jup
from astropy.units import cm
from src.ApolloFunctions import GetScaOpac
from dataclasses import dataclass, field, InitVar
from functools import partial
import numpy as np
from scipy.interpolate import interp1d
from src.wrapPlanet import PyPlanet
from typing import Callable, NamedTuple

from user.TP_models import TP_models


# REarth_in_cm = R_earth.to(cm).value
REarth_in_cm = 6.371e8

# RJup_in_REarth = (R_jup/R_earth).decompose().value
RJup_in_REarth = 11.2


class CPlanet_constructor(NamedTuple):
    observation_mode_index:       int
    cloud_model_index:            int
    hazetype_index:               int
    number_of_streams:            int
    TP_model_switch:              int
    model_wavelengths:            list[float]
    Teff_calculation_wavelengths: list[float]
    gas_species_indices:          list[int]
    gas_opacity_directory:        str
    opacity_catalog_name:         str
    Teff_opacity_catalog_name:    str


class MakeModel_constructor(NamedTuple):
    model_wavelengths:             list[float]
    gas_species:                   list[str]
    minimum_model_pressure:        float
    maximum_model_pressure:        float
    model_pressures:               list[float]
    number_of_atmospheric_layers:  int
    stellar_effective_temperature: float
    stellar_radius:                float
    semimajor_axis:                float
    distance_to_system:            float
    gas_start_index:               int
    gas_end_index:                 int
    TP_model_type:                 str
    TP_model_function:             Callable[[...], list[float]]
    is_TP_profile_gray:            bool
    isothermal_temperature:        float
    is_TP_profile_verbatim:        bool
    TP_start_index:                int
    TP_end_index:                  int
    number_of_TP_arguments:        int
    number_of_params1_arguments:   int
    basic_parameter_names:         list[str]
    basic_start_index:             int
    cloud_model_index:             int
    hazetype_index:                str
    cloud_parameter_names:         list[str]
    cloud_start_index:             int
    calibration_parameter_names:   list[str]
    calibration_start_index:       int
    APOLLO_use_mode:               str


#emission_flux:                               list[float]
#radius_parameter:                            float
#radius_case:                                 str
#wavelength_offset:                           float
#model_parameters:                            list[float]
#planet_parameters:                           list[float]
#free_variable_indices:                       list[int]
#basic_parameter_names:                       list[str]
#basic_start_index:                           int
#end_parameter_names:                         list[str]
#end_start_index:                             int


class ObserveModel_constructor(NamedTuple):
    data_wavelengths:                            list[float]
    model_wavelengths:                           list[float]
    model_indices_by_opacity_bins:               list[float]
    distance_to_system:                          float
    bin_low_wavelength_indices_by_opacity_bins:  list[float]
    bin_high_wavelength_indices_by_opacity_bins: list[float]
    bin_low_wavelength_offset_by_index:          list[float]
    bin_high_wavelength_offset_by_index:         list[float]
    binned_data_flux:                            list[float]
    data_downbinning_factor:                     float
    data_convolution_factor:                     float
    data_band_bounding_indices:                  list[list[int]]
    observation_mode_index:                      int
    normalize_spectrum:                          bool
    use_default_radius:                          bool
    data_errors_binned_to_model:                 list[float]


@dataclass
class Planet:
    """Class for wrapping the Planet class as defined in C++."""
    name: str
    cclass_constructor: InitVar[CPlanet_constructor]
    cclass: PyPlanet = field(init=False)

    def __post_init__(self, cclass_constructor: CPlanet_constructor):
        self.cclass = PyPlanet()
        cclass_args = [
            list(cclass_constructor[:5]),
            *list(cclass_constructor[5:8]),
            *[arg.encode('utf-8') for arg in cclass_constructor[8:]]
            ]
        self.cclass.MakePlanet(*cclass_args)

    def MakeModel(self, model_constructor: MakeModel_constructor):
        def ModelFunction(
                model_parameters:              list[float],
                ### GetModel_kwargs below
                model_wavelengths:             list[float],
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
                is_TP_profile_verbatim:        bool,
                TP_start_index:                int,
                TP_end_index:                  int,
                number_of_TP_arguments:        int,
                number_of_params1_arguments:   int,
                basic_parameter_names:         list[str],
                basic_start_index:             int,
                cloud_model_index:             int,
                hazetype_index:                str,
                cloud_parameter_names:         list[str],
                cloud_start_index:             int,
                calibration_parameter_names:   list[str],
                calibration_start_index:       int,
                APOLLO_use_mode:               str
                ):
            ##########################################
            # Start by relabeling all the arguments. #
            ##########################################
            x = model_parameters
            modwave = model_wavelengths
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
            cloudmod = cloud_model_index
            hazetype = hazetype_index
            clouds = cloud_parameter_names
            c1 = cloud_start_index
            end = calibration_parameter_names
            e1 = calibration_start_index
            task = APOLLO_use_mode

            ##########################################
            make_gaussian_profile = lambda x: x - 0.5*(np.arange(npress:=np.shape(x)[-1])-npress/2)**2/((npress)**2)
            
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
                # abund = make_gaussian_profile(np.asarray([abund]*number_of_atmospheric_layers).T)
                abund[0] = 1. - asum
                mmw,rxsec = GetScaOpac(gases, abund[1:])
            print(f"{abund=}")
            self.cross_sections = rxsec

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

                self.cclass.set_Params(params1,abund,rxsec,tplong)
                if cloud_fraction == 1:
                    specflux = self.cclass.get_Spectrum()
                else:
                    specflux = cloud_fraction*np.asarray(self.cclass.get_Spectrum()) + (1-cloud_fraction)*np.asarray(self.cclass.get_ClearSpectrum())

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
                    self.cclass.set_Params(params1,abund,rxsec,tplong)
                    if cloud_fraction == 1:
                        specflux = self.cclass.get_Spectrum()
                    else:
                        specflux = cloud_fraction*np.asarray(self.cclass.get_Spectrum()) + (1-cloud_fraction)*np.asarray(self.cclass.get_ClearSpectrum())
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
                    self.cclass.set_Params(params1,abund,rxsec,tplong)
                    if cloud_fraction == 1:
                        specflux = self.cclass.get_Spectrum()
                    else:
                        specflux = cloud_fraction*np.asarray(self.cclass.get_Spectrum()) + (1-cloud_fraction)*np.asarray(self.cclass.get_ClearSpectrum())
                if atmtype == 'Parametric':
                    # Compute spectrum
                    self.cclass.set_Params(params1,abund,rxsec,tpprof)
                    if cloud_fraction == 1:
                        specflux = self.cclass.get_Spectrum()
                    else:
                        specflux = cloud_fraction*np.asarray(self.cclass.get_Spectrum()) + (1-cloud_fraction)*np.asarray(self.cclass.get_ClearSpectrum())

            if atmtype == "Piette":
                monotonic = (np.diff(x[a1:a2]) >= 0)
                if not np.all(monotonic):
                    print("Profile isn't purely monotonic. The array: {}".format(x[a1:a2]))
                    print("Truth array: {}".format(monotonic))

            teff = self.cclass.get_Teff()
            if task!='Ensemble':
                print('M/Mj: ',mass)
                print('C/O: ',ctoo)
                print('[Fe/H]: ',fetoh)
                print('Teff: ',teff)

            if 'scaleJ' in end:
                pos = end.index('scaleJ')
                scaleJ = x[e1+pos]
            else:
                scaleJ = 1.0
            if 'scaleH' in end:
                pos = end.index('scaleH')
                scaleH = x[e1+pos]
            else:
                scaleH = 1.0
            if 'scaleK' in end:
                pos = end.index('scaleK')
                scaleK = x[e1+pos]
            else:
                scaleK = 1.0
            if 'scaleG395' in end:
                pos = end.index('scaleG395')
                scaleG395 = x[e1+pos]
            else:
                scaleG395 = 1.0
            if 'scaleG395_ch1' in end:
                pos = end.index('scaleG395_ch1')
                scaleG395_ch1 = x[e1+pos]
            else:
                scaleG395_ch1 = 1.0
            if 'scaleG395_ch2' in end:
                pos = end.index('scaleG395_ch2')
                scaleG395_ch2 = x[e1+pos]
            else:
                scaleG395_ch2 = 1.0

            specflux = np.asarray(specflux)
            # ada: scale band data
            J_boundaries = [1.10, 1.36]
            H_boundaries = [1.44, 1.82]
            K_boundaries = [1.94, 2.46]
            G395_boundaries = [2.8, 5.3]
            G395_ch1_boundaries = [2.8, 4.05]
            G395_ch2_boundaries = [4.15, 5.3]

            wavelengths = modwave
            specflux = np.where(np.logical_and(J_boundaries[0]<=wavelengths, wavelengths<=J_boundaries[1]), specflux*scaleJ, specflux)
            specflux = np.where(np.logical_and(H_boundaries[0]<=wavelengths, wavelengths<=H_boundaries[1]), specflux*scaleH, specflux)
            specflux = np.where(np.logical_and(K_boundaries[0]<=wavelengths, wavelengths<=K_boundaries[1]), specflux*scaleK, specflux)
            specflux = np.where(np.logical_and(G395_boundaries[0]<=wavelengths, wavelengths<=G395_boundaries[1]), specflux*scaleG395, specflux)
            specflux = np.where(np.logical_and(G395_ch1_boundaries[0]<=wavelengths, wavelengths<=G395_ch1_boundaries[1]), specflux*scaleG395_ch1, specflux)
            specflux = np.where(np.logical_and(G395_ch2_boundaries[0]<=wavelengths, wavelengths<=G395_ch2_boundaries[1]), specflux*scaleG395_ch2, specflux)
            # Plot the results of get_Spectrum() directly for testing purposes.
            '''
            figtest = plt.figure()
            ax = figtest.add_subplot()
            ax.plot(modwave,specflux)
            plt.savefig('plots/test.png')
            sys.exit()
            '''

            return specflux, [mass, ctoo, fetoh, teff]

        return partial(ModelFunction, **model_constructor._asdict())
