import numpy as np
from numpy.typing import ArrayLike
from os import PathLike
import pickle
from scipy.interpolate import PchipInterpolator as monotonic_interpolation
from scipy.ndimage import gaussian_filter1d as gaussian_smoothing
from typing import Callable, Union


Jupiter_radius_in_Earth_radii = 11.2

parameter_names = [
    'Rad',
    'Log(g)',
    'h2o',
    'co',
    'co2',
    'ch4',
    'Lupu_alk',
    'h2s',
    'nh3',
    'T_m4',
    'T_m3',
    'T_m2',
    'T_m1',
    'T_0',
    'T_0p5',
    'T_1',
    'T_1p5',
    'T_2',
    'T_2p5',
    'deltaL',
    'scaleG395_ch1',
    'scaleG395_ch2',
    'logf', # end of free parameters
    'Mass',
    'C/O',
    '[Fe/H]',
    'Teff'
]

basic_index = slice(0, 2)
gas_index = slice(2, 9)
TP_index = slice(9, 19)
calibration_index = slice(19, 23)

# log(P/bar)
TP_resolution = dict(number_of_pressures=71, minimum_log_pressure=-4, maximum_log_pressure=2.5)


def process_dynesty_samples(
        dynesty_results,
        derived_parameters,
        parameter_names,
        importance_percentile=0.95
        ):
    log_importance_weights = dynesty_results.logwt
    importance_weights = np.exp(log_importance_weights - np.max(log_importance_weights))
    cumsum_logwt = np.cumsum(importance_weights)
    important_iterations = cumsum_logwt/cumsum_logwt[-1] >= (1-importance_percentile)

    important_samples = dynesty_results.samples[important_iterations]
    important_derived = derived_parameters[important_iterations]
    parameters = np.concatenate((important_samples, important_derived), axis=-1)

    log_likelihoods = dynesty_results.logl[important_iterations]

    sample_weights = importance_weights[important_iterations]

    radius_index = parameter_names.index("Rad")
    parameters[radius_index] /= Jupiter_radius_in_Earth_radii 

    return dict(parameters=parameters, log_likelihoods=log_likelihoods, weights=sample_weights)


def piette(T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2, T_2p5,
           number_of_pressures, minimum_log_pressure, maximum_log_pressure):
    logP_nodes = np.array([-4, -3, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5])
    T_nodes = np.array([T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2, T_2p5])

    interpolated_function = monotonic_interpolation(logP_nodes, T_nodes)
    input_Ps = np.linspace(minimum_log_pressure, maximum_log_pressure, num=number_of_pressures)

    output_Ts = gaussian_smoothing(interpolated_function(input_Ps), sigma=0.3)

    return output_Ts


def calculate_CtoO_and_metallicity(gases, gas_logabundances):
    carbon = 0.
    oxygen = 0.
    metals = 0.
    ccompounds = ['ch4','co','co2','hcn']
    cmult = [1.,1.,1.,1.]
    ocompounds = ['h2o','co','co2','tio','vo']
    omult = [1.,1.,2.,1.,1.]
    zcompounds = ['h2o','ch4','co','co2','nh3','h2s','Burrows_alk','Lupu_alk','crh','feh','tio','vo','hcn','n2','ph3']
    zmult = [16.,12.,28.,44.,14., 32.,24.,24.,52.,56., 64.,67.,26.,28.,31.]

    for i in np.arange(0,len(ccompounds)):
        if ccompounds[i] in gases:
            j = gases.index(ccompounds[i])
            carbon = carbon + cmult[i]*(10**gas_logabundances[j]) # -1 because of hydrogen
    for i in np.arange(0,len(ocompounds)):
        if ocompounds[i] in gases:
            j = gases.index(ocompounds[i])
            oxygen = oxygen + omult[i]*(10**gas_logabundances[j])
    for i in np.arange(0,len(zcompounds)):
        if zcompounds[i] in gases:
            j = gases.index(zcompounds[i])
            metals = metals + zmult[i]*(10**gas_logabundances[j])

    ctoo = carbon/oxygen
    fetoh = np.log10(metals/0.0196)

    return np.array([ctoo, fetoh])


def extract_dynesty_samples(
        sampler_results_filepath: Union[PathLike, str],
        derived_parameters_filepath: Union[PathLike, str],
        TP_function: Callable[list[float], ArrayLike] = piette,
        TP_resolution_kwargs: dict[str, Union[int, float]] = TP_resolution
        ) -> dict[str, ArrayLike]:
    with open(sampler_results_filepath, "rb") as results_file:
        dynesty_results = pickle.load(results_file)

    with open(derived_parameters_filepath, "rb") as derived_file:
        derived_parameters = np.asarray(pickle.load(derived_file))

    parameter_dict = process_dynesty_samples(dynesty_results, derived_parameters, parameter_names)

    # max_likelihood_index = np.argmax(parameter_dict["log_likelihoods"])
    TP_parameters = parameter_dict["parameters"][:, TP_index].T
    TP_profiles = TP_function(*TP_parameters, **TP_resolution_kwargs)

    return parameter_dict | dict(TP=TP_profiles)