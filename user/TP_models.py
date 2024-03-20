import numpy as np
from user.defaults import minP, maxP, minT, maxT, tgray, vres
from scipy.interpolate import PchipInterpolator as monotonic_interpolation
from scipy.ndimage import gaussian_filter1d as gaussian_smoothing

from user.model_module import ModelModule

from user.priors import priors


###############################################################################
######################### Define model functions here. ########################
###############################################################################


########## ISOTHERMAL Profile ##########
# A gray (isothermal) function.
def gray(T=tgray, num_layers_final=vres, P_min=minP, P_max=maxP):
    profile = np.ones(vres) * T

    return profile


########## VERBATIM Profile ##########
# Placeholder for a direct layer-by-layer profile. Will place the input
# temperatures evenly spaced in log-pressure between min and max pressure.
def verbatim(*input_Ts, num_layers_final=vres, P_min=minP, P_max=maxP):
    profile = input_Ts

    return profile


########## SIMPLE INTERPOLATION Profile ##########
# Straightforward layer-by-layer interpolation. Uses linear interpolation
# unless enough interpolation points warrant a cubic.
def interpolation(*input_Ts, num_layers_final=vres, P_min=minP, P_max=maxP):
    num_layers_initial = len(input_Ts)
    input_Ps = np.linspace(P_max, P_min, num=num_layers_initial)
    output_Ps = np.linspace(P_max, P_min, num=num_layers_final)

    interpolation_level = "cubic" if num_layers_final > 4 else "linear"

    output_Ts = (np.interp1d(input_Ps, input_Ts, kind=interpolation_level))(output_Ps)
    output_Ts = np.where(output_Ts < minT, minT, output_Ts)
    profile = np.where(output_Ts > maxT, maxT, output_Ts)

    return profile


########## POWER LAW+LINEAR Profile ##########
# A parametrization I tried to use for L dwarf retrievals. Loosely based on
# the Madhusudhan & Seager (2009) prescription with a free power law exponent
# for the shallower of 2 portions, connected via smooth 1st derivative to a
# deeper linear component. I now prefer to use something more flexible, like
# Piette's.
def power_linear(
    exponent, logPiso, T0, T_mid, T_max, num_layers_final=vres, P_min=minP, P_max=maxP
):
    logP0 = P_min
    T_frac = exponent * (T_mid - T0) / (T_max - T_mid)
    logPmid = (logP0 + T_frac * logPiso) / (1 + T_frac)
    alpha2 = (logPiso - logPmid) / (T_max - T_mid)
    alpha1 = exponent * (T_mid - T0) ** (1 - 1 / exponent) * alpha2

    input_Ps = np.linspace(P_min, P_max, num=num_layers_final)

    output_Ts = np.ones_like(input_Ps) * T0
    output_Ts = np.where(
        (input_Ps > logP0) & (input_Ps < logPmid),
        T0 + ((input_Ps - logP0) / alpha1) ** exponent,
        output_Ts,
    )
    output_Ts = np.where(
        (input_Ps > logPmid) & (input_Ps < logPiso),
        T_mid + ((input_Ps - logPmid) / alpha2),
        output_Ts,
    )
    profile = np.where(input_Ps >= logPiso, T_max, output_Ts)

    return profile


########## ANJALI PIETTE et. al. Profile ##########
# A specific interpolation with its own smoothing. Should be flexible for
# most retrieval cases without thermal inversions. See for reference
# Piette, Anjali A. A., and Nikku Madhusudhan. “Considerations for Atmospheric
# Retrieval of High-Precision Brown Dwarf Spectra” 19, no. July (July 29, 2020):
# 1–19. http://arxiv.org/abs/2007.15004.
def piette(
    T_m4,
    T_m3,
    T_m2,
    T_m1,
    T_0,
    T_0p5,
    T_1,
    T_1p5,
    T_2,
    T_2p5,
    num_layers_final=vres,
    P_min=minP,
    P_max=maxP,
):
    logP_nodes = np.array([-4, -3, -2, -1, 0, 0.5, 1, 1.5, 2])
    T_nodes = np.array([T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2])

    interpolated_function = monotonic_interpolation(logP_nodes, T_nodes)
    input_Ps = np.linspace(P_min, P_max, num=num_layers_final)

    profile = gaussian_smoothing(interpolated_function(input_Ps), sigma=0.3)

    return profile


###############################################################################


###############################################################################
########### Define any special prior functions here, if applicable. ###########
###############################################################################

# Given a list of prior functions, one per parameter of your profile function,
# are there dependencies between the priors? This function should set the
# parameters of the prior distributions given the specific profile's needs.


def sample_dependent_prior(
    distribution_functions,
    distribution_parameter_sets,
    input_values,
    sampler_type,
    logP_index,
    bound_bounds,
):
    # If there are user-specified priors in individual temperature nodes that aren't
    # just the min/max possible temperatures, check to see if they are within the bounds
    # of the dependent priors. If not, discard them.
    set_bound = lambda initial, minimum, maximum: np.min(
        [np.max([initial, minimum, minT]), maximum, maxT]
    )

    # The lower and upper bounds are expected to be the final 2 parameters.
    lower_bounds, upper_bounds = np.asarray(distribution_parameter_sets)[:, -2:].T
    lower_bound = set_bound(lower_bounds[logP_index], *bound_bounds)
    upper_bound = set_bound(upper_bounds[logP_index], *bound_bounds)
    distribution_parameter_set = np.r_[
        np.asarray(distribution_parameter_sets)[logP_index, :-2],
        lower_bound,
        upper_bound,
    ]

    prior = priors[distribution_functions[logP_index]]
    prior.generate_prior(distribution_parameter_set)
    T = prior.generate_sample_from_input(input_values[logP_index], sampler_type)

    return prior, T


########## ANJALI PIETTE et. al. Prior ##########
# For isolated brown dwarfs or distant planets/BDs where thermal inversions
# are not expected. I assume you will use uniform priors for each of the
# temperature nodes but in principle any bounded prior should work.
#################################################
def piette_monotonic_prior(
    distribution_functions, distribution_parameter_sets, input_values, sampler_type
):
    # Assume 10 temperature nodes: log(P/bar) = [-4, -3, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5]

    # logP=0.5 starts with the priors input by the user.
    logP_0p5 = 5
    prior_logP_0p5, T_logP_0p5 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_0p5,
        [minT, maxT],
    )

    # logP=-4 has an upper prior bound set by logP=0.
    logP_m4 = 0
    prior_logP_m4, T_logP_m4 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_m4,
        [minT, T_logP_0p5],
    )

    # logP=-2's bounds are themselves bounded by logP=-4 and logP=0.
    logP_m2 = 2
    prior_logP_m2, T_logP_m2 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_m2,
        [T_logP_m4, T_logP_0p5],
    )

    # logP=-3's bounds are themselves bounded by logP=-4 and logP=-2.
    logP_m3 = 1
    prior_logP_m3, T_logP_m3 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_m3,
        [T_logP_m4, T_logP_m2],
    )

    # logP=-1's bounds are themselves bounded by logP=-2 and logP=0.5.
    logP_m1 = 3
    prior_logP_m1, T_logP_m1 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_m1,
        [T_logP_m2, T_logP_0p5],
    )

    # logP=0's bounds are themselves bounded by logP=-1 and logP=0.5.
    logP_0 = 4
    prior_logP_0, T_logP_0 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_0,
        [T_logP_m1, T_logP_0p5],
    )

    # logP=2 has a lower prior bound set by logP=0.5.
    logP_2p5 = 9
    prior_logP_2p5, T_logP_2p5 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_2p5,
        [T_logP_0p5, maxT],
    )

    # logP=1.5's bounds are themselves bounded by logP=0.5 and logP=2.5.
    logP_1p5 = 7
    prior_logP_1p5, T_logP_1p5 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_1p5,
        [T_logP_0p5, T_logP_2p5],
    )

    # logP=2's bounds are themselves bounded by logP=1 and logP=2.
    logP_2 = 8
    prior_logP_2, T_logP_2 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_2,
        [T_logP_1p5, T_logP_2p5],
    )

    # logP=1's bounds are themselves bounded by logP=0.5 and logP=1.5.
    logP_1 = 6
    prior_logP_1, T_logP_1 = sample_dependent_prior(
        distribution_functions,
        distribution_parameter_sets,
        input_values,
        sampler_type,
        logP_1,
        [T_logP_0p5, T_logP_1p5],
    )

    proposed_profile = [
        T_logP_m4,
        T_logP_m3,
        T_logP_m2,
        T_logP_m1,
        T_logP_0,
        T_logP_0p5,
        T_logP_1,
        T_logP_1p5,
        T_logP_2,
        T_logP_2p5,
    ]

    priors = [
        prior_logP_m4,
        prior_logP_m3,
        prior_logP_m2,
        prior_logP_m1,
        prior_logP_0,
        prior_logP_0p5,
        prior_logP_1,
        prior_logP_1p5,
        prior_logP_2,
        prior_logP_2p5,
    ]

    return np.array(
        [
            prior.evaluate_prior_from_sample(proposed_node, sampler_type)
            for proposed_node, prior in zip(proposed_profile, priors)
        ]
    )


###############################################################################


###############################################################################
##################### The dictionary to index the models. #####################
###############################################################################

TP_models = {
    "gray": ModelModule(gray),
    "verbatim": ModelModule(verbatim),
    "interpolation": ModelModule(interpolation),
    "power_linear": ModelModule(power_linear),
    "piette": ModelModule(piette, piette_monotonic_prior),
}

###############################################################################
