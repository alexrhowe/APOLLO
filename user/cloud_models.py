import numpy as np
from user.defaults import minP, maxP, minT, maxT, tgray, vres

from user.model_module import ModelModule

from user.priors import priors, evaluate_default_priors


###############################################################################
######################### Define model functions here. ########################
###############################################################################


########## VERBATIM Profile ##########
# Placeholder for a direct layer-by-layer profile. Will place the input
# temperatures evenly spaced in log-pressure between min and max pressure.
def verbatim(parameters):

    return parameters


###############################################################################
########### Define any special prior functions here, if applicable. ###########
###############################################################################

# Given a list of prior functions, one per parameter of your profile function,
# are there dependencies between the priors? This function should set the
# parameters of the prior distributions given the specific profile's needs.


def sample_dependent_prior(distribution_functions, distribution_parameter_sets, input_values, sampler_type,
                           logP_index, bound_bounds):
    # If there are user-specified priors in individual temperature nodes that aren't
    # just the min/max possible temperatures, check to see if they are within the bounds
    # of the dependent priors. If not, discard them.
    minP_bar = -4
    maxP_bar = 2.5
    set_bound = lambda initial, minimum, maximum: np.min([np.max([initial, minimum, minP_bar]), maximum, maxP_bar])

    # The lower and upper bounds are expected to be the final 2 parameters.
    lower_bounds, upper_bounds = np.asarray(distribution_parameter_sets)[:, -2:].T
    lower_bound = set_bound(lower_bounds[logP_index], *bound_bounds)
    upper_bound = set_bound(upper_bounds[logP_index], *bound_bounds)
    distribution_parameter_set = np.r_[np.asarray(distribution_parameter_sets)[logP_index, :-2],
                                       lower_bound, upper_bound]

    prior = priors[distribution_functions[logP_index]]
    prior.generate_prior(distribution_parameter_set)
    T = prior.generate_sample_from_input(input_values[logP_index], sampler_type)

    return prior, T


########## ANJALI PIETTE et. al. Prior ##########
# For isolated brown dwarfs or distant planets/BDs where thermal inversions
# are not expected. I assume you will use uniform priors for each of the
# temperature nodes but in principle any bounded prior should work.
#################################################
def cloud_base_in_model(distribution_functions, distribution_parameter_sets, input_values, sampler_type):
    # Set the maximum thickness of the cloud layer based on the base
    # (deepest pressure) of the model.

    cloud_default_priors = evaluate_default_priors(distribution_functions, distribution_parameter_sets, input_values, sampler_type)

    # Hard coding the pressure limits in bars for now.
    minP_bar = -4
    maxP_bar = 2.5
    cloud_top_index = 1
    cloud_thickness_index = 2
    maximum_thickness = maxP_bar - cloud_default_priors[cloud_top_index]
    prior, logP_thickness = sample_dependent_prior(distribution_functions, distribution_parameter_sets, input_values, sampler_type, cloud_thickness_index, [0, maximum_thickness])

    cloud_default_priors[cloud_thickness_index] = logP_thickness

    return cloud_default_priors

###############################################################################


###############################################################################
##################### The dictionary to index the models. #####################
###############################################################################

cloud_models = {

    "verbatim": ModelModule(verbatim, cloud_base_in_model)

}

###############################################################################
