import numpy as np
from scipy.stats import norm, truncnorm, uniform


###############################################################################
############### Define distribution functions for priors here. ################
###############################################################################

# APOLLO currently has 4 slots for prior parameters, which are nominally mean,
# standard deviation, and lower and upper bounds. But you can use whatever
# you need in these 4 slots.

def uniform_distribution(distribution_parameters):
    _, _, uniform_min, uniform_max = distribution_parameters
    uniform_width = uniform_max - uniform_min

    return uniform(uniform_min, uniform_width)


def normal_distribution(distribution_parameters):
    mean, sigma, _, _ = distribution_parameters

    return norm(mean, sigma)


def truncated_normal_distribution(distribution_parameters):
    mean, sigma, lower_bound, upper_bound = distribution_parameters
    a, b = (lower_bound - mean) / sigma, (upper_bound - mean) / sigma

    return truncnorm(a, b, mean, sigma)

###############################################################################


###############################################################################
############ A class to load priors given a distribution function. ############
###############################################################################

class Prior:
    def __init__(self, prior_distribution):
        self.prior_distribution = prior_distribution


    def generate_prior(self, distribution_parameters):
        self.prior = self.prior_distribution(distribution_parameters)

        
    def generate_sample_from_input(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee sends a list of proposed parameters, so no need to do
            # anything to generate a proposal.
            return input_value

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior given a
            # uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")


    def evaluate_prior_from_sample(self, input_sample, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the log of the PDF for its prior calculation.
            return self.prior.logpdf(input_sample)

        elif sampler_type == "dynesty":

            # dynesty uses the proposed parameter values, so no need to do
            # anything if the proposed sample is already generated.
            return input_sample

        else:
            raise ValueError("Sampler type not valid.")

    def evaluate_prior_from_input(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the log of the PDF for its prior calculation.
            return self.prior.logpdf(input_value)

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior given a
            # uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")

###############################################################################


###############################################################################
##################### The dictionary to index the priors. #####################
###############################################################################

priors = {

    "Uniform":          Prior(uniform_distribution),

    "Normal":           Prior(normal_distribution),

    "Truncated_Normal": Prior(truncated_normal_distribution)

}

###############################################################################


###############################################################################
################ A function to calculate priors independently. ################
###############################################################################

def evaluate_default_priors(distribution_functions, distribution_parameter_sets, input_values, sampler_type):
    # Just take the given prior distributions and calculate everything independently.
    evaluated_priors = []

    for distribution_function, distribution_parameter_set, input_value in zip(distribution_functions, distribution_parameter_sets, input_values):
        
        prior_distribution = priors[distribution_function]
        prior_distribution.generate_prior(distribution_parameter_set)

        evaluated_prior = prior_distribution.evaluate_prior_from_input(input_value, sampler_type)
        evaluated_priors.append(evaluated_prior)

    return np.asarray(evaluated_priors)
