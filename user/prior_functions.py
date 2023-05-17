import numpy as np
from scipy.stats import norm, truncnorm, uniform

'''
class prior:
    def __init__(self, prior_function, prior_parameters):
        self.prior = prior_function(prior_parameters)


    def generate_prior(self, prior_function):
        self.prior = prior_function


    def sample_prior(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the PDF for its prior calculation.
            return self.prior.logpdf(input_value)

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior from a uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")


def uniform_prior_function(prior_parameters):
    _, _, uniform_min, uniform_max = prior_parameters

    return uniform(uniform_min, uniform_max)


def normal_prior_function(prior_parameters):
    mean, sigma, _, _ = prior_parameters

    return norm(loc=mean, scale=sigma)


def truncated_normal_prior(prior_parameters):
    mean, sigma, lower_bound, upper_bound = prior_parameters
    a, b = (lower_bound - mean) / sigma, (upper_bound - mean) / sigma

    return truncnorm(a=a, b=b, loc=mean, scale=sigma)
'''


class uniform_prior:
    def __init__(self, prior_parameters):
        _, _, uniform_min, uniform_max = prior_parameters
        uniform_width = uniform_max - uniform_min
        self.prior = uniform(uniform_min, uniform_width)

    
    def sample_prior(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the PDF for its prior calculation.
            return self.prior.logpdf(input_value)

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior from a uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")


class normal_prior:
    def __init__(self, prior_parameters):
        mean, sigma, _, _ = prior_parameters
        self.prior = norm(loc=mean, scale=sigma)

    
    def sample_prior(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the PDF for its prior calculation.
            return self.prior.logpdf(input_value)

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior from a uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")


class truncated_normal_prior:
    def __init__(self, prior_parameters):
        mean, sigma, lower_bound, upper_bound = prior_parameters
        a, b = (lower_bound - mean) / sigma, (upper_bound - mean) / sigma
        self.prior = truncnorm(a=a, b=b, loc=mean, scale=sigma)

    
    def sample_prior(self, input_value, sampler_type):
        if sampler_type == "emcee":

            # emcee uses the PDF for its prior calculation.
            return self.prior.logpdf(input_value)

        elif sampler_type == "dynesty":

            # dynesty uses the inverse CDF (ppf) to sample the prior from a uniform sample on [0, 1].
            return self.prior.ppf(input_value)

        else:
            raise ValueError("Sampler type not valid.")


prior_dictionary = {

    "Uniform": uniform_prior,

    "Normal": normal_prior,

    "Truncated_Normal": truncated_normal_prior

}
