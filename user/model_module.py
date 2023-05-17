from user.priors import evaluate_default_priors


###############################################################################
############# A class to load model functions and special priors. #############
###############################################################################

class ModelModule:
    def __init__(self, model_function, prior_function=evaluate_default_priors):
        self.model_function = model_function
        self.prior_function = prior_function


    def evaluate_model(self, *model_args, **model_kwargs):
        return self.model_function(*model_args, **model_kwargs)


    def evaluate_prior(self, distribution_functions, distribution_parameter_sets, input_values, sampler_type):
        return self.prior_function(distribution_functions, distribution_parameter_sets, input_values, sampler_type)

###############################################################################
