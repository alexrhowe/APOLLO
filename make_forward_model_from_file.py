import numpy as np
import sys

sys.path.append("../newAPOLLO")

from apollo.Apollo_components import (
    ReadInputsfromFile,
    ProcessInputs,
    MakeObservation,
    MakePlanet,
)

OPACITY_DIRECTORY = "/Volumes/ResearchStorage/Opacities_0v10/"


def generate_model_spectrum(input_filepath):
    inputs = ReadInputsfromFile(input_filepath)
    inputs["opacdir"] = OPACITY_DIRECTORY

    processed_inputs = ProcessInputs(**inputs)

    planet = MakePlanet(processed_inputs["MakePlanet_kwargs"])

    get_model = planet.MakeModel(processed_inputs["MakeModel_initialization_kwargs"])

    exact_parameters = processed_inputs["parameters"]

    model = get_model(exact_parameters)

    observation = MakeObservation(
        processed_inputs["ModelObservable_initialization_kwargs"]
    )

    obs_args = [
        np.asarray(model[0]),
        processed_inputs["parameters"][0],
        "Rad",
        exact_parameters[-4],
    ]

    observed_binned_model, observed_full_resolution_model = observation(*obs_args)

    return observed_binned_model["data"]
