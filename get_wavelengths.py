from ast import literal_eval as guess_type
from src.ApolloFunctions import SliceModel
from collections import namedtuple
from glob import glob
import numpy as np
import pandas as pd
from warnings import warn

from create_opacity_catalog import OpacityEntry

SpectralParameters = namedtuple(
    "SpectralParameters",
    "minimum_wavelength " +\
    "maximum_wavelength " +\
    "resolution"
    )


def build_list_of_wavelengths(
        list_of_spectral_parameters: list[SpectralParameters]
        ) -> np.array:

    number_of_spectral_elements = \
        lambda parameters: parameters.resolution * \
            np.log(parameters.maximum_wavelength/parameters.minimum_wavelength)

    wavelengths = np.concatenate([
        parameter_set.minimum_wavelength *
        np.exp(
            np.arange(number_of_spectral_elements(parameter_set)) /
            parameter_set.resolution
            )
        for parameter_set in list_of_spectral_parameters
        ])
    return wavelengths


x = build_list_of_wavelengths([SpectralParameters(1.10, 1.34, 6500),
                               SpectralParameters(1.44, 1.85, 8000),
                               SpectralParameters(1.94, 2.56, 10000)])


def construct_model_wavelengths(
        list_of_band_parameters: list[SpectralParameters],
        opacity_entry: OpacityEntry,
        wavelength_calibration_range: list[float],
        model_resolution: float
        ) -> np.array[float]:

    band_lower_bounds = np.array([
        band.minimum_wavelength for band in list_of_band_parameters
        ])
    band_upper_bounds = np.array([
        band.maximum_wavelength for band in list_of_band_parameters
        ])

    if model_resolution < opacity_entry.effective_resolution:
        raise ValueError(
            "Your model resolution can't be higher than the resolution of the "+
            "input opacities (R={}).".format(opacity_entry.effective_resolution)
            )

    opacity_wavelengths = build_list_of_wavelengths(
        opacity_entry.minimum_wavelength,
        opacity_entry.maximum_wavelength,
        model_resolution
        )

    band_indices, model_indices, model_wavelengths = SliceModel(
        band_lower_bounds,
        band_upper_bounds,
        opacity_wavelengths,
        *wavelength_calibration_range
        )

    return band_indices, model_indices, model_wavelengths
