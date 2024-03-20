from dataclasses import dataclass
from functools import partial
import numpy as np
from numpy.typing import NDArray
import polars as pl
import polars.selectors as cs
from typing import Callable, NamedTuple

from general_protocols import Spectral, Measured
from src.ApolloFunctions import ConvSpec, BinSpec


@dataclass
class InputFileReader:
    """Docstring goes here."""

    pass


@dataclass
class ModelBuilder:
    """Docstring goes here."""

    """
    You could have
    - a core function, that translates a set of values and some core component of the
        physical structure (pressure, wavelength) that are arguments.
        You return some physical representation of a part of the model.

    - a method that takes the parameters and partially applies them to
        the function, leaving a function of just pressure, for example.

    - a container for plotting information, like the long name, units, etc.

    - flexibility to calculate 
    """
    model_function: Callable
    model_attributes: NamedTuple

    def __call__(self, *args: np.Any, **kwargs: np.Any) -> np.Any:
        return self.model_function(*args, **kwargs)

    def load_function(self, *parameter_args, **parameter_kwargs):
        return partial(self.model_function, *parameter_args, **parameter_kwargs)


def read_APOLLO_data(
    filepath,
    band_names=None,
    data_column_names=[
        "wavelength_bin_starts",
        "wavelength_bin_ends",
        "spectrum",
        "lower_errors",
        "upper_errors",
        "spectrum_with_noise",
    ],
) -> dict[str, DataSpectrum]:
    return read_Polars_to_Spectrum(
        read_APOLLO_data_into_Polars(filepath, band_names, data_column_names)
    )


def read_APOLLO_data_into_Polars(
    filepath,
    band_names=None,
    data_column_names=[
        "wavelength_bin_starts",
        "wavelength_bin_ends",
        "spectrum",
        "lower_errors",
        "upper_errors",
        "spectrum_with_noise",
    ],
) -> dict[str, pl.DataFrame]:
    data = pl.read_csv(filepath, separator=" ", has_header=False).drop(~cs.numeric())

    column_renaming_dict = {
        default_column_name: data_column_name
        for default_column_name, data_column_name in zip(
            data.columns, data_column_names
        )
    }

    data = data.rename(column_renaming_dict)

    band_start_indices = (
        (
            data.select(
                [
                    pl.arg_where(
                        ~pl.col("wavelength_bin_ends").is_in(
                            pl.col("wavelength_bin_starts")
                        )
                    )
                ]
            )
            + 1
        )
        .to_numpy()
        .squeeze()
    )
    band_start_indices = np.hstack((0, band_start_indices))
    band_slices = [
        slice(start, end)
        for start, end in zip(band_start_indices[:-1], band_start_indices[1:])
    ]
    number_of_bands = len(band_slices)

    if band_names is None:
        band_names = [f"Band_{n+1}" for n in range(number_of_bands)]

    banded_data = {
        band_name: data[band_slice]
        for band_name, band_slice in zip(band_names, band_slices)
    }
    return banded_data


def read_Polars_to_Spectrum(banded_data: dict[str, pl.DataFrame]):
    return {
        band_name: DataSpectrum.from_binned_data(
            **{
                attribute: np.asarray(value.to_numpy())
                for attribute, value in data.to_dict().items()
                if attribute in DataSpectrum.__slots__
            }
        )
        for band_name, data in banded_data.items()
    }
