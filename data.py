import numpy as np
import polars as pl
import polars.selectors as pls
from spectral import DataSpectrum
from toolz.functoolz import compose_left
from typing import Any, Callable, Sequence
from warnings import warn

from general_protocols import Filelike

APOLLO_COLUMN_NAMES = [
    "wavelength_bin_starts",
    "wavelength_bin_ends",
    "spectrum",
    "lower_errors",
    "upper_errors",
    "spectrum_with_noise",
]
# NOTE TO SELF: look into "suffix" to add a "_with_noise" or whatever else you want.

APOLLO_COLUMN_NAMES_WITHOUT_NOISE = APOLLO_COLUMN_NAMES[:-1]


def read_APOLLO_data(
    filepath,
    band_names: Sequence[str],
    data_column_names: Sequence[str] = APOLLO_COLUMN_NAMES,
    post_polars_processing_functions: Sequence[Callable] | None = None,
) -> dict[str, DataSpectrum]:
    """_summary_

    Args:
        filepath (_type_): _description_
        band_names (list[str], optional): _description_. Defaults to None.
        data_column_names (list[str], optional): _description_. Defaults to APOLLO_COLUMN_NAMES.
        post_polars_processing_functions (list[Callable], optional): _description_. Defaults to None.

    Returns:
        dict[str, DataSpectrum]: _description_
    """
    if post_polars_processing_functions is not None:
        apply_post_polars_functions = compose_left(*post_polars_processing_functions)

        polars_data = apply_post_polars_functions(
            read_APOLLO_data_into_Polars(filepath, data_column_names)
        )

    else:
        polars_data = read_APOLLO_data_into_Polars(filepath, data_column_names)

    polars_data_by_bands = slice_Polars_data_into_bands(
        polars_data,
        band_names,
    )
    return {
        band_name: DataSpectrum.from_wavelength_bins(dataframe=polars_data_by_band)
        for band_name, polars_data_by_band in polars_data_by_bands.items()
    }


def read_APOLLO_data_into_Polars(
    filepath: Filelike,
    data_column_names: Sequence[str] = APOLLO_COLUMN_NAMES,
) -> dict[str, pl.DataFrame]:
    """_summary_

    Args:
        filepath (Filelike): _description_
        data_column_names (Sequence[str], optional): _description_. Defaults to APOLLO_COLUMN_NAMES.

    Returns:
        dict[str, pl.DataFrame]: _description_
    """
    data = pl.read_csv(filepath, separator=" ", has_header=False).drop(~pls.numeric())

    column_renaming_dict = {
        default_column_name: data_column_name
        for default_column_name, data_column_name in zip(
            data.columns, data_column_names
        )
    }

    data = data.rename(column_renaming_dict)

    return data


def find_band_slices_in_Polars(dataframe: pl.DataFrame) -> list[slice]:
    """_summary_

    Args:
        dataframe (pl.DataFrame): _description_

    Returns:
        list[slice]: _description_
    """
    band_start_indices = (
        (
            dataframe.select(
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

    return band_slices


def slice_Polars_data_into_bands(
    dataframe: pl.DataFrame, band_names: Sequence[str]
) -> dict[str, pl.DataFrame]:
    """_summary_

    Args:
        dataframe (pl.DataFrame): _description_
        band_names (Sequence[str]): _description_

    Returns:
        dict[str, pl.DataFrame]: _description_
    """
    band_slices = find_band_slices_in_Polars(dataframe)
    number_of_bands = len(band_slices)

    if band_names is None:
        band_names = [f"Band_{n+1}" for n in range(number_of_bands)]
    else:
        number_of_bands_expected = len(band_names)
        if not number_of_bands == number_of_bands_expected:
            warn(
                f"The number of bands found ({number_of_bands}) "
                + f"does not match the number of band names you provided ({number_of_bands_expected}). "
                + "This could be a number precision issue."
            )

    return {
        band_name: dataframe[band_slice]
        for band_name, band_slice in zip(band_names, band_slices)
    }


def mock_spectrum_with_noise_column(dataframe: pl.DataFrame):
    return dataframe.with_columns(
        spectrum_with_noise=dataframe.select("spectrum").to_numpy().squeeze()
    )


def merge_bands_into_single_dataframe(
    dataframes: Sequence[pl.DataFrame],
    shared_data_entries: Sequence[str] = APOLLO_COLUMN_NAMES_WITHOUT_NOISE,
) -> pl.DataFrame:
    dataframes_with_common_entries = [
        dataframe.select(*shared_data_entries) for dataframe in dataframes
    ]
    return pl.concat(dataframes_with_common_entries)


def write_Polars_data_to_APOLLO_file(
    dataframe: pl.DataFrame,
    file: Filelike,
    write_csv_kwargs: dict[str, Any] = dict(include_header=False, separator=" "),
    pre_save_processing_functions: list[Callable] | None = None,
):
    if pre_save_processing_functions is not None:
        apply_pre_save_functions = compose_left(*pre_save_processing_functions)
        dataframe_for_file = apply_pre_save_functions(dataframe)
    else:
        dataframe_for_file = dataframe

    if not "spectrum_with_noise" in dataframe_for_file.columns:
        dataframe_with_mock_noise = mock_spectrum_with_noise_column(dataframe_for_file)
        return dataframe_with_mock_noise.select(*APOLLO_COLUMN_NAMES).write_csv(
            file, **write_csv_kwargs
        )

    else:
        return dataframe_for_file.select(*APOLLO_COLUMN_NAMES).write_csv(
            file, **write_csv_kwargs
        )


def reduce_float_precision_to_correct_band_finding(
    dataframe: pl.DataFrame, number_of_decimal_places=5
):
    """_summary_

    Args:
        dataframe (pl.DataFrame): _description_
        number_of_decimal_places (int, optional): _description_. Defaults to 5.

    Returns:
        _type_: _description_
    """
    return dataframe.with_columns(
        pl.col("wavelength_bin_starts").round(number_of_decimal_places),
        pl.col("wavelength_bin_ends").round(number_of_decimal_places),
    )
