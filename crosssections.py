from dataclasses import dataclass
import numpy as np
from os import PathLike
from pathlib import Path
from typing import Protocol
from warnings import warn
from xarray import DataArray, Dataset


@dataclass(slots=True)
class CrossSection(Protocol):
    """
    Defines the header entries and properties that are needed for a collection
    of cross-section tables. Attributes follow the APOLLO-style header.
    """

    number_of_pressure_layers: int
    minimum_log_pressure: float
    maximum_log_pressure: float
    number_of_temperatures: int
    minimum_log_temperature: float
    maximum_log_temperature: float
    number_of_spectral_elements: int
    minimum_wavelength: float
    maximum_wavelength: float
    effective_resolution: float

    @property
    def pressures(self) -> None: ...

    @property
    def temperatures(self) -> None: ...

    @property
    def wavelengths(self) -> None: ...


@dataclass
class CrossSectionDirectory(CrossSection):
    name: str
    species: tuple[str]
    filepaths: dict[str, PathLike | str]

    @property
    def pressures(self):
        return np.linspace(
            self.minimum_log_pressure,
            self.maximum_log_pressure,
            num=self.number_of_pressure_layers,
        )

    @property
    def temperatures(self):
        return np.logspace(
            self.minimum_log_temperature,
            self.maximum_log_temperature,
            num=self.number_of_temperatures,
        )

    @property
    def wavelengths(self):
        return self.minimum_wavelength * np.exp(
            np.arange(self.number_of_spectral_elements) / self.effective_resolution
        )


def get_file_headers(filepaths):
    file_headers = []
    for i, filepath in enumerate(filepaths):
        with open(filepath, "r") as file:
            header = tuple(
                [
                    float(entry) if "." in entry else int(entry)
                    for entry in file.readline().split()
                ]
            )
        file_headers.append(header)

    return tuple(file_headers)


def get_valid_files_and_header(filepaths, header_protocol=CrossSection):
    file_headers = get_file_headers(filepaths)

    expected_number_of_header_entries = len(header_protocol.__slots__)

    valid_filepaths = []
    potentially_valid_headers = []
    for filepath, file_header in zip(filepaths, file_headers):
        if len(file_header) != expected_number_of_header_entries:
            warn(
                f"The file {filepath} does not have "
                + "a header parsable as a gas file. "
                + "Excluding from gas cross section catalog."
            )
            continue
        valid_filepaths.append(filepath)
        potentially_valid_headers.append(file_header)

    file_headers_all_match = len(set(potentially_valid_headers)) == 1
    if not file_headers_all_match:
        warn(
            "There is apparently more than one set of header info. "
            + "Check to see all have the same pressure, temperature, "
            + "and wavelength ranges, and/or resolutions."
        )

    directory_header = potentially_valid_headers[0] if potentially_valid_headers else []

    return valid_filepaths, directory_header


def create_crosssection_catalog(
    file_directory: PathLike | str, filename_match_string: str = r"*.*.dat"
) -> dict[str, CrossSection]:
    crosssection_filepaths = [
        filepath
        for filepath in Path(file_directory).glob(filename_match_string)
        if filepath.is_file()
    ]
    crosssection_filenames = [filepath.name for filepath in crosssection_filepaths]
    species, crosssection_names, _ = list(
        zip(
            *[
                crosssection_filename.split(".")
                for crosssection_filename in crosssection_filenames
            ]
        )
    )

    directory_names = set(crosssection_names)
    directories = {}
    for directory_name in directory_names:
        candidate_filepaths = tuple(
            filepath
            for name, filepath in zip(crosssection_names, crosssection_filepaths)
            if name == directory_name
        )
        directory_species = tuple(
            molecule
            for name, molecule in zip(crosssection_names, species)
            if name == directory_name
        )

        directory_filepaths, directory_header = get_valid_files_and_header(
            candidate_filepaths
        )
        if not directory_header:
            warn(f"No files under the name {directory_name} have valid headers.")
            continue

        directories[directory_name] = CrossSectionDirectory(
            *directory_header,
            directory_name,
            directory_species,
            {
                species: filepath
                for species, filepath in zip(directory_species, directory_filepaths)
            },
        )

    return directories


def load_crosssections_into_array(
    directory: CrossSectionDirectory, excluded_species: list[str] = None
):
    def load_array(filepath, loadtxt_kwargs=dict(skiprows=1)):
        return np.flip(np.loadtxt(filepath, **loadtxt_kwargs)).reshape(
            directory.number_of_pressure_layers,
            directory.number_of_temperatures,
            directory.number_of_spectral_elements,
        )

    return {
        species: load_array(filepath)
        for species, filepath in directory.filepaths.items()
        if species not in excluded_species
    }


def load_crosssections_into_dataset(
    directory: CrossSectionDirectory, excluded_species: list[str] = None
):
    if excluded_species is None:
        excluded_species = []

    crosssection_array_dict = load_crosssections_into_array(directory, excluded_species)

    DataArray_dict = {
        species: DataArray(
            crosssection_array,
            coords=dict(
                pressures=directory.pressures,
                temperatures=directory.temperatures,
                wavelengths=directory.wavelengths,
            ),
        )
        for species, crosssection_array in crosssection_array_dict.items()
    }

    return Dataset(DataArray_dict)
