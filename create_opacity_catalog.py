from ast import literal_eval as guess_type
# from collections import namedtuple
from dataclasses import dataclass, field, InitVar
from glob import glob
import pandas as pd
from typing import NamedTuple
from warnings import warn


OpacityEntry = namedtuple(
    'OpacityEntry',
    'wavelength_range_name '   +\
    'species '                 +\
    'number_of_pressures '     +\
    'minimum_log_pressure '    +\
    'maximum_log_pressure '    +\
    'number_of_temperatures '  +\
    'minimum_log_temperature ' +\
    'maximum_log_temperature ' +\
    'table_size '              +\
    'minimum_wavelength '      +\
    'maximum_wavelength '      +\
    'effective_resolution '    +\
    'path_to_file'
    )

@dataclass
class OpacityCatalog:
    name:                        str
    species:                     list[str]
    number_of_pressure_layers:   int
    minimum_log_pressure:        float
    maximum_log_pressure:        float
    number_of_temperatures:      int
    minimum_log_temperature:     float
    maximum_log_temperature:     float
    number_of_spectral_elements: int
    minimum_wavelength:          float
    maximum_wavelength:          float
    effective_resolution:        float
    paths_to_files: 


def create_gas_opacity_catalog(gas_opacity_directory: str) -> pd.DataFrame:
    opacity_filepaths = glob(gas_opacity_directory+'/*.dat')

    expected_number_of_header_entries = len(OpacityEntry._fields)

    catalog_lines = []
    for i, filepath in enumerate(opacity_filepaths):
        with open(filepath, 'r') as opacity_file:
            header_line = opacity_file.readline().split(' ')

            raw_file_ids = filepath[(filepath.index(gas_opacity_directory) +
                                len(gas_opacity_directory)+1):].split('.')
            entry_ids = [raw_file_ids[-2], raw_file_ids[0]]

            if len(header_line)+len(entry_ids)+1 \
                != expected_number_of_header_entries:
                warn('The {} {} file'.format(*entry_ids) + ' does not have ' +
                     'a header parsable as a gas file. ' +
                     'Excluding from gas opacity catalog.')
                continue
            else:
                catalog_line = OpacityEntry(
                    *entry_ids,
                    *[guess_type(entry)for entry in header_line],
                    filepath
                    )
                catalog_lines.append(catalog_line)

    opacity_catalog = pd.DataFrame(catalog_lines).set_index(
        ['wavelength_range_name', 'species']
        )
    return opacity_catalog


def get_opacities_from_catalog(
        gas_opacity_directory:         str,
        opacity_wavelength_range_name: str,
        gas_species:           list[str]
        ) -> dict:

    opacity_catalog = create_gas_opacity_catalog(gas_opacity_directory)
    if opacity_catalog.empty:
        warn('Opacity catalog is empty! You might not have specified ' +
             'the correct directory, or the file names and/or headers ' +
             'might not be in the correct format.')

    opacities = opacity_catalog.\
        loc[opacity_wavelength_range_name].\
            loc[gas_species]

    return opacities


def load_test_catalog() -> pd.DataFrame:
    test_catalog = get_opacities_from_catalog(
        '/Volumes/Research Mobile/Opacities_0v10/gases',
        'nir',
        ['h2o', 'co', 'co2']
        )

    return(test_catalog)
