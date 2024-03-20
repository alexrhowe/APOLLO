from astropy import units as u
from functools import partial
from inspect import signature
import numpy as np
import os
from pathlib import Path
from picaso import justdoit as jdi
from typing import Callable
import xarray as xr
from xarray import DataArray

FILENAME_DB = os.path.join(
    os.environ.get("picaso_refdata"), "opacities", "opacities.db"
)
OPANNECTION = jdi.opannection(filename_db=FILENAME_DB)
SONORA_PROFILE_DB = "/Volumes/ResearchStorage/profile"
SAVE_DIRECTORY = Path.cwd() / "extensions"


def make_picaso_brown_dwarf():
    brown_dwarf = jdi.inputs(calculation="browndwarf")
    brown_dwarf.phase_angle(0)

    return brown_dwarf


def add_SONORA_TP_profile(
    picaso_brown_dwarf: jdi.inputs, log_g: float, effective_temperature: float
):
    picaso_brown_dwarf.gravity(gravity=10**log_g, gravity_unit=u.Unit("cm/(s**2)"))
    picaso_brown_dwarf.sonora(SONORA_PROFILE_DB, effective_temperature)
    picaso_brown_dwarf.spectrum(OPANNECTION, full_output=True)

    return picaso_brown_dwarf


def get_SONORA_TP_profile(
    picaso_brown_dwarf: jdi.inputs, log_g: float, effective_temperature: float
):
    SONORA_brown_dwarf = add_SONORA_TP_profile(picaso_brown_dwarf)

    log_pressures = np.asarray(
        SONORA_brown_dwarf.inputs["atmosphere"]["profile"]["pressure"].apply(np.log10)
    )
    temperatures = np.asarray(
        SONORA_brown_dwarf.inputs["atmosphere"]["profile"]["temperature"]
    )

    return dict(log_pressures=log_pressures, temperatures=temperatures)


def get_chemical_equilibrium_profile(
    picaso_brown_dwarf: jdi.inputs,
    c_to_o_relative_to_solar: float = 1,
    metallicity: float = 0,
    molecules=None,
):
    picaso_brown_dwarf.chemeq_visscher(c_to_o_relative_to_solar, metallicity)

    chemeq_dataframe = picaso_brown_dwarf.inputs["atmosphere"]["profile"]
    molecule_columns = chemeq_dataframe.columns.difference(["temperature"])
    chemeq_dataframe[molecule_columns] = chemeq_dataframe[molecule_columns].applymap(
        np.log10
    )

    chemeq_dataarray = (
        chemeq_dataframe.set_index(["pressure"]).to_xarray().set_coords("temperature")
    )
    return chemeq_dataarray if not molecules else chemeq_dataarray[molecules]


def find_photospheric_pressure(
    effective_temperature: float, pressure: list[float], temperature: list[float]
):
    index_upper_limit = np.nonzero(temperature >= effective_temperature)[0][0]

    if index_upper_limit > 0:
        index_lower_limit = index_upper_limit - 1

        temperature_lower_limit = temperature[index_lower_limit]
        temperature_upper_limit = temperature[index_upper_limit]
        fractional_index = (effective_temperature - temperature_lower_limit) / (
            temperature_upper_limit - temperature_lower_limit
        )

        pressure_lower_limit = pressure[index_lower_limit]
        pressure_upper_limit = pressure[index_upper_limit]
        photospheric_pressure = pressure_lower_limit + fractional_index * (
            pressure_upper_limit - pressure_lower_limit
        )

        return photospheric_pressure

    else:
        return pressure[0]


def provide_parameters_from_dataset(function: Callable, dataset: DataArray):

    def evaluate_from_dataset(function: Callable, dataset: DataArray, *args, **kwargs):
        dataset_parameters = {
            parameter.name: dataset.get(parameter.name).to_numpy()
            for parameter in signature(function).parameters.values()
            if parameter.name in list(dataset.keys())
            or parameter.name in list(dataset.coords)
        }
        return function(*args, **dataset_parameters, **kwargs)

    return partial(evaluate_from_dataset, function, dataset)


def estimate_significant_molecules(
    log_g: float,
    effective_temperature: float,
    ctoo_vs_solar: float,
    metallicity: float,
    significance_threshold=-8,
):
    picaso_model = add_SONORA_TP_profile(
        make_picaso_brown_dwarf(), log_g, effective_temperature
    )
    chemeq_profile = get_chemical_equilibrium_profile(
        picaso_model, ctoo_vs_solar, metallicity
    )

    get_photospheric_pressure_from_dataset = provide_parameters_from_dataset(
        find_photospheric_pressure, chemeq_profile
    )
    photospheric_pressure = get_photospheric_pressure_from_dataset(
        effective_temperature
    )

    photospheric_chemistry = chemeq_profile.interp(pressure=photospheric_pressure)
    significant_molecules = [
        molecule
        for molecule in photospheric_chemistry.keys()
        if (photospheric_chemistry[molecule] >= significance_threshold).all()
    ]

    # return photospheric_chemistry[significant_molecules]
    return chemeq_profile.get(significant_molecules)


LOG_G = 5.44
EFFECTIVE_TEMPERATURE = 1200
CTOO_VS_SOLAR = 0.80 / 0.48
METALLICITY = 0

save_filename = (
    "SONORA_chemical-equilibrium_profile_"
    + f"logg{LOG_G}_"
    + f"Teff{EFFECTIVE_TEMPERATURE}_"
    + f"ctoo-vs-solar{CTOO_VS_SOLAR}_"
    + f"metallicity{METALLICITY}"
    + ".nc"
)
save_filepath = SAVE_DIRECTORY / Path(save_filename)

significant_molecules = estimate_significant_molecules(
    LOG_G, EFFECTIVE_TEMPERATURE, CTOO_VS_SOLAR, METALLICITY
)
print(significant_molecules)

significant_molecules.to_netcdf(save_filepath)
