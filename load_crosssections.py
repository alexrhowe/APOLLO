import numpy as np
from pathlib import Path
import polars as pl
from pprint import pprint
from xarray import Dataset

from crosssections import (
    create_crosssection_catalog,
    load_crosssections_into_dataset,
)

OPACITY_DIRECTORY = Path("/Volumes/ResearchStorage/Opacities_0v10")
GAS_OPACITY_DIRECTORY = OPACITY_DIRECTORY / "gases"


def load_gas_opacity_data_into_dataset(
    gas_opacity_directory: Path,
    opacity_name: str,
    wavelength_units="micron",
    pressure_units="dex(mbar)",
    temperature_units="K",
) -> Dataset:
    crosssection_catalog = create_crosssection_catalog(gas_opacity_directory)

    crosssection_dataset = load_crosssections_into_dataset(
        crosssection_catalog[opacity_name]
    )

    crosssection_dataset.wavelengths.attrs.update(dict(units=wavelength_units))
    crosssection_dataset.pressures.attrs.update(dict(units=pressure_units))
    crosssection_dataset.temperatures.attrs.update(dict(units=temperature_units))

    return crosssection_dataset


opacities_interpd = opacity_dataset.interp(
    pressures=Ps, temperatures=Ts, wavelengths=waves, assume_sorted=True
)
