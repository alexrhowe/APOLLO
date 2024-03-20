from dataclasses import dataclass, field, InitVar
from numpy.typing import NDArray
import numpy as np
from pathlib import Path
from pint import Quantity, UnitRegistry
from pprint import pprint
from typing import Any, Callable, IO, Optional, Protocol, Sequence

Pathlike = Path | str
Filelike = Pathlike | IO

LIST_OF_MEASURED_VARIABLES = [
    "lower_errors",
    "upper_errors",
]


class Measured(Protocol):
    """A thing that has errors.

    Args:
        Protocol (_type_): _description_

    Returns:
        _type_: _description_
    """

    lower_errors: NDArray
    upper_errors: NDArray

    @property
    def errors(self):
        return (self.lower_errors + self.upper_errors) / 2


class Parametrized(Protocol):
    """
    A sub-model for some component (for example, T-P profile, clouds) that
    has a function that takes parameters and returns some physical ``profile''
    function, i.e. a temperature function that takes pressures directly.
    """

    model_function: Callable

    def generate_profile_function(self, *args, **kwargs) -> None: ...


@dataclass
class Parameter:
    name: str
    print_name: str
    value: InitVar[float | int]
    unit: InitVar[str]
    print_formatter: InitVar[str]
    _quantity: Quantity = field(init=False)

    def __post_init__(self, value, unit, print_formatter):
        default_unit_registry = UnitRegistry()
        self._quantity = value * default_unit_registry(unit)
        self._quantity.default_format = f"{print_formatter}P"

    def __getattr__(self, __name: str) -> Any:
        return self._quantity.__getattribute__(__name)

    def __repr__(self):
        return f"{self.name}: {self._quantity}"


test_parameter = Parameter(
    name="radius",  # must match the argument name in the function
    print_name="R",
    value=np.array([2.33485734, 4.29457829, 6.48438237]),
    unit="meters",
    print_formatter=".2f",
)

pprint(test_parameter)


@dataclass
class ThermalModel(Parametrized):
    """Returns a profile of T(P) when called."""

    pass


@dataclass
class VerticalModel(Parametrized):
    """Returns a profile of z(P) when called."""

    pass


@dataclass
class CompositionModel(Parametrized):
    """
    Holds models that can return absorption and scattering coefficients
    as functions of T, P or z. And quantities of specific elements???
    I think you want a general object that can give you dTau/dP, and
    then maybe an extension for objects with chemical compositions.
    """

    pass


@dataclass
class AtmosphericModel(Parametrized):
    """
    Returns a profile of optical parameters
    (optical depth: tau, scattering asymmetry: g, single-scattering albedo: omega_0)
    as functions of P/z.
    """

    thermal_model: ThermalModel
    vertical_model: VerticalModel
    composition_models: dict[str, CompositionModel]
