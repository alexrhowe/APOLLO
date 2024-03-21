from dataclasses import dataclass
from typing import Protocol


@dataclass
class RadiativeProfile:
    absorption_coefficients: list[float]
    scattering_coefficients: list[
        list[float]
    ]  # for now, 2 elements: forward and backward


class Absorbs(Protocol):
    """Thing that absorbs. Profile is equivalent to dTau/dP."""

    def calculate_absorption_profile(self, *args, **kwargs) -> None: ...


class Scatters(Protocol):
    """Thing that scatters, possibly backwards. Profile is equivalent to dTau/dP."""

    def calculate_forward_scattering_profile(self, *args, **kwargs) -> None: ...

    def calculate_backward_scattering_profile(self, *args, **kwargs) -> None: ...


class Material(Absorbs, Scatters, Protocol):
    """Placeholder for a thing that both absorbs and scatters."""

    pass


class UsesCrossSections(Material, Protocol):
    """A material that is built with cross-section data."""

    def calculate_number_density(self, *args, **kwargs) -> None: ...

    def calculate_absorption_cross_section(self, *args, **kwargs) -> None: ...

    def calculate_forward_scattering_cross_section(self, *args, **kwargs) -> None: ...

    def calculate_backward_scattering_cross_section(self, *args, **kwargs) -> None: ...
