from dataclasses import dataclass


@dataclass
class VerticalProfile:
    pressure: list[float]
    altitude: list[float]
    temperature: list[float]
