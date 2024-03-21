from enum import StrEnum, auto


class ObservationMode(StrEnum):
    RESOLVED = auto()
    ECLIPSE = auto()
    TRANSIT = auto()
