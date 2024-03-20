import numpy as np
from scipy.interpolate import PchipInterpolator as monotonic_interpolation
from scipy.ndimage import gaussian_filter1d as gaussian_smoothing


########## ANJALI PIETTE et. al. Profile ##########
# A specific interpolation with its own smoothing. Should be flexible for
# most retrieval cases without thermal inversions. See for reference
# Piette, Anjali A. A., and Nikku Madhusudhan. “Considerations for Atmospheric
# Retrieval of High-Precision Brown Dwarf Spectra” 19, no. July (July 29, 2020):
# 1–19. http://arxiv.org/abs/2007.15004.
def piette(
    T_m4: float,
    T_m3: float,
    T_m2: float,
    T_m1: float,
    T_0: float,
    T_0p5: float,
    T_1: float,
    T_1p5: float,
    T_2: float,
    T_2p5: float,
    pressures: list[float],
):
    logP_nodes = np.array([-4, -3, -2, -1, 0, 0.5, 1, 1.5, 2])
    T_nodes = np.array([T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2])

    interpolated_function = monotonic_interpolation(logP_nodes, T_nodes)

    TP_profile = gaussian_smoothing(interpolated_function(pressures), sigma=0.3)

    return TP_profile


###############################################################################
