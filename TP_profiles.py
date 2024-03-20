import numpy as np
from user.defaults import minP, maxP, minT, maxT, tgray, vres
from scipy.interpolate import PchipInterpolator as monotonic_interpolation
from scipy.ndimage import gaussian_filter1d as gaussian_smoothing


# def gray(T=tgray, num_layers_final=vres, P_min=minP, P_max=maxP):
#     return np.ones(vres)*T


# def verbatim(*input_Ts, num_layers_final=vres, P_min=minP, P_max=maxP):
#     return layers

'''
def interpolate(*input_Ts, num_layers_final=vres, P_min=minP, P_max=maxP):
    num_layers_initial = len(input_Ts)
    input_Ps = np.linspace(P_max, P_min, num=num_layers_initial)
    output_Ps = np.linspace(P_max, P_min, num=num_layers_final)

    if num_layers <= 4:
        interpolation_level = linear
    else:
        interpolation_level = cubic

    output_Ts = (np.interp1d(input_Ps, input_Ts, kind=interpolation_level))(output_Ps)
    output_Ts = np.where(output_Ts < minT, minT, output_Ts)
    output_Ts = np.where(output_Ts > maxT, maxT, output_Ts)

    return output_Ts
'''

def power_linear(exponent, logPiso, T0, T_mid, T_max,
                 num_layers_final=vres, P_min=minP, P_max=maxP):
    logP0 = P_min
    T_frac = exponent * (T_mid-T0)/(T_max-T_mid)
    logPmid = (logP0+T_frac*logPiso) / (1+T_frac)
    alpha2 = (logPiso-logPmid) / (T_max-T_mid)
    alpha1 = exponent * (T_mid-T0)**(1-1/exponent) * alpha2

    input_Ps = np.linspace(P_min, P_max, num=num_layers_final)

    output_Ts = np.ones_like(input_Ps) * T0
    output_Ts = np.where((input_Ps>logP0) & (input_Ps<logPmid),
                         T0 + ((input_Ps-logP0)/alpha1)**exponent,
                         output_Ts)
    output_Ts = np.where((input_Ps>logPmid) & (input_Ps<logPiso),
                         T_mid + ((input_Ps-logPmid)/alpha2),
                         output_Ts)
    output_Ts = np.where(input_Ps >= logPiso,
                         T_max,
                         output_Ts)

    return output_Ts


def Piette(T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2, T_2p5,
           num_layers_final=vres, P_min=minP, P_max=maxP):
    logP_nodes = np.array([-4, -3, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5])
    T_nodes = np.array([T_m4, T_m3, T_m2, T_m1, T_0, T_0p5, T_1, T_1p5, T_2, T_2p5])

    interpolated_function = monotonic_interpolation(logP_nodes, T_nodes)
    input_Ps = np.linspace(P_min, P_max, num=num_layers_final)

    output_Ts = gaussian_smoothing(interpolated_function(input_Ps), sigma=0.3)

    return output_Ts
