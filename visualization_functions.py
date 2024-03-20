import numpy as np
from colormath.color_objects import sRGBColor, HSLColor
from colormath.color_conversions import convert_color
from matplotlib.colors import LinearSegmentedColormap, to_rgb


def create_linear_colormap(color_arg,
                           resolution=256,
                           lightness_minimum=0.05,
                           lightness_maximum=0.95,
                           saturation_minimum=0.2,
                           saturation_maximum=0.8):
    base_colors = np.atleast_1d(color_arg)
    base_rgbs = [to_rgb(color) for color in base_colors]
    lightnesses = np.linspace(lightness_maximum, lightness_minimum, num=resolution)
    saturations = np.linspace(saturation_minimum, saturation_maximum, num=resolution)

    if len(base_colors) > 1:
        unscaled_colormap = LinearSegmentedColormap.from_list("base_cmap", base_rgbs, N=resolution)
        unscaled_colors = unscaled_colormap(np.linspace(0, 1, num=resolution))
        unscaled_HSLs = [convert_color(sRGBColor(*color), HSLColor)
                         for color in unscaled_colors]
        scaled_HSLs = [HSLColor(unscaled_HSL.hsl_h, saturation, lightness)
                       for unscaled_HSL, lightness, saturation
                       in zip(unscaled_HSLs, lightnesses, saturations)]

    else:
        unscaled_color = to_rgb(color_arg)
        unscaled_HSL = convert_color(sRGBColor(*unscaled_color), HSLColor)
        scaled_HSLs = [HSLColor(unscaled_HSL.hsl_h, saturation, lightness)
                       for lightness, saturation in zip(lightnesses, saturations)]

    scaled_colors = [convert_color(HSL, sRGBColor).get_value_tuple()
                     for HSL in scaled_HSLs]
    scaled_colormap = LinearSegmentedColormap.from_list("scaled_cmap", scaled_colors, N=resolution)

    return scaled_colormap
