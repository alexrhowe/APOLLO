from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import cnames
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path
import pickle
from plot_functions import *
from TP_profiles import Piette

from visualization_functions import create_linear_colormap

plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("axes", labelsize=19)
plt.rc("xtick", labelsize=15)
plt.rc("ytick", labelsize=15)

filetypes = ["pdf", "png"]

DIRECTORY_RESULTS_2M2236 = Path.cwd() / "results" / "2M2236"
RESULTS_FILE_2M2236 = (
    DIRECTORY_RESULTS_2M2236
    / "2M2236.Piette.G395H.cloud-free.2024-02-27.continuation.dynesty_results.dat"
)
DERIVED_FILE_2M2236 = (
    DIRECTORY_RESULTS_2M2236
    / "2M2236.Piette.G395H.cloud-free.2024-02-27.continuation.dynesty_derived.dat"
)
OBSERVED_SPECTRUM_FILE_2M2236 = "2M2236b_NIRSpec_G395H_R500_APOLLO.dat"

MLE_SPECTRUM_FILE_2M2236 = "2M2236.Piette.G395H.cloud-free.2024-02-27.continuation.retrieved.Spectrum.binned.dat"
CONTRIBUTIONS_FILE_2M2236 = (
    "2M2236.Piette.G395H.cloud-free.2024-02-27.continuation.retrieved_contributions.p"
)

MODEL_FILENAME_TEMPLATE = "2M2236.Piette.G395H.cloud-free.2024-02-27.continuation"
OUTPUT_FILENAME_TEMPLATE = "2M2236.G395H.cloud-free.2024-02-27.continuation"

with open(RESULTS_FILE_2M2236, "rb") as results_file:
    dynesty_results = pickle.load(results_file)

results = dynesty_results.samples
log_likelihoods = dynesty_results.logl
importance_weights = np.exp(
    (log_importance := dynesty_results.logwt) - np.max(log_importance)
)
cumsum_logwt = np.cumsum(importance_weights)

# combined retrievals on 2M2236
directories = ["2M2236"]
object_label = "2M2236.G395H.cloud-free.2024-01-22"
object_headers = ["2M2236.Piette.G395H.cloud-free.2024-01-22"]
data_files = ["2M2236b_NIRSpec_G395H_R500_APOLLO.dat"]
plotting_colors = ["lightcoral"]

MLE_color = "gold"

cloudy_parameter_names = [
    "Rad",
    "Log(g)",
    "h2o",
    "co",
    "co2",
    "ch4",
    "Lupu_alk",
    "h2s",
    "nh3",
    "T_m4",
    "T_m3",
    "T_m2",
    "T_m1",
    "T_0",
    "T_0p5",
    "T_1",
    "T_1p5",
    "T_2",
    "T_2p5",
    "deltaL",
    "scaleG395_ch1",
    "scaleG395_ch2",
    "logf",
    "Mass",
    "C/O",
    "[Fe/H]",
    "Teff",
]

temp_index_bounds = (9, 18)
gas_index_bounds = (2, 9)
# cloud_index_bounds = (19, 24)

temps = slice(*temp_index_bounds)
gas = slice(*gas_index_bounds)
gases = cloudy_parameter_names[gas]


def calculate_CtoO_and_metallicity(gas_logabundances):
    carbon = 0.0
    oxygen = 0.0
    metals = 0.0
    ccompounds = ["ch4", "co", "co2", "hcn"]
    cmult = [1.0, 1.0, 1.0, 1.0]
    ocompounds = ["h2o", "co", "co2", "tio", "vo"]
    omult = [1.0, 1.0, 2.0, 1.0, 1.0]
    zcompounds = [
        "h2o",
        "ch4",
        "co",
        "co2",
        "nh3",
        "h2s",
        "Burrows_alk",
        "Lupu_alk",
        "crh",
        "feh",
        "tio",
        "vo",
        "hcn",
        "n2",
        "ph3",
    ]
    zmult = [
        16.0,
        12.0,
        28.0,
        44.0,
        14.0,
        32.0,
        24.0,
        24.0,
        52.0,
        56.0,
        64.0,
        67.0,
        26.0,
        28.0,
        31.0,
    ]

    for i in np.arange(0, len(ccompounds)):
        if ccompounds[i] in gases:
            j = gases.index(ccompounds[i])
            carbon = carbon + cmult[i] * (
                10 ** gas_logabundances[j]
            )  # -1 because of hydrogen
    for i in np.arange(0, len(ocompounds)):
        if ocompounds[i] in gases:
            j = gases.index(ocompounds[i])
            oxygen = oxygen + omult[i] * (10 ** gas_logabundances[j])
    for i in np.arange(0, len(zcompounds)):
        if zcompounds[i] in gases:
            j = gases.index(zcompounds[i])
            metals = metals + zmult[i] * (10 ** gas_logabundances[j])

    ctoo = carbon / oxygen
    fetoh = np.log10(metals / 0.0196)

    return np.array([ctoo, fetoh])


sample_dicts = {}
results_dicts = {}
derived_dicts = {}
importance_dicts = {}
iteration_dicts = {}
string_format_dicts = {}

for directory, object_header, parameter_names in zip(
    directories, object_headers, [cloudy_parameter_names]
):
    results_file = RESULTS_FILE_2M2236
    with open(results_file, "rb") as file:
        results = pickle.load(file)
    results_dicts[directory] = results
    print(results.niter)

    derived_file = DERIVED_FILE_2M2236
    with open(derived_file, "rb") as file:
        derived = np.asarray(pickle.load(file))
    derived_dicts[directory] = derived

    important_percentile = 0.95
    important_iterations = cumsum_logwt / cumsum_logwt[-1] >= (1 - important_percentile)
    iteration_dicts[directory] = important_iterations
    samples = results.samples[important_iterations]
    samples = np.concatenate(
        (results.samples[important_iterations], derived[important_iterations]), axis=-1
    )
    likelihoods = log_likelihoods[important_iterations]
    importance_dicts[directory] = importance_weights[important_iterations]

    # Correct for radius offset.
    """
    scale_J, scale_H, scale_K = [-4, -3, -2]
    scale_factor_K = samples[:, scale_K]
    samples[:, scale_J] /= scale_factor_K
    samples[:, scale_H] /= scale_factor_K
    samples[:, scale_K] /= scale_factor_K
    samples[:, 0] *= np.sqrt(scale_factor_K)
    # samples[:, -4] *= scale_factor_K**1.5
    """

    samples[:, 0] /= 11.2

    # samples = np.append(samples, calculate_CtoO_and_metallicity(samples[:, gas].T).T, axis=-1)
    # parameter_names.append('C/O')
    # parameter_names.append('[Fe/H]')
    print(list(zip(parameter_names, samples[np.argmax(likelihoods)])))

    default_precisions = [2] * np.shape(samples)[1]
    default_precisions[temps] = [0] * (temp_index_bounds[1] - temp_index_bounds[0])
    # default_precisions[-1] = 0

    print()
    samples_low, samples_median, samples_high = np.nanpercentile(
        samples, [16, 50, 84], axis=0
    )
    lower_error = np.abs(samples_median - samples_low)
    upper_error = np.abs(samples_high - samples_median)
    max_precisions = np.max(
        (
            np.ceil(-np.log10(lower_error)).astype(int),
            np.ceil(-np.log10(upper_error)).astype(int),
        ),
        axis=0,
    )
    string_precisions = np.where(
        max_precisions > default_precisions, max_precisions, default_precisions
    )
    string_formats = [
        ".{}f".format(string_precisions) for string_precisions in string_precisions
    ]
    string_format_dicts[directory] = string_formats

    # samples_low, samples_median, samples_high = np.nanpercentile(samples, [2.5,50,97.5], axis=0)
    samples_low, samples_median, samples_high = np.nanpercentile(
        samples, [16, 50, 84], axis=0
    )
    lower_error = np.abs(samples_median - samples_low)
    upper_error = np.abs(samples_high - samples_median)
    max_precisions = np.max(
        (
            np.ceil(-np.log10(lower_error)).astype(int),
            np.ceil(-np.log10(upper_error)).astype(int),
        ),
        axis=0,
    )

    sample_dicts[directory], npars = parse_samples(
        samples, parameter_names, likelihoods
    )
    for group_name, sample_group in sample_dicts[directory].items():
        if sample_group[0]:
            print(group_name.capitalize())
            print(
                [
                    name
                    + ": {1:.5f} +{2:.5f} - {3:.5f}, MLE {0:.5f}".format(
                        value, median, high, low
                    )
                    for name, value, median, high, low in zip(
                        sample_group[1],
                        sample_group[3],
                        sample_group[5],
                        sample_group[7],
                        sample_group[6],
                    )
                ]
            )

group_string_formats = {
    run_name: {
        group_name: [
            string_format
            for string_format, parameter_name in zip(
                string_format_dicts[run_name], parameter_names
            )
            if parameter_name in group[1]
        ]
        for group_name, group in run.items()
    }
    for (run_name, run), parameter_names in zip(
        sample_dicts.items(), [cloudy_parameter_names]
    )
}

cmap_kwargs = {
    "lightness_minimum": 0.15,
    "lightness_maximum": 0.85,
    "saturation_minimum": 0.2,
    "saturation_maximum": 0.8,
}

cmap_H2O = create_linear_colormap(["#226666", "#2E4172"], **cmap_kwargs)
cmap_CO = create_linear_colormap(["#882D60", "#AA3939"], **cmap_kwargs)
cmap_CO2 = create_linear_colormap(["#96A537", "#669933"], **cmap_kwargs)
cmap_CH4 = create_linear_colormap(["#96A537", "#669933"], **cmap_kwargs)

cmap_cloudy = create_linear_colormap(
    [cnames["lightcoral"], cnames["lightcoral"]], **cmap_kwargs
)
cmap_clear = create_linear_colormap(
    [cnames["cornflowerblue"], cnames["cornflowerblue"]], **cmap_kwargs
)

cmap_cloud = plt.get_cmap("Greys")

plotted_components = ["h2o", "co", "co2", "ch4"]
plotted_titles = ["H$_2$O", "CO", "CO$_2$", "CH$_4$"]
cmaps = [cmap_H2O, cmap_CO, cmap_CO2, cmap_CH4]

JHK_contributions = {}
for directory, object_header in zip(directories, object_headers):
    # JHK_contribution_file = directory+"/"+object_header+"_contributions.p"
    JHK_contribution_file = DIRECTORY_RESULTS_2M2236 / CONTRIBUTIONS_FILE_2M2236
    # JHK_contribution_file = "mock-L/input/mock-L.2022-12-08.forward-model.PLO.JHK_contributions.p"
    with open(JHK_contribution_file, "rb") as pickle_file:
        JHK_contributions[directory] = pickle.load(pickle_file)

contributions_max = np.log10(
    np.nanmax(
        [
            np.nanmax(contribution)
            for (species, contribution) in JHK_contributions[directory].items()
            if species not in ["gas", "cloud", "total"]
        ]
    )
)
contributions_min = contributions_max - 3

wavelengths_low, wavelengths_high, data, data_error = (
    np.genfromtxt(DIRECTORY_RESULTS_2M2236 / OBSERVED_SPECTRUM_FILE_2M2236).T
)[:4]

band_breaks = np.r_[
    0,
    np.nonzero(
        (
            JHK_contributions[directory]["h2o"].index.to_numpy()[1:]
            - JHK_contributions[directory]["h2o"].index.to_numpy()[:-1]
        )
        > 0.05
    )[0]
    + 1,
    len(JHK_contributions[directory]["h2o"].index),
]
break_indices = np.where(wavelengths_low[1:] != wavelengths_high[:-1])[0]
mask = break_indices

padding = 0.025
number_of_bands = 2
wavelength_ranges = np.array(
    [
        JHK_contributions[directory]["h2o"].index[band_breaks[i + 1] - 1]
        - JHK_contributions[directory]["h2o"].index[band_breaks[i]]
        for i in range(number_of_bands)
    ]
)

fig = plt.figure(figsize=(40, 30))
gs = fig.add_gridspec(
    nrows=6,
    ncols=2,
    height_ratios=[4, 2, 3, 3, 3, 3],
    width_ratios=(1 + 2 * padding) * wavelength_ranges,
    wspace=0.1,
)

spectrum_lines = []
spectrum_axes = []
residual_axes = []
contribution_columns = []
for j, (directory, object_header, data_file, color, model_title) in enumerate(
    zip(
        directories,
        object_headers,
        data_files,
        plotting_colors,
        ["Cloudy model", "Clear model"],
    )
):
    spectrum_dict = extract_spectra(
        directory,
        object_header,
        DIRECTORY_RESULTS_2M2236 / MLE_SPECTRUM_FILE_2M2236,
        DIRECTORY_RESULTS_2M2236 / OBSERVED_SPECTRUM_FILE_2M2236,
    )

    datas = spectrum_dict["data"]
    models = spectrum_dict["MLE spectrum"]
    errors = spectrum_dict["data errors"]
    residuals = (models - datas) / errors

    # ymin = np.min([np.min(data), np.min(model)])
    # ymax = np.max([np.max(data), np.max(model)])
    # ymin = ymin - padding*np.abs(ymax-ymin)
    # ymax = ymax + padding*np.abs(ymax-ymin)

    residual_ymin = np.nanmin(residuals)
    residual_ymax = np.nanmax(residuals)
    residual_ymin = residual_ymin - padding * np.abs(residual_ymax - residual_ymin)
    residual_ymax = residual_ymax + padding * np.abs(residual_ymax - residual_ymin)

    concocted_boundaries = [[2.85, 4.01], [4.19, 5.30]]
    for i, band_boundaries in enumerate(concocted_boundaries):
        if j == 0:
            if i > 0:
                spectrum_ax = fig.add_subplot(gs[0, i])  # sharey=spectrum_axes[i-1])
                spectrum_axes.append(spectrum_ax)

                residual_ax = fig.add_subplot(
                    gs[1, i], sharex=spectrum_ax
                )  # sharey=residual_axes[i-1])
                residual_axes.append(residual_ax)
            else:
                spectrum_ax = fig.add_subplot(gs[0, i])
                spectrum_axes.append(spectrum_ax)

                residual_ax = fig.add_subplot(gs[1, i], sharex=spectrum_ax)
                residual_axes.append(residual_ax)
        else:
            spectrum_ax = spectrum_axes[i]
            residual_ax = residual_axes[i]

        # if i==len(concocted_boundaries)-1:
        #    spectrum_ax.text(
        #        0.975, 0.95-0.15*j,
        #        model_title,
        #        horizontalalignment="right",
        #        verticalalignment="top",
        #        transform=spectrum_ax.transAxes,
        #        fontsize=36,
        #        color=color
        #        )

        band_condition = (spectrum_dict["wavelengths"] > band_boundaries[0]) & (
            spectrum_dict["wavelengths"] < band_boundaries[1]
        )
        wavelengths = spectrum_dict["wavelengths"][band_condition]
        data = spectrum_dict["data"][band_condition]
        error = spectrum_dict["data errors"][band_condition]
        model = spectrum_dict["MLE spectrum"][band_condition]
        residual = residuals[band_condition]

        xmin = np.min(wavelengths)
        xmax = np.max(wavelengths)
        xmin = xmin - padding * np.abs(xmax - xmin)
        xmax = xmax + padding * np.abs(xmax - xmin)

        spectrum_ax.set_xlim([xmin, xmax])
        # spectrum_ax.set_ylim([ymin, ymax])

        linestyles = ["solid", "solid"]
        spectrum_ax.errorbar(
            wavelengths,
            data,
            error,
            color="#444444",
            fmt="x",
            linewidth=0,
            elinewidth=2,
            alpha=1,
            zorder=-3,
        )
        spectrum_ax.plot(
            wavelengths,
            model,
            color=color,
            linewidth=3,
            linestyle=linestyles[j],
            alpha=1,
            zorder=2 - j,
            label=model_title,
        )

        # if i==0 and j==0:
        # [spectrum_ax.axvline(line_position, linestyle="dashed", linewidth=1.5, zorder=-10, color="#888888")
        # for line_position in [1.139, 1.141, 1.169, 1.177, 1.244, 1.253, 1.268]]
        # y_text = spectrum_ax.get_ylim()[0] + 0.1*np.diff(spectrum_ax.get_ylim())
        # spectrum_ax.text((1.169+1.177)/2, y_text, "KI", fontsize=20, horizontalalignment="center")
        # spectrum_ax.text((1.244+1.253)/2, y_text, "KI", fontsize=20, horizontalalignment="center")
        # spectrum_ax.text(1.268, y_text, "NaI", fontsize=20, horizontalalignment="center")

        residual_ax.plot(
            wavelengths,
            residual,
            color=color,
            linewidth=3,
            linestyle=linestyles[j],
            alpha=1,
            zorder=2 - j,
        )
        residual_ax.axhline(
            0, color="#444444", linewidth=2, linestyle="dashed", zorder=-10
        )
        residual_ax.set_ylim([residual_ymin, residual_ymax])
        residual_ax.minorticks_on()

        if i == 0:
            spectrum_ax.set_ylabel("Flux (erg s$^{-1}$ cm$^{-3}$)", fontsize=36)
            residual_ax.set_ylabel(r"Residual/$\sigma$", fontsize=26)
        if (j == 1) and (i == len(concocted_boundaries) - 1):
            legend_handles, legend_labels = spectrum_ax.get_legend_handles_labels()
            legend_handles.append(Patch(facecolor="k", edgecolor="#DDDDDD"))
            legend_labels.append(r"Cloud $\tau > 0.1$")
            spectrum_ax.legend(
                handles=legend_handles, labels=legend_labels, fontsize=22, frameon=False
            )  # , loc="upper center")

        # residual = (model-data)/error
        reduced_chi_squared = np.sum(residual**2) / (np.shape(residual)[0] - npars)
        print("Reduced chi squared is {}".format(reduced_chi_squared))

        spectrum_ax.tick_params(axis="x", labelsize=26)
        spectrum_ax.tick_params(axis="y", labelsize=26)
        spectrum_ax.yaxis.get_offset_text().set_size(26)

        residual_ax.tick_params(axis="x", labelsize=26)
        residual_ax.tick_params(axis="y", labelsize=26)

        # CONTRIBUTION PLOTS

        contribution_axes = []
        for n, (comp, cmap, title) in enumerate(
            zip(plotted_components, cmaps, plotted_titles)
        ):
            if j == 0:
                contribution_ax = fig.add_subplot(gs[n + 2, i], sharex=spectrum_axes[i])
                contribution_axes.append(contribution_ax)
                if n == len(plotted_components) - 1:
                    contribution_columns.append(contribution_axes)
                if i == 0:
                    contribution_ax.text(
                        0.025,
                        0.05,
                        title + " contribution",
                        horizontalalignment="left",
                        verticalalignment="bottom",
                        transform=contribution_ax.transAxes,
                        fontsize=32,
                    )
            else:
                contribution_ax = contribution_columns[i][n]

            if j == 0:
                cont_cmap = cmap_cloudy
                outline_color = "lightcoral"
            else:
                cont_cmap = cmap_clear
                outline_color = "cornflowerblue"

            cf = JHK_contributions[directory][comp].to_numpy()
            x, y = np.meshgrid(
                JHK_contributions[directory][comp].index,
                JHK_contributions[directory][comp].columns,
            )
            contribution_ax.contourf(
                x[:, band_breaks[i] : band_breaks[i + 1] : 8],
                y[:, band_breaks[i] : band_breaks[i + 1] : 8],
                np.log10(cf).T[:, band_breaks[i] : band_breaks[i + 1] : 8],
                cmap=cont_cmap,
                levels=contributions_max - np.array([4, 2, 0]),
                alpha=0.66,
                zorder=0,
            )
            contribution_ax.contour(
                x[:, band_breaks[i] : band_breaks[i + 1] : 8],
                y[:, band_breaks[i] : band_breaks[i + 1] : 8],
                np.log10(cf).T[:, band_breaks[i] : band_breaks[i + 1] : 8],
                colors=outline_color,
                levels=contributions_max - np.array([2]),
                linewidths=3,
                alpha=1,
                zorder=1,
            )
            if j == 0:
                contribution_ax.invert_yaxis()
            if "cloud" in JHK_contributions[directory]:
                cloud_cf = JHK_contributions[directory]["cloud"]
                x, y = np.meshgrid(cloud_cf.index, cloud_cf.columns)
                contribution_ax.contourf(
                    x[:, band_breaks[i] : band_breaks[i + 1] : 8],
                    y[:, band_breaks[i] : band_breaks[i + 1] : 8],
                    cloud_cf.to_numpy().T[:, band_breaks[i] : band_breaks[i + 1] : 8],
                    colors="k",
                    # cmap=cmap_cloud,
                    alpha=0.75,
                    # levels=np.logspace(-1, 2, num=20),
                    levels=[0.1, 0.75],
                    zorder=2,
                )
                if i == 0:
                    contribution_ax.contour(
                        x[:, band_breaks[i] : band_breaks[i + 1] : 8],
                        y[:, band_breaks[i] : band_breaks[i + 1] : 8],
                        cloud_cf.to_numpy().T[
                            :, band_breaks[i] : band_breaks[i + 1] : 8
                        ],
                        colors="#DDDDDD",
                        linestyles="solid",
                        linewidth=3,
                        alpha=1,
                        levels=[0.1],
                        zorder=3,
                    )
                contribution_ax.tick_params(axis="x", labelsize=26)
                contribution_ax.tick_params(axis="y", labelsize=26)
                contribution_ax.minorticks_on()
            # contribution_ax.contour(x[:, band_breaks[i]:band_breaks[i+1]:8], y[:, band_breaks[i]:band_breaks[i+1]:8],
            #                       np.log10(cf).T[:, band_breaks[i]:band_breaks[i+1]:8],
            #                       cmap=cmap,
            #                       levels=contributions_max-np.array([3, 2, 1, 0]),
            #                       alpha=1,
            #                       zorder=0)

contributions_ax = fig.add_subplot(gs[2:, :])
contributions_ax.spines["top"].set_color("none")
contributions_ax.spines["bottom"].set_color("none")
contributions_ax.spines["left"].set_color("none")
contributions_ax.spines["right"].set_color("none")
contributions_ax.tick_params(
    labelcolor="none", top=False, bottom=False, left=False, right=False
)
contributions_ax.grid(False)
contributions_ax.set_xlabel(
    r"$\lambda\left(\mu\mathrm{m}\right)$", fontsize=36, y=-0.075
)
contributions_ax.set_ylabel(
    r"$\log_{10}\left(P/\mathrm{bar}\right)$", fontsize=36, labelpad=20
)

"""
overall_ax = fig.add_subplot(111)
overall_ax.spines['top'].set_color('none')
overall_ax.spines['bottom'].set_color('none')
overall_ax.spines['left'].set_color('none')
overall_ax.spines['right'].set_color('none')
overall_ax.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
overall_ax.grid(False)
overall_ax.set_xlabel(r"$\lambda\left(\mu\mathrm{m}\right)$", fontsize=36, y=-0.075)
"""


for filetype in filetypes:
    plt.savefig(
        object_label + ".fit-spectrum+contributions.{}".format(filetype),
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )

fig, ax = plt.subplots(figsize=(8, 6))
for (band, samples), color, label in zip(
    sample_dicts.items(), plotting_colors, ["cloudy", "clear"]
):
    T_samples = (samples["T-P"][2]).T
    MLE_Tsample = samples["T-P"][3]
    logP, [T_minus_1sigma, T_median, T_plus_1sigma] = generate_profiles(
        T_samples, Piette
    )
    # ax.plot(T_median, logP, color=color, linewidth=1.5)
    ax.fill_betweenx(
        logP,
        T_minus_1sigma,
        T_plus_1sigma,
        linewidth=0.5,
        color=color,
        facecolor=color,
        alpha=0.5,
        label="95\% confidence interval, {}".format(label),
    )
    logP, MLE_Tprofile = generate_profiles(MLE_Tsample, Piette)
    # ax.plot(T_median, logP, color=color, linewidth=1.5, label="Median profile")
    ax.plot(MLE_Tprofile, logP, color=color, linewidth=4)
    ax.plot(MLE_Tprofile, logP, color=MLE_color, linewidth=2, label="MLE profiles")
    ax.set_ylim([-4, 2.5])
    # reference_TP = true_values[temps]
    # logP, reference_Tprofile = generate_profiles(reference_TP, Piette)
    # ax.plot(reference_Tprofile, logP, color="#444444", linewidth=3)
    # ax.plot(reference_Tprofile, logP, color=reference_color, linewidth=2, label="True profile")
# ax.plot(sonora_T, sonora_P, linewidth=2, linestyle="dashed", color="sandybrown", label=r"SONORA, $\log g$=3.67, $T_\mathrm{eff}$=1584 K", zorder=-8)
# ax.plot(sonora_T, sonora_P, linewidth=4, color="saddlebrown", label="", zorder=-9)
cloud_MLEs = dict(
    list(zip(sample_dicts["2M2236"]["clouds"][1], sample_dicts["2M2236"]["clouds"][3]))
)

fig.gca().invert_yaxis()
ax.set_xlim(ax.get_xlim()[0], 4000)
if cloud_MLEs:
    ax.fill_between(
        np.linspace(ax.get_xlim()[0], 4000),
        cloud_MLEs["Haze_minP"],
        cloud_MLEs["Haze_minP"] + cloud_MLEs["Haze_thick"],
        color="#444444",
        alpha=0.5,
        label="Retrieved cloud layer",
        zorder=-10,
    )

ax.set_xlabel("Temperature (K)")
ax.set_ylabel("log$_{10}$(Pressure/bar)")
handles, labels = plt.gca().get_legend_handles_labels()
order = [1, 0]  # , 4, 5]
# order=[1, 2, 0]
ax.legend(
    [handles[idx] for idx in order],
    [labels[idx] for idx in order],
    fontsize=11,
    facecolor="#444444",
    framealpha=0.25,
    loc="center right",
)

for filetype in filetypes:
    plt.savefig(
        object_label + ".T-P_profiles.{}".format(filetype), bbox_inches="tight", dpi=300
    )

for group_name in sample_dicts[directory].keys():
    print("Generating corner plots for " + group_name.capitalize())
    for i, (band, sample_dict) in enumerate(sample_dicts.items()):
        if not sample_dict[group_name][0]:
            print("Empty group. Skipping corner plots.")
            continue
        print_names, names, samples, MLE_sample, range, _, lower_error, upper_error = (
            sample_dict[group_name]
        )
        print(samples.shape)
        reference_values = MLE_sample
        # reference_values = reference_dict[group_name][3]
        color = plotting_colors[directories.index(band)]
        plot_generic_legend_first = group_name == "clouds"
        # if i == 0:
        if False:
            base_weights = importance_dicts[band]
            fig, par_titles, title_color = generate_cornerplot(
                samples,
                weights=base_weights,
                group_name=group_name,
                parameter_names=print_names,
                parameter_range=range,
                confidence=0.95,
                color=color,
                MLE_values=MLE_sample,
                MLE_name="Cloudy model",
                MLE_color=color,
                overtext_y=0.8,
                string_formats=group_string_formats[band][group_name],
                # reference_values=reference_values,
                # reference_name=reference_name,
                # reference_markerstyle=reference_markerstyle,
                # reference_color=reference_color,
                reference_values=None,
                reference_name=None,
                reference_markerstyle=None,
                reference_color=None,
                plot_generic_legend_labels=plot_generic_legend_first,
            )
        else:
            print(group_string_formats[band][group_name])
            weights = importance_dicts[band]
            fig, _, _ = generate_cornerplot(
                samples,
                weights=weights,
                group_name=group_name,
                parameter_names=print_names,
                parameter_range=range,
                confidence=0.95,
                # existing_figure=fig,
                # existing_titles=par_titles,
                # existing_title_color=title_color,
                color=color,
                MLE_values=MLE_sample,
                MLE_name="Clear model",
                MLE_color=color,
                string_formats=group_string_formats[band][group_name],
                # reference_values=reference_values,
                # reference_name=reference_name,
                # reference_markerstyle=reference_markerstyle,
                # reference_color=reference_color,
                reference_values=None,
                reference_name=None,
                reference_markerstyle=None,
                reference_color=None,
                plot_generic_legend_labels=True,
            )

    for filetype in filetypes:
        fig.savefig(
            object_label + "." + group_name + ".{}".format(filetype),
            bbox_inches="tight",
            dpi=300,
        )
