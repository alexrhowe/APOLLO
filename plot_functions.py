import corner
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import cnames
import matplotlib.patheffects as pe
import numpy as np
import sys


def extract_spectra(
    directory,
    object_header,
    MLE_spec_filepath,
    data_filepath,
    steps_per_file=100,
    max_memory_size=50000,
):
    _, _, MLE_spectrum = np.genfromtxt(MLE_spec_filepath)[:, :3].T

    # spec_file = glob(directory+"/"+object_header+"*ensemble.dat")[0]
    # MCMC_all = np.genfromtxt(spec_file)
    # MCMC_waves = (MCMC_all[:, 0]+MCMC_all[:, 1])/2.
    # MCMC_spectra = MCMC_all[:, 2:].T

    wavens_hi, wavens_lo, data, error = np.genfromtxt(data_filepath)[:, :4].T
    waves = (wavens_hi + wavens_lo) / 2

    return {
        "MLE spectrum": MLE_spectrum,
        "wavelengths": waves,
        "data": data,
        "data errors": error,
    }
    # "MCMC waves": MCMC_waves,
    # "MCMC spectra": MCMC_spectra}


def generate_cornerplot(
    samples,
    weights,
    group_name,
    parameter_names,
    parameter_range,
    confidence=0.99,
    existing_figure=None,
    existing_titles=None,
    existing_title_color=None,
    color="#53C79B",
    MLE_color="gold",
    MLE_values=None,
    MLE_name="MLE value",
    overtext_y=0.5,
    string_formats=".4g",
    reference_values=None,
    reference_lowerbound=None,
    reference_upperbound=None,
    reference_name="Reference",
    reference_markerstyle="*",
    reference_color="gold",
    plot_generic_legend_labels=False,
):
    if existing_figure is None:
        fig = corner.corner(
            samples,
            weights=weights,
            bins=20,
            color=color,
            labels=parameter_names,
            quantiles=[0.16, 0.84],
            range=parameter_range,
            # smooth=0.3,
            # smooth1d=0.3,
            hist_kwargs={"facecolor": color, "alpha": 0.5},
            show_titles=True,
            title_fmt=string_formats[0],
            title_kwargs={"pad": 0, "color": color},
            title_quantiles=[0.16, 0.84],
            plot_datapoints=True,
            labelsize=11,
            use_math_text=True,
        )
        fig.subplots_adjust(left=0.10, bottom=0.10, wspace=0.075, hspace=0.075)
    else:
        fig = corner.corner(
            samples,
            fig=existing_figure,
            weights=weights,
            bins=20,
            color=color,
            labels=parameter_names,
            quantiles=[0.16, 0.84],
            range=parameter_range,
            # smooth=0.3,
            # smooth1d=0.3,
            hist_kwargs={"facecolor": color, "alpha": 0.5},
            show_titles=True,
            title_fmt=string_formats[0],
            title_kwargs={"pad": 10, "color": color},
            title_quantiles=[0.16, 0.84],
            plot_datapoints=True,
            labelsize=11,
            use_math_text=True,
        )
        fig.subplots_adjust(left=0.10, bottom=0.10, wspace=0.075, hspace=0.075)

    par_titles = []

    if reference_values is None:
        reference_values = np.full_like(MLE_values, np.nan)

    # Extract the axes
    ndim = len(parameter_names)
    axes = np.array(fig.axes).reshape((ndim, ndim))

    # if MLE_values is not None:
    #    retrieved_values = MLE_values
    # else:
    #    retrieved_values = np.nanpercentile(samples,50,axis=0)
    retrieved_values = np.nanpercentile(samples, 50, axis=0)

    cumulative_xlims = []
    # Loop over the diagonal
    for i in range(ndim):
        ax = axes[i, i]
        xlim = ax.xaxis.get_data_interval()
        if not np.isnan(reference_values[i]):
            xlim[0] = np.min(
                [xlim[0], reference_values[i] - 0.1 * np.abs(np.ptp(xlim))]
            )
            xlim[1] = np.max(
                [xlim[1], reference_values[i] + 0.1 * np.abs(np.ptp(xlim))]
            )
            formatted_syntax = r"{:" + string_formats[i] + r"}"
            reference_overtext = ax.text(
                reference_values[i],
                np.sum(ax.yaxis.get_data_interval()) / 5,
                formatted_syntax.format(reference_values[i]),
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=18,
                color=reference_color,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="#444444")],
            )
            if reference_lowerbound is not None:
                ax.axvline(
                    reference_lowerbound[i],
                    linestyle="dashed",
                    linewidth=2,
                    color=reference_color,
                )
                if reference_lowerbound[i] < xlim[0]:
                    xlim[0] = reference_lowerbound[i]
            if reference_upperbound is not None:
                ax.axvline(
                    reference_upperbound[i],
                    linestyle="dashed",
                    linewidth=2,
                    color=reference_color,
                )
                if reference_upperbound[i] > xlim[1]:
                    xlim[1] = reference_upperbound[i]

        formatted_syntax = r"{:" + string_formats[i] + r"}"
        reference_overtext = ax.text(
            MLE_values[i],
            overtext_y * np.sum(ax.yaxis.get_data_interval()),
            formatted_syntax.format(MLE_values[i]),
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=18,
            color=MLE_color,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="#444444")],
        )

        ax.set_xlim(xlim)
        cumulative_xlims.append(xlim)
        ylim = ax.yaxis.get_data_interval()
        ax.set_ylim(ylim)
        ax.margins(0.05)

        if group_name == "clouds":
            hist_title_fontsize = 12
        else:
            hist_title_fontsize = 16

        # USE CONDITIONAL FOR OVERPLOTTING
        if existing_figure is None:
            par_title = ax.get_title()
            ax.set_title(
                par_title,
                fontsize=hist_title_fontsize,
                color=color,
                path_effects=[pe.withStroke(linewidth=1, foreground="#444444")],
            )
            par_titles.append(par_title)

        else:
            par_overtitle = ax.get_title()
            ax.text(
                0.5,
                1.25,
                par_overtitle,
                fontsize=hist_title_fontsize,
                color=color,
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                path_effects=[pe.withStroke(linewidth=1, foreground="#444444")],
            )
            if existing_titles is not None:
                ax.set_title(
                    existing_titles[i],
                    fontsize=hist_title_fontsize,
                    color=existing_title_color,
                    path_effects=[pe.withStroke(linewidth=1, foreground="#444444")],
                )

        title_color = color

        # ax.axvline(retrieved_values[i], linewidth=3, color="#444444")
        # ax.axvline(retrieved_values[i], color=color)

        ax.axvline(MLE_values[i], linewidth=3, color="#444444")
        ax.axvline(MLE_values[i], color=MLE_color)

        if not np.isnan(reference_values[i]):
            ax.axvline(reference_values[i], linewidth=3, color="#444444")
            ax.axvline(reference_values[i], color=reference_color)

    # Loop over the histograms
    for yi in range(ndim):
        for xi, xlim in enumerate(cumulative_xlims[:yi]):
            ax = axes[yi, xi]
            ax.set_xlim(xlim)
            ylim = ax.yaxis.get_data_interval()
            if not np.isnan(reference_values[yi]):
                ylim[0] = np.min(
                    [ylim[0], reference_values[yi] - 0.1 * np.abs(np.ptp(ylim))]
                )
                ylim[1] = np.max(
                    [ylim[1], reference_values[yi] + 0.1 * np.abs(np.ptp(ylim))]
                )
            ax.set_ylim(ylim)
            ax.margins(0.05)

            # ax.axvline(retrieved_values[xi], linewidth=3, color="#444444")
            # ax.axhline(retrieved_values[yi], linewidth=3, color="#444444")
            # ax.axvline(retrieved_values[xi], color=color)
            # ax.axhline(retrieved_values[yi], color=color)

            ax.axvline(MLE_values[xi], linewidth=3, color="#444444")
            ax.axhline(MLE_values[yi], linewidth=3, color="#444444")
            ax.axvline(MLE_values[xi], color=MLE_color)
            # MLE_locator = ax.axhline(MLE_values[yi], color=MLE_color, label=MLE_name)
            MLE_locator = ax.axhline(MLE_values[yi], color=MLE_color)

            if not np.isnan(reference_values[yi]) and not np.isnan(
                reference_values[xi]
            ):
                ax.axvline(reference_values[xi], linewidth=3, color="#444444")
                ax.axvline(reference_values[xi], color=reference_color)
                ax.axhline(reference_values[yi], linewidth=3, color="#444444")
                ax.axhline(reference_values[yi], color=reference_color)

            # ax.plot(retrieved_values[xi], retrieved_values[yi], marker="o", markersize=10, color="#444444")
            # median_marker, = ax.plot(retrieved_values[xi], retrieved_values[yi], marker="o", markersize=8, color=color)
            ax.plot(
                MLE_values[xi],
                MLE_values[yi],
                marker="D",
                markersize=10,
                color="#444444",
            )
            (MLE_marker,) = ax.plot(
                MLE_values[xi],
                MLE_values[yi],
                marker="D",
                markersize=8,
                color=MLE_color,
            )

            if not np.isnan(reference_values[yi]) and not np.isnan(
                reference_values[xi]
            ):
                if xi == 0 and yi == 1 and reference_name is not None:
                    (reference_locator,) = ax.plot(
                        reference_values[xi],
                        reference_values[yi],
                        marker=reference_markerstyle,
                        markersize=10,
                        color="#444444",
                    )
                    (reference_marker,) = ax.plot(
                        reference_values[xi],
                        reference_values[yi],
                        marker=reference_markerstyle,
                        markersize=8,
                        color=reference_color,
                    )
                    print("We're setting the reference label!")
                else:
                    ax.plot(
                        reference_values[xi],
                        reference_values[yi],
                        marker=reference_markerstyle,
                        markersize=10,
                        color="#444444",
                    )
                    ax.plot(
                        reference_values[xi],
                        reference_values[yi],
                        marker=reference_markerstyle,
                        markersize=8,
                        color=reference_color,
                    )

    if reference_name is not None:
        reference_marker.set_label("{} value".format(reference_name))
    # MLE_locator.set_label(MLE_name)

    if plot_generic_legend_labels:
        legend_handles, legend_labels = axes[1, 0].get_legend_handles_labels()
        legend_handles.append(
            Line2D([0], [0], marker="D", markersize=10, color="#444444")
        )
        # legend_handles.append(Line2D([0], [0], marker="D", markersize=10, color=MLE_color))
        legend_labels.append("MLE value")
        legend_handles.append(
            Line2D([0], [0], linestyle="dashed", linewidth=3, color="#444444")
        )
        # legend_handles.append(Line2D([0], [0], linestyle="dashed", linewidth=3, color=MLE_color))
        legend_labels.append("68\% confidence interval")

    # median_marker.set_label("Median value")
    # USE CONDITIONAL FOR OVERPLOTTING
    if plot_generic_legend_labels:
        plt.figlegend(
            handles=legend_handles,
            labels=legend_labels,
            fontsize=28,
            loc="upper right",
            facecolor="#444444",
            framealpha=0.25,
        )
    return fig, par_titles, title_color


def generate_profiles(samples, profile_function, minP=-4, maxP=2.5, num_points=1000):
    logP = np.linspace(minP, maxP, num=num_points)
    if len(np.shape(samples)) > 1:
        # logP_array = logP.reshape(np.r_[logP.shape, [1]*(samples.ndim-1)])
        T_profiles = profile_function(
            *samples, num_layers_final=num_points, P_min=minP, P_max=maxP
        )
        return logP, np.nanpercentile(
            T_profiles, [2.5, 50, 97.5], axis=np.arange(T_profiles.ndim - 1) + 1
        )
    else:
        T_profile = profile_function(
            *samples, num_layers_final=num_points, P_min=minP, P_max=maxP
        )
        return logP, T_profile


def parse_samples(samples, pnames, likelihoods, confidence=0.99):
    likelihood_index = np.argmax(likelihoods)
    print("Maximum posterior probability is {}".format(likelihoods[likelihood_index]))
    MLE_sample = samples[likelihood_index]

    baselist = ["Rad", "RtoD", "RtoD2U", "Log(g)", "Mass", "C/O", "[Fe/H]", "Teff"]
    bnamelist = [
        r"$R/R_\mathrm{J}$",
        r"R/R$_\mathrm{J}$",
        r"R/R$_\mathrm{J}$",
        r"$\log g$",
        r"$M/M_\mathrm{J}$",
        "C/O",
        "[Fe/H]",
        r"T$_\mathrm{eff}$ (K)",
    ]
    cloudlist = [
        "Haze_abund",
        "Haze_expo",
        "Haze_size",
        "Haze_minP",
        "Haze_thick",
        "Haze_tau",
        "Haze_w0",
        "Cloud_Fraction",
    ]
    cnamelist = [
        r"log$_{10}$($n_\mathrm{cloud}$/cm$^{-3}$)",
        r"$\alpha$",
        r"log$_{10}$($r_\mathrm{part}$/$\mu$m)",
        r"log$_{10}$(P$_\mathrm{top}$/bar)",
        r"log$_{10}$($\Delta$P$_\mathrm{cloud}$/bar)",
        r"$\log_{10}\!\left[\tau\!\left(\lambda_0\right)\right]$",
        r"$\omega_0$",
        r"$f_\mathrm{cloud}$",
    ]
    gaslist = [
        "h2",
        "h2only",
        "he",
        "h-",
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
    gnamelist = [
        r"H$_2$+He",
        r"H$_2$",
        "He",
        "[H-]",
        r"[H$_2$O]",
        r"[CH$_4$]",
        "[CO]",
        r"[CO$_2$]",
        r"[NH$_3$]",
        r"[H$_2$S]",
        "[Na,K]",
        "[Na,K]",
        "[CrH]",
        "[FeH]",
        "[TiO]",
        "[VO]",
        "[HCN]",
        "[N2]",
        "[PH3]",
    ]
    templist = [
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
    ]
    tnamelist = [
        r"$T_{-4}$",
        r"$T_{-3}$",
        r"$T_{-2}$",
        r"$T_{-1}$",
        r"$T_{0}$",
        r"$T_{0.5}$",
        r"$T_{1}$",
        r"$T_{1.5}$",
        r"$T_{2}$",
        r"$T_{2.5}$",
    ]

    endlist = ["deltaL", "logf", "scaleG395_ch1", "scaleG395_ch2"]
    enamelist = [r"$\Delta \lambda$ (nm)", r"$\log f$", r"$\Delta_1$", r"$\Delta_2$"]

    # customlist = ['Rad', 'Log(g)', 'h2o', 'co', 'co2', 'C/O', 'Teff']
    # knamelist = ['$R/R_\mathrm{J}$', '$\log g$', '[H$_2$O]', '[CO]', '[CO$_2$]', 'C/O', 'T$_\mathrm{eff}$ (K)']

    # customlist = baselist + gaslist + templist + endlist
    # knamelist = bnamelist + gnamelist + tnamelist + enamelist

    customlist = endlist
    knamelist = enamelist

    basic = []
    bnames = []
    bnames2 = []
    cloud = []
    cnames = []
    cnames2 = []
    gases = []
    gnames = []
    gnames2 = []
    temps = []
    tnames = []
    tnames2 = []
    # A separate corner plot for custom comparisons across existing categories (e.g. H2O abundance vs. cloud deck pressure).
    custom = []
    knames = []
    knames2 = []

    for i in range(0, len(pnames)):
        # pnames[i] = pnames[i].decode('UTF-8')
        if pnames[i] in baselist:
            basic.append(i)
            j = baselist.index(pnames[i])
            bnames.append(bnamelist[j])
            bnames2.append(pnames[i])

        if pnames[i] in cloudlist:
            cloud.append(i)
            j = cloudlist.index(pnames[i])
            cnames.append(cnamelist[j])
            cnames2.append(pnames[i])

        if pnames[i] in gaslist:
            gases.append(i)
            j = gaslist.index(pnames[i])
            gnames.append(gnamelist[j])
            gnames2.append(pnames[i])

        if pnames[i] in templist:
            temps.append(i)
            j = templist.index(pnames[i])
            tnames.append(tnamelist[j])
            tnames2.append(pnames[i])

        if pnames[i] in customlist:
            custom.append(i)
            j = customlist.index(pnames[i])
            knames.append(knamelist[j])
            knames2.append(pnames[i])

        # Used for layer-by-layer T-P profiles.
        # if pnames[i][0]=='T' and not pnames[i]=='Teff': temps.append(i)

    bsamples = samples[..., basic]
    csamples = samples[..., cloud]
    gsamples = samples[..., gases]
    tsamples = samples[..., temps]
    ksamples = samples[..., custom]

    MLE_bsample = MLE_sample[basic]
    MLE_csample = MLE_sample[cloud]
    MLE_gsample = MLE_sample[gases]
    MLE_tsample = MLE_sample[temps]
    MLE_ksample = MLE_sample[custom]
    brange = np.zeros(len(bnames))
    for i in range(0, len(bnames)):
        brange[i] = confidence
    crange = np.zeros(len(cnames))
    for i in range(0, len(cnames)):
        crange[i] = confidence
    grange = np.zeros(len(gnames))
    for i in range(0, len(gnames)):
        grange[i] = confidence
    trange = np.zeros(len(tnames))
    for i in range(0, len(tnames)):
        trange[i] = confidence
    krange = np.zeros(len(knames))
    for i in range(0, len(knames)):
        krange[i] = confidence

    blow, bmedian, bhigh = np.nanpercentile(bsamples, [16, 50, 84], axis=0)
    if np.shape(csamples)[-1] > 0:
        clow, cmedian, chigh = np.nanpercentile(csamples, [16, 50, 84], axis=0)
    else:
        clow = 0
        cmedian = 0
        chigh = 0
    glow, gmedian, ghigh = np.nanpercentile(gsamples, [16, 50, 84], axis=0)
    tlow, tmedian, thigh = np.nanpercentile(tsamples, [16, 50, 84], axis=0)
    klow, kmedian, khigh = np.nanpercentile(ksamples, [16, 50, 84], axis=0)
    return [
        {
            "basic": [
                bnames,
                bnames2,
                bsamples,
                MLE_bsample,
                brange,
                bmedian,
                blow - bmedian,
                bhigh - bmedian,
            ],
            "clouds": [
                cnames,
                cnames2,
                csamples,
                MLE_csample,
                crange,
                cmedian,
                clow - cmedian,
                chigh - cmedian,
            ],
            "gases": [
                gnames,
                gnames2,
                gsamples,
                MLE_gsample,
                grange,
                gmedian,
                glow - gmedian,
                ghigh - gmedian,
            ],
            "T-P": [
                tnames,
                tnames2,
                tsamples,
                MLE_tsample,
                trange,
                tmedian,
                tlow - tmedian,
                thigh - tmedian,
            ],
            "custom": [
                knames,
                knames2,
                ksamples,
                MLE_ksample,
                krange,
                kmedian,
                klow - kmedian,
                khigh - kmedian,
            ],
        },
        np.shape(samples)[-1],
    ]


def TP_exponential_linear(
    exponent, logP_iso, T0, T_mid, T_max, logP=np.linspace(-4, 2.3, 1000), logP0=-4
):
    T_frac = exponent * (T_mid - T0) / (T_max - T_mid)
    logP_mid = (logP0 + T_frac * logP_iso) / (1 + T_frac)

    alpha2 = (logP_iso - logP_mid) / (T_max - T_mid)
    alpha1 = exponent * (T_mid - T0) ** (1 - (1 / exponent)) * alpha2
    T = np.where(logP > logP0, T0 + ((logP - logP0) / alpha1) ** exponent, T0)
    T = np.where(logP >= logP_mid, T_mid + (logP - logP_mid) / alpha2, T)
    T = np.where(logP >= logP_iso, T_max, T)

    return T
