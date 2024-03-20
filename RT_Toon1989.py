import numpy as np
from numpy.typing import ArrayLike

from useful_internal_functions import interleave

c = 2.99792458e10  # in CGS
hc = 1.98644568e-16  # in CGS
hc_over_k = 1.98644568 / 1.38064852

MAXIMUM_EXP_FLOAT = 35

stream_cosine_angles = np.array(
    [
        0.0446339553,
        0.1443662570,
        0.2868247571,
        0.4548133152,
        0.6280678354,
        0.7856915206,
        0.9086763921,
        0.9822200849,
    ]
)

stream_weights = np.array(
    [
        0.0032951914,
        0.0178429027,
        0.0454393195,
        0.0791995995,
        0.1060473494,
        0.1125057995,
        0.0911190236,
        0.0445508044,
    ]
)


###############################################################################
############################ Main callable function. ##########################
###############################################################################
def RT_Toon1989(
    wavelengths_in_cm: ArrayLike,
    temperatures_in_K: ArrayLike,
    optical_depth_per_layer: ArrayLike,
    single_scattering_albedo: ArrayLike,
    scattering_asymmetry: ArrayLike,
    stream_cosine_angles: ArrayLike = stream_cosine_angles,
    stream_weights: ArrayLike = stream_weights,
):
    thermal_intensity, delta_thermal_intensity = thermal_intensity_by_layer(
        temperatures_in_K, wavelengths_in_cm
    )

    terms_for_DSolver = calculate_terms_for_DSolver(
        optical_depth_per_layer,
        single_scattering_albedo,
        scattering_asymmetry,
        thermal_intensity,
        delta_thermal_intensity,
    )

    terms_for_DTRIDGL = DSolver_subroutine(*terms_for_DSolver)

    xki_terms = DTRIDGL_subroutine(*terms_for_DTRIDGL)

    total_flux = calculate_flux(
        optical_depth_per_layer,
        single_scattering_albedo,
        scattering_asymmetry,
        thermal_intensity,
        delta_thermal_intensity,
        stream_cosine_angles,
        stream_weights,
        xki_terms,
    )
    return total_flux


###############################################################################


def blackbody_intensity_by_wavelength(temperature_in_K, wavelength_in_cm):
    T, wave = temperature_in_K, wavelength_in_cm
    return (2 * hc * c / wave**5) / (np.exp(hc_over_k * (wave / T)) - 1)


def thermal_intensity_by_layer(
    temperatures_in_K: ArrayLike, wavelengths_in_cm: ArrayLike
):
    wavelength_grid, temperature_grid = np.meshgrid(
        wavelengths_in_cm, temperatures_in_K
    )
    thermal_intensity_bin_edges = blackbody_intensity_by_wavelength(
        temperature_grid, wavelength_grid
    )

    # mean across each layer bin
    thermal_intensity = (
        thermal_intensity_bin_edges[:, :-1] + thermal_intensity_bin_edges[:, 1:]
    ) / 2
    # change across each layer bin
    delta_thermal_intensity = (
        thermal_intensity_bin_edges[:, 1:] - thermal_intensity_bin_edges[:, :-1]
    )

    return thermal_intensity, delta_thermal_intensity


def calculate_terms_for_DSolver(
    optical_depth_per_layer: ArrayLike,
    single_scattering_albedo: ArrayLike,
    scattering_asymmetry: ArrayLike,
    thermal_intensity: ArrayLike,
    delta_thermal_intensity: ArrayLike,
    mu_1: float = 0.5,  # This is mu_1 in Toon et al. 1989
):

    number_of_wavelengths, number_of_layers = np.shape(optical_depth_per_layer)
    tau = optical_depth_per_layer
    w0 = single_scattering_albedo
    g = scattering_asymmetry
    tbfrac = 1  # INCOMPLETE IMPLEMENTATION
    # tbase = getT(hmin)       # INCOMPLETE IMPLEMENTATION
    thermal_intensity_at_TOA = thermal_intensity[0] - delta_thermal_intensity[0] / 2
    thermal_intensity_at_base = thermal_intensity[-1] + delta_thermal_intensity[-1] / 2

    alpha = np.sqrt(
        (1 - w0) / (1 - w0 * g)
    )  # sqrt( (1.-w0[i][j])/(1.-w0[i][j]*asym[i][j]) )
    lamda = alpha * (1 - w0 * g) / mu_1  # alpha[j]*(1.-w0[i][j]*cosbar[j])/mu_1
    gama = (1 - alpha) / (1 + alpha)  # (1.-alpha[j])/(1.+alpha[j])
    term = 1 / 2 / (1 - w0 * g)  # 0.5/(1.-w0[i][j]*cosbar[j])

    dti_x_term = delta_thermal_intensity * term
    dti_x_tau = delta_thermal_intensity * tau

    cpm1 = thermal_intensity + dti_x_term
    cp = cpm1 + dti_x_tau
    cmm1 = thermal_intensity - dti_x_term
    cm = cmm1 + dti_x_tau

    lamda_x_tau = np.clip(lamda * tau, a_min=None, a_max=MAXIMUM_EXP_FLOAT)
    ep = np.exp(lamda_x_tau)

    tautop = tau[:, 0]
    btop = (
        1 - np.exp(-tautop / mu_1)
    ) * thermal_intensity_at_TOA  # (1. - exp(-tautop/mu_1))*blackbodyL(tprof[0],wavelength)
    bsurf = thermal_intensity_at_base  # blackbodyL(tbase,wavelength)
    bottom = bsurf + delta_thermal_intensity[:, -1] * (
        mu_1 / tbfrac
    )  # Equivalent to multiplying taulayer by tbfrac.

    return cp, cpm1, cm, cmm1, ep, btop, bottom, gama


def DSolver_subroutine(cp, cpm1, cm, cmm1, ep, btop, bottom, gama, rsf=0):
    # rsf is "surface" reflectivity, can set to zero
    # Surface reflectivity should be zero for emission.
    # (It is used in the tridiagonal matrix in the bottom layer.)

    # DSolver subroutine to compute xk1 and xk2.
    # Computes a,b,c,d coefficients first, top to bottom
    # Then as and ds, *bottom to top*
    # Then xk coefficients, top to bottom
    # af, bd, cd, and df appear to be A_l, B_l, D_l, and E_l in Toon et al.
    # xk1 and xk2 appear to be Y_1n and Y_2n in Toon et al.
    # However, these do not match their formulae.
    top_layer = [slice(None), slice(0)]
    # bottom_layer = [slice(None), slice(-1)]
    upper_edges = [slice(None), slice(None, -1)]
    lower_edges = [slice(None), slice(1, None)]

    e1 = ep + gama / ep
    e2 = ep - gama / ep
    e3 = gama * ep + 1 / ep
    e4 = gama * ep - 1 / ep

    af_top = 0
    bf_top = gama[top_layer] + 1
    cf_top = gama[top_layer] - 1
    df_top = btop - cmm1[top_layer]

    # odd indices
    odd_afs = (e1[upper_edges] + e3[upper_edges]) * (
        gama[lower_edges] - 1
    )  # (e1[nn]+e3[nn])*(gama[nn+1]-1.)
    odd_bfs = (e2[upper_edges] + e4[upper_edges]) * (
        gama[lower_edges] - 1
    )  # e2[nn]+e4[nn])*(gama[nn+1]-1.)
    odd_cfs = 2 * (1 - gama[lower_edges] ** 2)  # 2.*(1.-gama[nn+1]*gama[nn+1])
    odd_dfs = (
        gama[lower_edges] - 1
    ) * (  # (gama[nn+1]-1.)*(cpm1[nn+1]-cp[nn]) + (1.-gama[nn+1])*(cm[nn]-cmm1[nn+1])
        (cpm1[lower_edges] - cp[upper_edges]) + (cmm1[lower_edges] - cm[upper_edges])
    )

    # even indices -- NOTE: even and odd have been switched from the
    # Fortran code and Toon et al. due to fencepost effects.
    even_afs = 2 * (1 - gama[upper_edges] ** 2)  # 2.*(1.-gama[nn]*gama[nn])
    even_bfs = (e1[upper_edges] - e3[upper_edges]) * (
        gama[lower_edges] + 1
    )  # (e1[nn]-e3[nn])*(1.+gama[nn+1])
    even_cfs = odd_afs  # (e1[nn]+e3[nn])*(gama[nn+1]-1.)
    even_dfs = e3[upper_edges] * (cpm1[lower_edges] - cp[upper_edges]) + e1[
        upper_edges
    ] * (
        cm[upper_edges] - cmm1[lower_edges]
    )  # e3[nn]*(cpm1[nn+1]-cp[nn]) + e1[nn]*(cm[nn]-cmm1[nn+1])

    af_base = e1[:, -1] - rsf * e3[:, -1]  # e1[nlayershort-1] - rsf*e3[nlayershort-1]
    bf_base = e2[:, -1] - rsf * e4[:, -1]  # e2[nlayershort-1] - rsf*e4[nlayershort-1]
    cf_base = 0
    # note: original C++ version says bsurf, but was called with bottom
    df_base = (
        bottom - cp[:, -1] + rsf * cm[:, -1]
    )  # bottom - cp[nlayershort-1] + rsf*cm[nlayershort-1]

    afs = np.stack([af_top, interleave(odd_afs, even_afs), af_base], axis=-1)
    bfs = np.stack([bf_top, interleave(odd_bfs, even_bfs), bf_base], axis=-1)
    cfs = np.stack([cf_top, interleave(odd_cfs, even_cfs), cf_base], axis=-1)
    dfs = np.stack([df_top, interleave(odd_dfs, even_dfs), df_base], axis=-1)

    return afs, bfs, cfs, dfs


# End of DSolver subroutine.


def DTRIDGL_subroutine(afs, bfs, cfs, dfs):
    # DTRIDGL subroutine to compute the necessary xki array
    # This matches the algorithm in Toon et al.
    af_base = afs[-1]
    bf_base = bfs[-1]
    df_base = dfs[-1]

    as_base = af_base / bf_base  # as[nl2-1] = af[nl2-1]/bf[nl2-1]
    ds_base = df_base / bf_base  # ds[nl2-1] = df[nl2-1]/bf[nl2-1]

    as_terms = np.empty_like(afs, dtype=np.float64)
    as_terms[-1] = as_base
    ds_terms = np.empty_like(afs, dtype=np.float64)
    ds_terms[-1] = ds_base

    twice_number_of_layers = np.shape(afs)[-1]

    for half_layer in reversed(range(twice_number_of_layers)):
        xx = 1 / (bfs[half_layer] - cfs[half_layer] * as_terms[half_layer])
        as_terms[half_layer - 1] = afs[half_layer] / xx
        ds_terms[half_layer - 1] = (
            dfs[half_layer] - cfs[half_layer] * ds_terms[half_layer]
        ) * xx

    xki_terms = np.empty_like(ds_terms)
    xki_terms[0] = ds_terms[0]

    for half_layer, (as_term, ds_term) in enumerate(zip(as_terms[1:], ds_terms[1:])):
        xki_terms[half_layer + 1] = ds_term - as_term * xki_terms[half_layer]

    return xki_terms


# End of DTRIDGL subroutine


def calculate_flux(
    optical_depth_per_layer: ArrayLike,
    single_scattering_albedo: ArrayLike,
    scattering_asymmetry: ArrayLike,
    thermal_intensity: ArrayLike,
    delta_thermal_intensity: ArrayLike,
    stream_cosine_angles: ArrayLike,
    stream_weights: ArrayLike,
    xki_terms: ArrayLike,
    mu_1: float = 0.5,  # This is mu_1 in Toon et al. 1989
):
    tau = optical_depth_per_layer
    w0 = single_scattering_albedo
    g = scattering_asymmetry
    thermal_intensity_at_base = thermal_intensity[-1] + delta_thermal_intensity[-1] / 2

    bsurf = thermal_intensity_at_base
    number_of_layers = np.shape(tau)[-1]

    # NOTE: there was a line in the original C++ for loop (index n3):
    # if(xk2[n3]!=0. && fabs(xk2[n3]/xk[2*n3] < 1.e-30)) xk2[n3] = 0.;
    # but note xk was only initialized, so all xk would be zero at this step?
    even_xki_terms = xki_terms[0::2]
    odd_xki_terms = xki_terms[1::2]
    xk1_terms = even_xki_terms + odd_xki_terms
    xk2_terms = even_xki_terms - odd_xki_terms

    # These are calculated just as they are in the setup function.
    # My goal is to decouple the RT components as much as possible, which leads
    # to this bit of redundant calculation (there's probably a better way!).
    alpha = np.sqrt(
        (1 - w0) / (1 - w0 * g)
    )  # sqrt( (1.-w0[i][j])/(1.-w0[i][j]*asym[i][j]) )
    lamda = alpha * (1 - w0 * g) / mu_1  # alpha[j]*(1.-w0[i][j]*cosbar[j])/mu_1
    lamda_x_tau = np.clip(lamda * tau, a_min=None, a_max=MAXIMUM_EXP_FLOAT)

    # These are the variables that are used to compute the flux.
    # They are all functions of the e-coefficient and blackbody fluxes via the matrix solver.
    gg_terms = xk1_terms * 2 * np.pi * w0 * (1 + (g * alpha)) / (1 + alpha)
    hh_terms = xk2_terms * 2 * np.pi * w0 * (1 - (g * alpha)) / (1 + alpha)

    blackbody_scattering_term = delta_thermal_intensity * (
        mu_1 * (w0 * g) / (1 - w0 * g)
    )
    alpha1 = 2 * np.pi * (thermal_intensity + blackbody_scattering_term)
    alpha2 = 2 * np.pi * delta_thermal_intensity

    stream_cosine_angles = np.expand_dims(
        stream_cosine_angles, axis=tuple(range(1, w0.ndim + 1))
    )
    stream_weights = np.expand_dims(stream_weights, axis=tuple(range(1, w0.ndim + 1)))

    epp_terms = np.exp(lamda_x_tau)
    em1_terms = np.exp(-lamda_x_tau)
    em2_terms = np.exp(-tau / stream_cosine_angles)
    em3_terms = em1_terms * em2_terms

    fpt_base = (
        2 * np.pi * (bsurf + delta_thermal_intensity[:, -1] * stream_cosine_angles)
    )
    fpt_terms = np.empty_like(delta_thermal_intensity, dtype=np.float64)
    fpt_terms[-1] = fpt_base

    for layer in reversed(range(number_of_layers)):
        fpt_terms[:, layer - 1] = (
            fpt_terms[:, layer] * em2_terms[:, layer]
            + gg_terms[:, layer]
            / (lamda[:, layer] * stream_cosine_angles - 1)
            * (epp_terms[:, layer] * em2_terms[:, layer] - 1)
            + hh_terms[:, layer]
            / (lamda[:, layer] * stream_cosine_angles + 1)
            * (1 - em3_terms[:, layer])
            + alpha1[:, layer] * (1 - em2_terms[:, layer])
            + alpha2[:, layer]
            * (stream_cosine_angles * (em2_terms[:, layer] - 1) + tau[..., layer])
        )

    total_flux = stream_weights * fpt_terms[..., ::-1]
    return total_flux
