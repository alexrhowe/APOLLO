USER GUIDE TO THE APOLLO ATMOSPHERE RETRIEVAL CODE

APOLLO v0.11.3 test build
30 September 2020

Developers:
Alex R. Howe, NASA Goddard Space Flight Center
Arthur Adams, University of Michigan, Department of Astronomy

Included files:
Apollo.tar
	README.txt
	maketar.sh
	APOLLO source code
	APOLLO subroutines and setup code
	MakeHaze source code
	Plot Apollo source code
	Plot Converge source code
	Plot Opac source code
	JWST filter list
	JWST filter functions
	Example files
Apollo.opacs.tar.gz
	Molecular cross section tables
	Aerosol optical properties tables
	Aerosol optical cross section tables

CONTENTS
1. Description
2. Functionality
3. Installation
4. Running
5. Input File Format: Settings
6. Input File Format: Parameters
7. Spectrum File Format
8. Cross Section Table Formats
9. Output File Formats

1. DESCRIPTION

APOLLO is an atmosphere code for extrasolar planets designed to fit
forward models to both transit and emission spectra and to retrieve atmospheric
properties. APOLLO uses an MCMC algorithm to find the best fit to a spectrum
for a given parameter set. The code is designed to be modular enough to have
many options that may be turned on and off, including the type of
observations, a flexible molecular composition, multiple cloud prescriptions,
multiple temperature-pressure profile prescriptions, multiple priors, and
continuum normalization. It uses cgs units unless otherwise specified.

APOLLO is written in C++ and Python 3.8 using the Cython extension. (It should
still be backward-compatible with Python 2.) Anaconda is the recommended build
of Python given that scipy and astropy arei required. The emcee, schwimmbad, and
corner packages are needed in addition to Anaconda. The original build
environment includes:

Python 3.8.3
Anaconda 2020.07 for Python 3.8
emcee 3.0.2
schwimmbad 0.3.1
corner 2.1.0
Cython 0.29.21
GCC 4.9.4

The build environment needed to run APOLLO on a computing cluster will vary
with the particular cluster. You may need to load Anaconda and the C++
compiler, and possibly create a Python virtual environment. Contact your local
support team for assistance.

DEVELOPER'S NOTE: to convert from the local version to run on the Discover
cluster, do the following (some of this may apply to your local cluster):
1. In setup.py, change "gcc-mp-9" to "gcc" and "g++-mp-9" to "g++".
2. In Apollo.py, comment out "import corner".
3. Comment out the entire plotting section except for the definition of "first"
   at the beginning. (Ending at "fig3.savefig(fig3name)".)
4. Run setup.py.
5. Change the email address and account name in the run file.
6. Change the file paths in the example input file.

2. FUNCTIONALITY

PENDING TESTING

APOLLO version 0.11.2.1 is a test build. v0.10.4 was verified to be accurate and
stable with the one-stream radiative transfer algorithm for transit, secondary
eclipse, and direct imaging; for layered and 5-parameter T-P profiles; for
cloud-free models, an opaque cloud deck, a cloud slab of uniform particle size
and density, and a cloud of gaussian particle density and uniform particle size.

The two-stream radiative transfer algorithm is verified only for cloud-free and
opaque cloud deck eclipse spectra.

APOLLO is bundled with the following programs:

Plot Apollo, a standalone version of the plotting routines for APOLLO.

Plot Converge plots the motion of the walkers in each dimension to the screen.
This can be useful to troubleshoot if the walking is malfunctioning.

Plot Opac plots two sets of cross sections from different tables to the screen
for a particular pressure and temperature. This is meant to compare similar
tables to check for any any differences.

MakeHaze computes new haze opacity tables using Mie scattering from a table of
indices of refraction for the haze material. Note that MakeHaze is a fully C++
program, unlike the others.

3. INSTALLATION

To install APOLLO, unpack the Apollo.tar tarfile, and compile the code IN THE
SRC DIRECTORY using the Cython setup file:

python setup.py build_ext --inplace

NOTE: You will need to set the C and C++ compiler names in setup.py to the
compilers used by your system.

This script contains all of the setup needed to run APOLLO. However, it also
requires tables of molecular and haze particle cross sections for NIR and MIR
observations. Due to the size of these tables they are available only via USB
drive or FTP. Contact the developers for more information or create your own.

By default, these tables are read from a directory named Opacities in the
parent directory to your work directory. This allows multiple builds on the
code to be run in different directories. However, the memory allocation of a
cluster is often such that you may need to store the opacity tables in an
archival or scratch directory. The filepath to this directory is specified with
the Opacities parameter in the input file.

Use the included Makefile to build the MakeHaze program.

4. RUNNING

APOLLO has two operational modes: one to compute a planetary spectrum given a
set of parameters, and one to retrieve those parameters based on an observed
spectrum. For most purposes, the code is called with two command line arguments:
an input filename and the operational mode. Thus:

python Apollo.py <input file> <Spectrum OR Retrieval>

Note that the "Retrieval" mode will not work correctly from the command line
unless the code is placed in "Serial" mode.

If the operational mode is not specified, APOLLO will default to "Spectrum"
mode, and if no command line arguments are given, it will default to the
bundled example file, examples/example.resolved.dat.

For "Retrieval" mode, if a third command line argument, "Serial", is added, it
will override parallel operation and initiate a minimal-length serial run for
testing purposes. (This should take less than 10 minutes on a modern processor.)
APOLLO may also be run in serial mode by specifying "Parallel False" in the
input file.

For "Retrieval" mode, the default for APOLLO is parallel processing. For
typical computing clusters this can be done using a PBS, SLURM, or similar
script incorporating the following command:

mpirun python Apollo.py <input file> >> <scratch directory>/Apollo.log

Sending the log output to a scratch directory is often needed on a cluster to
reduce memory usage in the work directory. If it is not included, the log file
will be sent to the work directory. Other MPI arguments may also be needed,
depending on the cluster.

The example files include the parameters for GJ 570D for resolved spectra, and
for HD 189733b for transit and eclipse spectra.

NOTE: Different clusters may require different settings to run APOLLO
efficiently. For examples:

Discover at NCCS requires nodes=8 (28 cores per node), ntasks=216, and
cpus-per-task=1.

The old Flux cluster at U of M (which used PBS) required nodes=14 and ppn=16
(16 cores per node).

Great Lakes at U of M requires nodes=8 (36 cores per node), ntasks=216, and
cpus-per-task=1.

Plot Apollo is run from the command line as Plot.Apollo.py. It can replicate
the output plots for both APOLLO and Makespectrum. Plot Apollo takes up to 4
command line arguments:
1. Input file. This is a sample file or a spectrum file.
2. Object name. This will determine the output file names.
3. Data type. "Samples" for a samples file, which will create the retrieval
   plots such as the waterfall plots and will also print some of the 1-sigma
   confidence intervals to the screen. "Spectrum" for a spectrum file, which
   will create the spectrum plots only.
4. Reference spectrum. Required for SPECTRUM files only. The input (observed)
   spectrum to plot against the spectral fit. Residuals will be calculated by
   binning the spectral fit to the reference spectrum.
4. Confidence interval. Optional for SAMPLES files. Sets the bounds on the
   histograms such that the specified fraction of samples fall within the
   waterfall plots. Default: 0.99

Plot Converge is run from the command line as Plot.Converge.py. It plots the
motion of all of the walkers over a single dimension in each frame and will
show if APOLLO has converged. It takes a samples file name as a command line
argument. By default, it displays the plots in a window instead of saving them.

Plot Opac is run from the command line as Plot.Opac.py. It plots the molecular
cross sections from two specified tables at a particular pressure and
temperature. By default, it displays the plot in a window instead of saving it.
Plot Opac takes 4 command line arguments:
1. First cross section file name.
2. Second cross section file name.
3. Pressure in mbar.
4. Temperature in kelvin.
NOTE: Plot Opac rounds down to the nearest pressure-temperature point instead
of interpolating. This is because it is meant to compare the actual values in
the tables for similar tables, and interpolation could mask the actual
differences.

MakeHaze is run from the command line as ./MakeHaze. It takes up to 3 command
line arguments:
1. Opacities directory for both input and output files. Default: ../Opacities
2. Input file name, a table of real and complex indices of refraction versus
   wavelength. Default: blank.r.dat, a table with all indices of refraction set
   to 1, representing no refraction.
3. Output file name. Default: haze.table.dat

Note that MakeHaze.cpp should be able to compute >100 wavelengths x 161
particle sizes per minute. If it is taking significantly longer, there may
be a problem with the indices file.

5. INPUT FILE FORMAT: SETTINGS

The input file for APOLLO consists of two sections: a list of settings for a
particular configuration of the code, and a list of parameters for the forward
model to be fit. The settings may be listed in any order, but they all must
come before the parameters. The parameters must be grouped into blocks with
specific labels.

Defaults are set so that any parameters may be omitted from the input file. The
header line for the parameter section must always be included. In "Spectrum"
mode, the code will execute without any other lines present. However, in
"Retrieval" mode as least one variable parameter must be included to be
retrieved upon.

The SETTINGS block of the input file may contain any of the following settings.
(Settings not included in this list will be ignored.)

Object, 1 field:
1. The name of the target object. This field may be any string and is used only
   to set the name of the output file. Default: "example"

Mode, 1 field:
1. Type of spectrum. "Resolved" for resolved objects in emission. "Transit"
   for transit spectra. "Eclipse" for secondary eclipse. Default: "Resolved"

Parallel, 1 field:
1. Mode of operation. "Parallel" for parallel operation using the MPI Pool.
   "Serial" for serial operation from the command line. Default: "Parallel"
NOTE: The "Serial" command line argument overrides this setting.

Data, 3 fields:
1. Spectrum file name. This file contains the observations to be fit. Default:
   "example.obs.dat" (the bundled example file)
2. Polynominal fitting flag. If the codeword "Polyfit" is included, the code
   will normalize a subset of observations to a polynomial fit based on another
   subset of observations. The code selects the first increase in the wavenumber
   of the observations as the boundary between the two subsets, so all bands to
   be fit should be listed before the reference bands.

Convolve, 1 field:
1. Convolution factor. This is the factor by which to convolve the data to a
   lower resolution to account for the true resolving power of the instrument,
   or if it is to be blurred, e.g. to simulate a different instrument.
   Default: 1

Binning, 1 field:
1. Binning factor. This is the factor by which to bin the observations to a
   lower resolution. This may be useful because it is recommended for the
   forward model to be 5-10 times higher resolution than the observations for
   reliable results. It may also be useful to reduce the number of datapoints
   and mitigate underfitting. Default: 1

Degrade, 1 field:
1. Degradation factor. This is an integer factor by which to degrade the
   resolution of the forward model. This may be done to speed up the code
   proportionally. However, it should only be used with low-resolution
   observations, or for testing purposes. Default: 1

Prior, 1 field:
1. Type of prior to be used. Currently supported priors are "Uniform" and
   "Normal". Default: "Uniform"

N_Walkers, 1 field:
1. Number of walkers to be initialized by emcee. Default: 8*N, where N is the
   number of parameters.
NOTE 1: if this number is set too low, it will reset to the minimum allowed
   value of 2*N+2.
NOTE 2: N_Walkers must be even. If it is odd, the code will add 1.
NOTE 3: The "Serial" command line argument overrides this setting and sets it to
   the minimum value.

N_Steps, 1 field:
1. The number of steps for the MCMC algorithm to take. This must be enough for
   the walkers to converge and may require adjustment for any given object.
   Default: 30000
NOTE 1: The "Serial" command line argument overrides this setting and sets
   nsteps to 2.
NOTE 2: autocorrelation convergence checking is not yet implemented.

Star, 3 fields:
1. Stellar temperature in kelvin. Default: 5770 (Solar temperature)
2. Stellar radius in Solar radii. Default: 1.0
3. Semi-major axis of the planet in AU. Default: 1.0
NOTE: These parameters are not used in "Resolved" mode.

Star_Spec, 1 field:
1. Stellar spectrum file name. Currently used only in "Spectrum" mode. Defines a
   stellar spectrum to use in the noise model in place of a blackbody. If no
   file is specified, it will default to a blackbody spectrum.

Location, 3 fields:
1. Distance in parsecs. Default: 10.0
2. Right ascension in decimal degrees. Default: 0.0
3. Declination in decimal degrees. Default: 0.0

Minimum_Mass, 1 field:
1. Minimum mass of the planet in Jupiter masses. This sets an additional prior
   on the radius and gravity of the planet. This is intended to prevent
   atmospheres too extended for the code to find a structural solution.
   Default: 0.5

Tables, 2 fields:
1. Set of opacity tables to use for computing the spectrum. For example,
   listing "nir" will use files of the form "<molecule>.nir.dat".
   NOTE: The default for this field is for the code to automatically select
   "nir", "mir", or "wide" depending on the wavelengths covered by the observed
   spectrum.
2. Set of opacity tables to use for computing the effective temperature. The
   effective temperature must be computed over a wide wavelength range, so a
   much lower spectral resolution is needed to do it efficiently.
   Default: "lores"

Standard tables:
"nir": 0.6-5.0 microns, R=10,000.
"mir": 5.0-30.0 microns, R=10,000.
"wide": 0.6-30.0 microns, R=10,000.
"lores": 0.6-30.0 microns, R=200.
NOTE: High-resolution tables with R=50,000 are available on an as-needed basis.

Pressure, 2 fields:
1. Minimum log-pressure in bars over which to compute the radiative transfer.
   Default: -3.0
2. Maximum log-pressure in bars over which to compute the radiative transfer.
   Default: 2.5

Gray, 2 fields:
1. Boolean switch to use a gray atmosphere model. Default: False
2. Temperature of the gray atmosphere model. Default: 1500 K.

Vres, 1 field:
1. Vertical resolution: number of layers to use for the radiative transfer
   calculation. Default: 71
NOTE: This applies to both the "Layers" and "Parametric" T-P profiles.

Streams, 1 field:
1. Which radiative transfer function to use. 1 for 1-stream, 2 for 2-stream.
   Default: 1
[PENDING]
NOTE: the 2-stream algorithm is currently implemented only for cloudless
   atmospheres and opaque cloud decks in emission.

Output_Mode, 1 field:
1. Output spectroscopic or photometric mode. Used only in "Spectrum" mode. If
   one of JWST's spectroscopic modes is specified, Makespectrum will create a
   spectrum with noise using a pipeline for that mode in addition to the
   full-resolution spectrum. If a file containing a filter throughput function
   is specified, Makespectrum will compute the weighted mean flux through the
   filter. Default: none
NOTE: These filter functions must have two columns:
1. Wavenumber in cm^-1
2. Fractional throughput.

Supported JWST modes:
NIRISS
NIRCam
G140H-F070LP (NIRSpec)
G140H-F100LP (NIRSpec)
G235H-F170LP (NIRSpec)
G395H-F290LP (NIRSpec)
G140M-F070LP (NIRSpec)
G140M-F100LP (NIRSpec)
G235M-F170LP (NIRSpec)
G395M-F290LP (NIRSpec)
Prism-Clear (NIRSpec)
LRS_Slit (MIRI)
LRS_Slitless (MIRI)
MRS_A (MIRI)
MRS_B (MIRI)
MRS_C (MIRI)

Output, 3 fields:
1. Output directory for the samples file. This may need to be a scratch or
   archival directory on a cluster. This field is included only for the samples
   file in APOLLO because of the large sizes of these files. All other outputs
   automatically go to one of the local subdirectories included in Apollo.tar,
   either "modelspectra" or "plots". Default: "samples"
2. Short filename flag. If the codeword "Short" is included, the code will
   create output files with short names incorporating only the object name
   rather than the default names, which incorporate other properties of the
   particular calculation.
3. Print full output flag. If the codeword "Full" is included, the code will
   output the full sample array to the samples file instead of the last 10% of
   the array. This allows analysis of the burn-in phase of the run.
NOTE 1: The "Serial" command line argument automatically sets "Full" to True.
NOTE 2: "Short" and "Full" may be included in any order.

Opacities, 1 field:
1. Opacity files directory. This will often need to be a scratch or archival
   directory on a cluster. Default: "../Opacities"

6. INPUT FILE FORMAT: PARAMETERS

The PARAMETERS block of the input file MUST begin with the following header.
This header shows the format of the parameters section:

"Parameter   Initial   Mu    Sigma    Min     Max"

The model parameters are listed in blocks with headers to separate them in the
code's handling. There are 5 blocks currently implemented, which may be in any
order, and any individual block may be omitted. Each parameter is listed with
6 fields:
1. Parameter name.
2. Initial guess.
3. Mean of the prior. (Only used for normal priors.)
4. Standard deviation of the prior. If set to 0.0, the parameter is implemented
   as a constant and not retrieved upon. (Otherwise used for normal priors.)
5. Minimum cutoff of the prior.
6. Maximum cutoff of the prior.

The currently implemented blocks are as follows:

Basic
Bulk parameters of the planet. Currently supported parameters:
Rad: Radius in EARTH radii
RtoD: Log of radius-to-distance ratio
RtoD2U: Log of squared radius-to-distance ratio
   Default radius: 1 JUPITER radius (non-log)
Log(g): log(g) in cm/s^2. Default: 4.1 (Jupiter)
Cloud_Base: log-pressure in bars at the cloud deck that serves as the base of
   the atmosphere for the forward model. Default: 2.5
P_cl: alternate name for Cloud_Base.

NOTE 1: If the radius is not specified, the default is set to 1 Jupiter radius.
   However, the code will automatically normalize the spectrum to the total flux
   of the observations.
NOTE 2: APOLLO includes an internal prior that the mass must be less than 80
   Jupiter masses. It is important to set the initial radius and gravity so that
   the mass falls below this limit.
NOTE 3: Another internal prior requires the scale height of the atmosphere to be
   less than 5% of the radius of the planet. This is another measure to avoid
   atmospheres too extended to solve. This prior may need to be bypassed for
   "super-puff" planets.
NOTE 4: Cloud_Base/P_cl may be used for an opaque cloud deck in addition
   to another cloud model. It also works to model the limits imposed by
   refraction on transit spectra.

Gases, 1 field:
1. Filler gas. This gas will comprise the balance of the atmosphere after the
   listed gases in the model are accounted for, subtracting their abundances
   from 1. Any of the supported molecules may be used as the filler gas.
   Default: "h2only".
Parameters are a list of molecules included in the model in log-abundance.
Currently supported molecules:
h2 (H2 + He)
h2only
he
h-
h2o
ch4
co
co2
nh3
h2s
Burrows_alk (Burrows Alkali Metals)
Lupu_alk (Lupu Alkali Metals)
crh
feh
tio
vo
hcn
n2
ph3

Atm, 2 fields:
1. Type of T-P profile for the atmosphere. "Layers" for an arbitrary number of
   temperature points spaced evenly in log-pressure space. "Parametric" for a
   5-parameter analytic T-P profile from Line et al. (2012).
2. Override vertical resolution flag. If the code word "Verbatim" is included,
   the code will use a layered profile directly from the input file instead of
   the standard (usually higher) "Vres" vertical resolution.
NOTE: If the T-P profile is omitted entirely, the code will default to an
   isothermal atmosphere with a temperature of 1500 K.

Parameters for the "Layers" profile:
T#: an arbitrary number of temperature points in kelvin. 15 is recommended.
NOTE 1: A single T# point will result in an isothermal atmosphere.
NOTE 2: The default cross section tables have boundaries of 75 K and 4000 K,
   which are the default limits on the priors.
gamma: a smoothing parameter that penalizes wide swings in the T-P profile in
   the likelihood function. If included, the program will automatically
   implement a smoothed T-P profile. Defined using an Inverse Gamma Function as
   in Line et al. (2015).
NOTE 1: Although the peak in the distribution for gamma is set to 5.e-5, the
   upper bound must be set to 1000 (possibly more for high temperatures) because
   this is the approximate point where the roughness penalty approaches its
   asymptotic value.
NOTE 2: It may be necessary to set the initial guess for gamma high, >10, to
   avoid the code getting stuck in a local minimum with too much smoothing.

Parameters for the "Parametric" profile:
Tint: Effective temperature due to internal heat only.
log(kappaIR): Mean opacity in the infrared.
gamma1: gamma1 parameter from Line et al. (2012).
gamma2: gamma2 parameter from Line et al. (2012).
alpha: alpha parameter from Line et al. (2012).

NOTE: If "Parametric" is specified, and the wrong number of parameters are
listed, the code will halt with an error message.

Clouds, 2 fields:
1. Cloud model number. Default: 0
2. Aerosol material. Default: "None"

Currently supported aerosol materials:
H2SO4
Polyacetylene
Tholin
Corundum
Enstatite
Forsterite
Iron
KCl
Na2S
NH3Ice
Soot
H2OCirrus
H2OIce
ZnS

Currently supported cloud (aerosol) models:
0: No clouds (0 parameters).
1: Opaque cloud deck (1 parameter).
2: Slab cloud profile with uniform particle size (4 parameters).
3: Gaussian cloud profile with uniform particle size (4 parameters).

Model 1:
Cloud_Base: log-pressure in bars at the cloud deck that serves as the base of
   the atmosphere for the forward model.
P_cl: alternate name for "Cloud_Base".
NOTE 1: The cloud pressure level should NOT be specified in both "Basic" and
   "Clouds". The "Basic" setting overrides the "Clouds" setting, and unreliable
   results may occur.
NOTE 2: If the cloud pressure level is specified in "Clouds", it will NOT be
   included in the "Basic" waterfall plot.

Model 2:
Haze_abund: log-number density of haze particles in the cloud in cm^-3.
Haze_size: log-radius of uniform haze particles in microns.
Haze_minP: log-pressure at the cloud top in bars.
Haze_thick: thickness of the slab cloud expressed as the difference in dex of
   pressure from the cloud top to the cloud base.

Model 3:
Haze_Pabund: log-peak number density of haze particles in the cloud in cm^-3.
Haze_size: log-radius of uniform haze particles in microns.
Haze_meanP: log-pressure of the peak density of the cloud.
Haze_scale: standard deviation of the cloud density in dex of pressure.

End:
Statistical parameters needed for an MCMC fit:
deltaL: wavelength offset in nm to correct for calibration errors of the 
   observations. Default: 0.0
logf: log of error bar correction factor (multiplies the size of the error
   bars of the input spectrum). Default: 1.0

7. SPECTRUM FILE FORMAT

The spectrum file format has 6 columns and no header:
1. Wavenumber of short-wave end of wavelength bin in cm^-1
2. Wavenumber of long-wave end of wavelength bin in cm^-1
3. The model spectrum without noise added. Used only for outputs of specific
   spectroscopic modes.
4. Lower error bar of flux
5. Upper error bar of flux
6. For emission spectra, observed flux in erg s^-1 cm^-3; for transit spectra,
   the observed planet-star radius ratio. In specific mode output files, this
   column includes the computed noise.

The spectrum file must be ordered from highest wavenumber to lowest wavenumber
within each band because this can used to truncate the computed spectrum for
speed and to separate out normalization bands.

If the input file is give as one continuous band (with no gaps between bins),
the code will compute and fit the spectrum over the full width of the
wavelength range. If there are gaps between the bins, it will interpret the
spectrum as a set of wavelength bands or filters and fit only on those bands.
It will not compute the flux between the bands (except for Teff calculation) in
order to reduce run times.

It is also possible to do a fit with continuum normalization with APOLLO. To do
this, add the codeword "Polyfit" to the "Data" line in the input file. Then,
list the bands to be fit in order from highest to lowest wavenumber, and then,
list the bands to normalize upon from highest to lowest wavenumber. The code
will interpret everything after a band increases in wavenumber as the set to be
normalized upon. APOLLO will average each of these bands and do a polynomial
fit of the input spectrum with the same degree as the number of normalization
bands. The bands to be fit will be divided by the polynomial fit to normalize
out any pseudocontinuum that would complicate the fit.

8. CROSS SECTION TABLE FORMATS

The molecular cross section files should be in a subdirectory of the Opacities
directory called "gases."
Hazes (aerosols) have three subdirectories in the Opacities directory:
"absorption" for absorption cross sections, "scattering" for scattering cross
sections, and "asymmetry" for scattering asymmetry parameters.

The molecular cross section files created by MakeHaze consist of a header row
followed by a grid of cross sections. The header row contains 10 fields:
1. Number of pressure points. Packaged as 18.
2. Log-minimum pressure in bar. Packaged as -6.0.
3. Log-maximum pressure in bar. Packaged as +2.5.
4. Number of temperature points. Packaged as 36.
5. Log-minimum temperature. Packaged as 1.875061 (75 K).
6. Log-maximum temperature. Packaged as 3.599199 (4000 K).
7. Number of wavelength points.
8. Minimum wavelength in microns.
9. Maximum wavelength in microns.
10. Spectral resolution.

Each row of the grid corresponds to a single wavelength. Within each row, cross
sections are arranged from lowest to highest pressure as the major axis and
lowest to highest temperature as the minor axis.
NOTE: APOLLO uses the boundary values of the cross section tables if the
temperature or pressure values go outside the table.

The haze absoprtion, scattering, and asymmetry tables both have a header
consisting of 6 fields:
1. Number of frequency points (ideally the same as the cross section tables).
2. Minimum wavelength in microns.
3. Maximum wavelength in microns.
4. Number of particle sizes. Packaged as 161.
5. Minimum log-particle size in microns. Packaged as -4.0 (100 pm).
6. Maximum log-particle size in microns. Packaged as 4.0 (1 cm).

Each line of the table begins with the reference wavenlength in microns and is
followed by the the cross sections per particle or asymmetry parameters at each
particle size.

In addition to the pre-computed cross sections, the included MakeHaze.cpp
program can compute new cross section tables from a list of real and complex
indices of refraction of the haze material, which are included in the "indices"
subdirectory of the Opacities directory. The index of refraction tables have
three columns:
1. Wavelength in microns.
2. Real index of refraction.
3. Complex index of refraction.

NOTE: Atmosphere.cpp contains functions to compute haze (aerosol) cross sections
directly from the indices of refraction via Mie scattering, but these are not
supported in v0.11 because they much slower. Use MakeHaze to compute new haze
cross sections as heeded.

10. OUTPUT FILE FORMATS

APOLLO prints a list of MCMC samples to a file in the specified output
directory. Be sure there is enough memory to store the output in that
directory, especially on a computing cluster. The output file is about 6 MB per
1000 steps, using the standard format of printing only the last 10% of samples.

The first line of the output file is a header containing 3-6 fields:
1. Number of walkers.
2. Number of steps PRINTED (by default, 10% of the total).
3. Number of parameters LISTED (including derived parameters).

If a Layered T-P profile is used, it also includes the following. (These are
needed to set the bounds on the T-P plot.)
4. Number of layers.
5. Log-minimum pressure in bars.
6. Log-maximum pressure in bars.

The second line is the parameter list as specified in the input file, followed
by the derived parameters:
1. Mass in Jupiter masses.
2. C/O ratio.
3. [Fe/H] in dex.
4. Effective temperature.

Subsequent lines list all of the paramters for each sample. By default, only the last 10% of samples are printed, after the burn-in, to reduce file sizes.
However, in serial mode, all of the samples are printed. The full sample array
can be also printed by adding the flag "Full" to the "Output" line in the input
file.

APOLLO also produces a new input file for the object with the best fit
parameters from the retrieval. This makes it much easier to feed it back into
APOLLO to make new plots, such as specific mode plots. to plot it. This file
goes to the work directory and ends with "retrieved.dat"
NOTE: This file uses the input values for the statistical "End" parameters
instead of the best fit values because the best fit values cause APOLLO to
crash if fed back in directly.

In "Retrieval" mode, APOLLO produces five plots, which go to the "plots"
subdirectory:
1. A waterfall plot of the "basic" parameters of radius, gravity, and cloud
   base pressure (if included under Basic), plus the four derived parameters:
   mass, C/O ratio, metallicity, and effective temperature.
2. A waterfall plot of the molecular abundances.
3. A plot of the retrieved temperature-pressure profile with 1-sigma limits.
4. A plot of the best fit retrieved spectrum at the code's full resolution, with
   the observations and residuals.
5. A plot of the best fit retrieved spectrum binned to match the observations,
   with the observations themselves and the residuals.

In "Spectrum" mode, APOLLO produces the two spectrum plots, plus an optional
plot with simulated observations for a specific spectroscopic mode, along with
the un-noised model spectrum.

NOTE: By default, the waterfall plots set the boundaries such that 99% of the
samples fall within the histograms, depending on your results, you may need to
change this value. This can be done by remaking the plots with Plot Apollo and
using the optional final command line argument to specify a different
confidence interval.

In "Retrieval" mode, APOLLO creates two output spectrum files with the same
format as the input spectrum, which go to the "modelspectra" subdirectory. These
are the full-resolution best fit spectrum and the best fit spectrum binned to
match the observations. In "Spectrum" mode, it creates these two files and an
optional additional file for a specific spectroscopic mode spectrum.

If a photometric mode is specified instead of a spectroscopic mode, APOLLO will
print three fields to the command line:
1. The band center wavelength in cm^-1.
2. The integrated flux over the band in erg s^-1 cm^-2.
3. The predicted uncertainty of the observation in erg s^-1 cm^-2, accounting
   for the throughput function of the filter and instrumental noise.

Plot Apollo replicates all of the plots produced by APOLLO depending on the type
of input file. It also prints 1-sigma bounds on the four derived parameters to
the command line.

*******************************************************************************

TO-DO LIST

1. Reorganize the header files to handle public and private correctly.
2. Am I actually using the different upper and lower error bars?
3. Do I like the "norad" normalization?
4. Do I like the default "h2only" filler?
4. Add cloud models and parametric T-P profiles to Plot Apollo. Also consider
   other useful combinations like log(g) with abundances.
5. Add a check to Atmosphere.cpp to ensure the absorption and scattering cross
   section files have identical binnings. (Actually, I think I can just give
   them different table sizes since they do get interpolated independently.)
6. Give each cloud model a name to create a switch for them.
7. Add options to test for convergence on the fly for both emcee's
   autocorrelation and an additional stopping criterion to max_steps.
8. Add options for reflection spectra.
9. Consider combining Planet_layer and Planet_auto by having Apollo call
   getProfile directly.
10. Add an option to lnlike() to bypass the calculation of binmod and s2, and
    call addnoise() to find goodness of fit for specific JWST modes.
11. Consider adding more command line arguments to MakeHaze and checking for
    bad inputs.
12. Figure out how to run several spectra in a row from the command line in
    the Python environment without rereading the cross section tables.
13. Implement 2-stream radiative transfer for transits.
14. Let transits use an input stellar spectrum.
15. Consider changing some parameter names to be less confusing.
16. Allow the default hire table setting to detect the resolution of the
    observations and respond appropriately.
17. Edge effects cause errors when the specific mode is the same width as the
    mode spectrum.

CHANGE LOG

v0.11.3
Fixed an index error in the derived parameter calculation.
Updated Apollo.run for the current Python environment on Discover.
Fixed a bug that returned the wrong derived parameters due to a mismatch in the
shapes of "chain" and "blobs".
Removed the defunct "usebands" flag.
Updated the documentation.
Fixed one of the error statements in Plot Opac.
Revised getTauProf() with correct formulas for w0 and asym and more intuitive
calculations.
Corrected the definition of asf in MakeHaze.
Fixed several problems with the haze handling in Planet_layer.cpp and
Planet_auto.cpp.
Added the correct rstar to the specific mode plot calculation.
Fixed several problems with the handling of lores in getFlux().

v0.11.2
Fixed an error that reversed the function of BinBands.
Allowed Lupu_alk into longer-wavelength retrievals.
NOTE: the Burrows_alk cross sections include zero opacity past 20 microns.
Added a minimum mass parameter. Default: 0.5 Jupiters.
Added a prior that requires the scale height of the atmosphere (based on T_eff)
to be less than 0.05 times the radius. This prevents gravity from going so low
that the code is unable to find a T-P profile that converges.
NOTE: This prior may need to be adjusted based on testing.
NOTE: All of the mass, radius, and gravity priors may need to be adjusted for
"super-puff" planets that are unusually low-density.
Tested the 2-stream T_eff calculation.
Added the functionality to fix the value of any input parameter by setting
sigma = 0.
Fixed a fencepost error with the abundance array.
Updated and verified the functionality of all of the plots.
Added the wavelength calibration parameter back into the spectral plots.
Changed the names of the specific JWST modes to the correct grism+filter
combinations.
NOTE: slashes are replaced with hyphens and spaces with underscores.
Cleaned up the JWST mode plot and fixed the normalization.
Made the JWST mode plot compute the full-width spectrum so that all of the JWST
modes will be covered.
NOTE: this is now limited to Spectrum mode only because Retrieval mode doesn't
compute the full opacity table for efficiency.
Added the correct Right Ascension and Declination into the JWST mode calculator.
Fixed the normalization of the spectrum plots.
Changed the specific mode names to the correct JWST gratings and filters.
Restored the functionality for an empty parameter file.
NOTE: Retrieval mode does not work without at least one variable parameter.
Fixed a bug that set nsteps wrong in Serial mode if it wasn't in the input file.
Fixed the plotting of the derived parameters in the waterfall plot.  (I think
this went wrong because I copied it from Plot Apollo, which works differently.)
Made sure the plotting functions execute with a nearly-empty input file.
NOTE: The plots may look weird or empty is one or more blocks is empty.
Cleaned up and updated the output parameter file.
Converted the remaining .py files to be Python 3 compatible.
Updated transFlux to work with the lores table.
Fixed a bug in transTauProf and updated it to work with the lores table.
Moved databin and dataconv to separate lines in the input file.
Allowed the 'Polyfit' and 'Bands' flags to be in either order.
Clarified the usage of the 3 degrading parameters.
Changed the list order of the parameters and matched the example files to it in
order to make it a more logical progression.
Updated the example files.
Removed unused variables from the C++ files.
Removed unused variables from the header files.
Tightened up and clarified the default hires table setting.
Trimmed the ends of the input spectrum if it goes outside the hires table.
Switched to astropy's "convolve" function, which handles boundaries better than
numpy's function.
Removed the now-defunct "fixed abunds" code.
Cleaned up the comments.
Removed Makespectrum from the .tar file.
Removed the defunct "Spectrum" parameter.

v0.11.1.1
Updated Apollo.run to work with the new Python environment on Discover.
Fixed the calls to binBands, binModel, and get_Teff with the correct parameters
for the new formats.
Fixed the importation of the MPIPool.
Fixed the order of imports.
Changed the default alkalis in the example files to Lupu_alk.
Fixed a bug in the handling of the radius in the plotting section.
Fixed a bug in the specific mode plotting setup.

v0.11.1
Updated the SLURM run file for the new login server.
Fixed a bug in that caused the waterfall plotting in APOLLO to fail to read the
derived parameters.
Added databin and dataconv to the output parameter file (rather than setting
then to 1).
Changed the output filenames in Spectrum mode to prevent overwriting Retrieved
mode results.
Added the hires and lores tables to the output parameter file.
Fixed a bug that overrode the hires parameter in the input file.
Incorporated Arthur Adams's edits to allow for fixed gas abundances that are
not retrieved upon.
Fixed some missed variable name changes.
Incorporated Arthur Adams's edits to include separate absorption and scattering
cross sections for aerosols.
Added the "streams" variable to the "switches" array, per Arthur Adams.
Fixed the "watercirrus" and "waterice" aerosol names.
Corrected a bug that made the plotting program not use the full chain in Serial
mode.
Upgraded from emcee v2 to emcee v3.
At some point, I futurized all the print statements, but I don't remember doing
it. This should still work with both Python 2 and Python 3.
Fixed a bug that used the wrong radius in some cases for the spectrum plots.
Added basic garbage collection to reading in the opacity files.
Removed "streams" as a function parameter.
Added the "Star_Spec", "Spectrum", and "Output_Mode" input parameters from
Makespectrum to APOLLO.
NOTE: "datain" replaced "obsfile".
NOTE: "outmode" is used for the specific JWST mode plots, not yet implemented.
Modified MakeHaze to compute all three haze tables more easily.
Reduced the command line arguments for MakeHaze to just opacities directory and
haze name.
Implemented Arthur Adams's "absorption", "scattering", and "asymmetry"
subdirectories in "Opacities."
Moved the molecular cross sections to a "gases" subdirectory in "Opacities."
Moved the incides of refraction to an "indices" subdirectory in "Opacities."
The commented-out "exit" command before the Planet object is created
disappeared at some point. I added it back in.
Switched the order of the dataconv and databin parameters in the input file to
reflect the order in which they're used in the code.
Completely rebuilt the band handling to be simpler and more effective. The code
will not automatically detect discrete bands in the observations, separate
them, and individually convolve and bin them before reassembling them. It will
also compute the model spectrum in these bands only to save time and will
continue to do polynomial fitting when required.
Removed the "usebands" parameter, as it will now be done automatically.
Reorganized the To-Do List.
Fixed a bug that caused readopac to fail to read the header line of the
opacity files in some instances.
Got the polynomial fitting working with the new system.
Incorporated the specific mode plotting from Makespectrum. It currently does
not work because 

v0.11.0
Heavy overhaul of the main code to fix the convolution and binning system,
simplify the code structure, and reconcile several different builds currently
in use.
Moved the non-MCMC functions at the beginning of the code to a separate file,
ApolloFunctions.py.
Split the likelihood function into two functions, one to compute the forward
model spectrum and one to compute the likelihood. This removes the need for
separate APOLLO and Makespectrum codes.
Removed the Apollo.Makespectrum.py code and incorporated its functionality
into APOLLO with the 'Spectrum' command line argument.
Added the required 'Retrival' command line argument to do retrievals. (Default:
'Spectrum'.)
Added functionality for APOLLO to automatically plot the best fit spectrum
along with the other plots.
Replaced the various instances of binning and interpolating spectra with the
correct procedure of convolving the spectrum to the desired resolving power,
then binning it to the desired sampling resolution, which may be different.
Changed a lot of the variable names to make more sense and be more consistent
in context.
Fixed a bug where the Teff calculation would use the hi-res opacity tables.
Updated the Readme file.

v0.10.5
Changed the T_eff calculation to use the algorithm specified in the input file
instead of automatically 1-stream. (It had been that way for speed, but this
was not needed.)
Changed the 2-stream calculation with an opaque cloud deck to split the layer
in which the cloud tops reside--setting the temperature of the "bottom" layer
equal to T(hmin) and scaling the optical depth to the fraction of the layer
above the cloud tops.

v0.10.4.4
Corrected a couple of array sizes in Planet_auto.cpp.
Modified the getProfile subroutine in Planet_layer.cpp to cut off the planet's
radius at 5 R_J if low surface gravity makes the profile fail to converge.
Made minor corrections to the Emcee 3 conversion.

v0.10.4.3
Completely rewrote the GetScaOpac subroutine to include he, h-, h2s, hcn, n2,
and ph3, to make it simpler overall, and to make it easier to add new
molecules. Verified outputs.
NOTE: h2 and h2only have molecular weights that differ by 15%, which may affect
the retrievals.

v0.10.4.2
Changed wavelength to wavenumber in the specific mode output file of
Makespectrum.
Wrote the "Streams" parameter correctly in the output parameter file.
Changed the "binned" spectrum in Makespectrum from interpolating to actually
binning.
Added the wavelength offset to the binned spectrum and specific mode plots in
Makespectrum. (Not the fullres.)
Corrected the mean molecular weight for h2only, h2, and he filler gases.
Corrected the labelling of the base pressure in the waterfall plot.
Prevented the imin, imax, etc. arrays from going outside the limits of specwave
due to rounding errors.

v0.10.4.1
Corrected the right ascension in the output file.
Fixed the handling of mu for the cloud deck.
Removed extraneous print statements from Planet_layer.cpp.
Removed the logbeta line from examples.resolved.dat.
Removed the wrong definition of wstart from Planet_layer.cpp.
Fixed several problems in the 2-stream radiative transfer algorithm.

v0.10.4
Added the parameter name to the "Out of Bound" message.
Fixed a bug in setting the bounds on the opaque cloud deck level.
Corrected reading of the filler gas.
Clarified the "Mass out of Bound" message and removed the blank line.
Fixed the formula for the smoothness prior, more strongly penalizing large
values of gamma.
Fixed a bug that caused the walkers to get stuck for parameters where the
initial guess is zero.
Fixed a bug that was failing to account for the wavelength offset.
Fixed a bug that specifically set the wavelength offset to zero.
Changed the headers in the opacity files from millibars to bars.
Cleaned up the distinction between the cloud deck bounds and the integration
range.
Cleaned up the conversion from bars to cgs.
Allowed the cloud deck pressure level to be in both "Basic" and "Clouds".
NOTE: The declaration in "Clouds" overrides "Basic".
Added the opaque cloud deck to the two-stream calculation.
Fixed a bug in the data binning that was causing bins to stretch across gaps
between bands.
Fixed a bug that made the error bars on binned data too large.
Fixed a bug in the data binning that wrongly divided the flux by the bin width.
Added the code to better handle the wavelength bins at the ends of bands to
APOLLO from Makespectrum.
Made some minor improvements to the documentation.
Corrected right ascension handling in APOLLO.
Corrected minimum pressure handling in Makespectrum.

v0.10.3
Added the filler gas to the new input file produced by APOLLO.
Cleaned up and greatly simplified the spectroscopic mode calculation in Apollo
Makespectrum.
Added a plot of the spectroscopic mode to Apollo Makespectrum.
Added the remaining JWST modes to Makespectrum.
Added a check in Makespectrum to ensure the spectroscopic modes match the
computed wavelengths.
Fixed a bunch of problems related to getting the spectroscopic modes to work
right.
Simplified the handling of the specwave array in Makespectrum.
Harmonized the GetBins and GetBinnedSpec routines in APOLLO and Makespectrum.
Removed the minimum mass from the internal mass prior.
Changed the standard log(g) limits in the example files to 2.0 and 6.0. This
may exclude some "super-puffs" but will otherwise include all currently
observable planets without being unnecessarily wide.
Clarified the h- label on the command line when reading opacity tables.
Fixed a bug that caused the effective temperature calculation to crash with
cloud models 2 and 3.
Added an option to discard wavelengths between observation bands to reduce run
times. (Does not work together with the polynomial fitting.)

v0.10.2.1
Fixed the labeling of the x-axis in Plot Opac.
Fixed the Lupu_alk opacity files and verified them.

v0.10.2
Added the Plot Opac script to plot cross sections to check for differences.
Added documentation for Plot Opac.
Added the Output_Mode parameter that specifies a JWST mode or a filter function
for which to compute a spectrum or flux.
Updated and tested the JWST mode calculations in Makespectrum.
Tweaked the labelling in the spectrum fit plot.
Revised the To-Do list.
Recorded derived parameters for all samples, not just those with nonzero
probability.
Included all samples after the burn-in in the waterfall plots, correcting the
previous incorrect method.
Changed the radius output from Earth radii to Jupiter radii.
Corrected the notes on blob handling.
Added a section to Plot Apollo to print the uncertainties for the derived
parameters.
Added a "confidence interval" command line parameter to Plot Apollo.
Added "binned" and "fullres" spectrum files and best fit plots to the tar
archive.
Updates examples files with the "Opacities" parameter.
Added new Lupu_alk opacity tables (though still not fully verified).

v0.10.1
Finished a few missing haze cross section tables.
Added the pressure profile to the samples output file for layered T-P profiles
to allow them to be plotted correctly in Plot Apollo.
Heavy overhaul of the plotting subroutines in APOLLO to improve the waterfall
plots.
The probability function now computes the mass, C/O ratio, metallicity, and
temperature of the planet on the fly and passes them to the sampler as a blob.
The "Basic" waterfall plot now includes the four derived quantities with
correct labels. The "Gases" plot also has corrected labels.
Harmonized the plotting functions in Plot Apollo with those in APOLLO and
Makespectrum.
Added material to the README file about the various configurations needed to
run APOLLO on a cluster.
Modified the spectral fit plots to draw the observations as points with error
bars on top of the computed spectrum.
Created separate plots for the full-resolution computed spectrum and the
spectrum binned to the observations.
Added an option to choose a filler gas.
Set the limits in the corner plots to 99% confidence levels to cut out
outliers.
Corrected the hydrogen scattering opacity.
Fixed a bug that prevented the gas list from being passed to the Planet object
correctly.
Made the opacity table read-in subroutine print full filenames instead of
species.
Verified that the output of the code is identical to version 0.9.x when the
same format of opacity tables is used.
Added notes to the README file that include how to switch between Emcee 2 and
Emcee 3.

v0.10
Heavy overhaul of the opacity managing system with four sets of opacity tables:
1. NIR tables (nir), 0.6-5.0 microns at R=10000
2. MIR tables (mir), 5.0-30.0 microns at R=10000
3. Wide tables (wide), 0.6-30.0 microns at R=10000
4. Lo-res tables (lores), 0.6-30.0 microns at R=200, used to calculate
   effective temperature from bolometric flux.
APOLLO now reads in two sets of tables: the lo-res tables for effective
temperature, and automatically selecting the appropriate tables for the
spectral retrieval. The lo-res tables have been designed so that they add about
10% to the workload on a typical ground-based (e.g. JHK) spectrum.
This means that while APOLLO computes a hi-res spectrum over only the
wavelength range of the observations, it automatically computes a lo-res
spectrum over the full 0.6-30 micron range.
Effective tempearture is always computed with one stream.
The hi-res and lo-res tables must have the same dimensions in temperature and
pressure (as well as all molecules within each set).
The haze opacity is interpolated automatically instead of read directly from
the (internal) table, so it always uses the hi-res tables.
The new tables also extend to 4000 K instead of 3500 K.
The opacity tables now include a header listing the size and limits of the grid
in pressure, temperature, and wavelength space (assuming all of them are
logarithmically distributed).
Added H2-only, He, HCN, N2, and PH3 to the list of supported molecules.
The standard "h2" opacity tables are now H2+He instead of H2-only.
Added Mass, C/O, [Fe/H], and the corrected Teff to the "basic" waterfall plots
in Plot Apollo (direct APOLLO plots pending).
Fixed a bug that set the radius wrong in the best fit output file.
Increased the number of steps to 11 in Serial mode in case of rounding errors
to ensure that some samples output is produced.
Removed the "Spectrum" setting and wrange variable, which are replaced by the
new system.
Combined the integer "switch" parameters for MakePlanet into a single int array
argument.
Made the "Serial" override set printfull to True.
Stopped printing "Out of Bound" statements for the statistical "End" parameters
to make the output less unwieldy.
Set the default number of steps in Serial mode to 2. With printfull set to
True, this is all that is needed for testing.
Added gamma into the best fit output file where it had been previously omitted.
Added the C++ cutoff switch to Makespectrum.
Excluded samples that were out of bounds from the waterfall plots.
NOTE: This may not be the correct way to do it. I just needed to do it to
remove nonphysical Teff calculations. The full samples are still computed as
intermediate products. EDIT: It was wrong and has been changed in v0.10.2.
Cleaned up the output from the likelihood function test in APOLLO.
Added another cutoff switch to APOLLO, which causes it to test the likelihood
function once and then halt.
Reset the statistical "End" parameters in the best fit output file because the
parameters found by emcee prevented it from running. (Reason unknown. There
might be a stray log somewhere.)

v0.9.13.2
Cleaned up maketar.sh to prevent it from including extra files.
Corrected the H2 scattering cross section in Makespectrum.
Fixed the alkali metal opacities, which were wrongly multiplied by 10,000.
Removed the effective temperature calculation because it was a shortcut that
didn't work. A script to compute it correctly from the luminosity is pending.
Verified that the fencepost error in gamma appears to have resolved the hot top
of the atmosphere problem.
Fixed a couple bugs in the best fit output format.
Fixed a bug in the TP profile plotting in Plot Apollo.
Removed the minus sign in the best fit output format.
Fixed the order of indices in Plot Converge.
Changed the output of APOLLO to list the number of steps actually recorded in
the samples file. (This makes the Plot programs work better.)
Changed Plot Apollo to match the new format.
Corrected the labelling in Plot Apollo.

v0.9.13.1
Changed number of walkers to account for adding Teff to the table.
Added the parameter list to the samples file produced by APOLLO.
Made APOLLO create a new input file with the best fit parameters.
Added an option to print the full sample array instead of the last 10%.
Added code to Plot Apollo to print waterfall plots and improved the interface.
Added a line in the input file to specify the Opacities directory instead of
hard-coding it.
Added the Opacities directory to the arguments for MakeHaze.
Added the Plot Converge script to troubleshoot the movement of the walkers.
Modified APOLLO to print to the sample file step-by-step instead of all at once
to print partial data if a run crashes.
Used emcee blobs to correct and simplify the handling of Teff.
Corrected the handling of the radius prior.
Allowed the "Short" output flag to be in any position in Makespectrum.
Updated and cleaned up the README file.
Fixed a bug that caused walkers to initialize with gamma being too large.
Corrected the upper bound on gamma from 1 to 1000.
Corrected the header of the samples output file.
Corrected the radius handling in Plot Apollo and the waterfall plot in APOLLO
so that it uses the radius in Earth radii and not the cross section.

v0.9.13
Fixed a bug in the h-minus handling.
Added a switch to turn the polynomial fitting on and off.
Warning: polynomial fitting and binning of observations will not work at the
same time.
Changed the code to use a spectrum truncated by the extreme wavelength of the
observations so that the wavelength bins do not have to be in order for a
non-normalized spectrum.
Fixed a major bug in how the polynomial fitting is handled related to passing
it to the Planet object.
Streamlined the polynomial fitting process in general.
Created "short" versions of the eclipse and transit observation files that fall
entirely within the NIR range, and changed the example input files to match.
Fixed a bug that was causing the code not to evalute all of the observations if
they were out of order.
Made the code able to accept arbitrary observations and ignore those that are
outside of the model's wavelength range.
Changed the model binning to actual binning instead of just interpolation.
Added new and updated Rayleigh scattering cross sections.
Removed the extraneous GetFluxPerBin subroutine.
Changed the transit and eclipse examples back to the regular observation files.

v0.9.12.1
Fixed a bug that caused a crash when databin is set to 1.
Corrected output directory in example files.
Commented out print statements in the spectrum calculation loop to make the
output less unwieldy.
Fixed several bugs in the plot generation in APOLLO.
Verified APOLLO end to end in serial operation.
Fixed a fencepost error in the index for gamma.
Verified the observation binning.
Fixed an error in the output file names in Makespectrum.
Verified the short filenames.
Fixed a wavelength-wavenumber mismatch in the noise model.
Fixed output for the JWST modes.
Verified input stellar spectrum functionality.
Updated addNoise function calls in Makespectrum.
Added an example photometric filter option to Makespectrum with the
corresponding thruput function in the Filter directory.
Fixed a bunch of stuff to make the photometric filters work right.
Verified MIR functionality.

v0.9.12
Enhanced the code to work on MIR data (5-28 microns).
Added new MIR opacity tables for gases and aerosols.
Note: only uses one set of tables at a time. Automatically removes alkalis
from MIR calculations because the tables don't go that far.
Corrected documentation of the spectrum file format.
Added an optional parameter to import a stellar spectrum to the noise model.
Added an option to create short output file names.
Added an optional paramter to call AddNoise for a photometric filter with mode
number -1 and an optional parameter for the filter throughput file (required
for this mode). Returns the central frequency, flux, and noise for the whole
filter.
Created a source code directory.
Created a directory for the example files.
Corrected the example files to reflect the change in directory.

v0.9.11.1
New version to work on new computer and cluster at GSFC, and other updates.
Added bound-bound and bound-free H-minus opacity.
Modified setup.py to use the gcc-mp-9 and g++-mp-9 compilers used at GSFC.
Began working on Python 3 compatibility. Used futurize on Apollo.py, but it
appears that only the print statements have changed.
Fixed a fencepost error in defining the master opacity table.
Fixed a bug that was causing Teff to evaluate to zero.
Added the Plot.Apollo.py script, which bins and plots a retrieved spectrum
against observations.
Added the code for Plot.Apollo.py to Apollo Makespectrum.
Commented out the JWST filters in Apollo Makespectrum.
Added code for APOLLO to automatically create waterfall plots and a plot of
the T-P profile.
Added an option to APOLLO to bin down the observations.
Changed output file names to work in the new system and to be cleaner.
Updated example files.
Made MakeHaze able to take input and output files as command line arguments.

v0.9.11
Changed the hard-coded Opacities directory to match the cluster.
Removed defunct and unnecessary lines from the example files.
Specified cloud model 2 in the example files.
Added right ascension and declination to the example files, changing the
"Distance" line to "Location."
Updated the README file.

v0.9.10
Added a gaussian function to constants.cpp.
Changed the number of cloud parameters to a numerical cloud model in the input
file format.
Implemented a gaussian cloud profile model.
Fixed a bug in the haze opacity calculation in getTauProf.
Completely rebuilt the polynomial normalization based on experiments on the
Flux cluster.
Fixed various mistakes along the way.
Cleaned up the outputs throughout.
Made the section header comments more informative.
Made some of the error messages more informative.
Rearranged a few lines to better harmonize APOLLO and Makespectrum.
Changed the output of APOLLO to only output the last 10% of the length of the
chain to reduce file sizes. (The original is left, but commented out.)

v0.9.9
Limited ability of haze calculation to go outside the size range due to the
granularity of the table.
Fixed a bug in the haze calculation that caused a crash in the cloud layer was
narrower than one layer in the T-P profile.
Changed the haze in the example files from tholin to enstatite.
Fixed some format problems with the indicies of refraction files.
Doubled the width and resolution in particles size of the haze opacity tables.
Made changes to Atmosphere.cpp to prevent crashing due to rounding errors.
Added more decimal places to haze opacity tables to prevent rounding errors.
Added automatic flux normalization if no radius variable is included in the
input file
Added an option for polynomial continuum normalization.
Updated the README file with the current input format.
Added another tweak to Atmosphere.cpp to prevent crashing due to rounding
errors in the interpolation.

v0.9.8
Fixed handling of the gray atmosphere in Makespectrum.
Repaired the functionality of a parametric T-P profile in transit.
Removed a test output filename that I kept forgetting to remove.
Fixed a bug that caused an error if the pressure limits were too wide.
Updated example files without the Hazetype line.
Added an example file for transits with a parametric T-P profile.
Fixed a bug in getH that broke Planet_auto.cpp.
Fixed the noise models for eclipse and transit.
Harmonized minor differences between Planet_layer.cpp and Planet_auto.cpp.
Fixed handling of gray boolean input parameter.
Adjusted T-P profile to correct for variation in gravity with altitude.
Fixed the Teff calculation with an opaque cloud deck included.
Change the transit and eclipse example files to published data for HD 189733b.
Added an eclipse, parametric profile example file.
Fixed a sign error in the transit spectrum calculation.
Rechecked all of the example files by plotting the results.

v0.9.7
Changed the code to count the number of parameters instead of reading it from
an input file entry.
Fixed a bug where bounds were not set for teff0.
Made the number of layers of the radiative transfer calculation variable.
Allowed the layered T-P profile to be as short as 1 entry (isothermal).
Added a check to make sure the parametric T-P profile is the right length.
Added a default for a T-P profile being omitted.
Made the code check for the smoothing parameter automatically.
Updated example files with the new format.
Added checks to handle missing fields in the input files.
Allowed the parameter blocks to be in any order.
Fixed some indentation errors.
Fixed some bugs in the radius handling.
Added an option to compute a gray atmosphere.
Implemented the other available haze materials.
Simplified the molecule and haze name handling.
Removed the Hazetype line from the input file format.

v0.9.6
Cleaned up the example files a bit.
Made one-stream radiative transfer the default.
Added a try-except to make sure the input file actually exists.
Made a small change to the output filename format.
Fixed a bug in setting the bounds for the cloud deck.
Added a crude method in getTauProf in Planet_layer.cpp to compute T_eff.
Added a check in case the probability function returns NaN.
Fixed a bug that didn't send haze parameters to setParams.
Moved the Atmosphere initialization into the Planet constructor to speed up
execution.
Implemented computing the effective temperature of the atmosphere for emission
spectra and outputting it with the retrieved parameters.
Checked all the example input files for errors.
Cleaned up comments and extra print statements in the code.
Fixed the cloudless atmosphere base in APOLLO.

v0.9.5
Harmonized Apollo.pbs with the other example files.
Changed the binning from cubic interpolation to linear because cubic crashes
on flux.
Automatically cut the computed spectrum to the length of the input spectrum to
speed up execution.
Improved the 1-stream algorithm to split a layer at a cloud deck boundary.
Fixed the opaque cloud deck in emission and some stuff related to it.
Added the needed header line to the haze cross section tables.
Modified MakeHaze to add the header line to new tables.
Fixed the 4-parameter haze model in emission and some stuff related to it.
Fixed both clouds and haze in transit.
Added cloud and haze example files to the tar ball.
Updated parts of the README file that I'd been neglecting.

v0.9.4.1
Turned the low optical depth correction in getFlux back on.
Changed the getSpectrum function to include a switch between the 1-stream and
2-stream radiative transfer algorithms.
Cleaned up the mess the new getFlux was in.
Fixed the getFluxOneStream function.
Fixed a misnumbering of haze parameters in Planet_auto.
Fixed name of APOLLO output file.
Verified the functionality of the parametric T-P profile with the two-stream
algorithm.
Verified the two-stream algorithm for an isothermal atmosphere.
Corrected an error in the Exponential Integral solver in constants.cpp.
Fixed a duplicate use of mu in the one-stream radiative transfer.
Verified both radiative transfer functions against the analytic solution for
a gray atmosphere, and double-checked the normalization.

v0.9.4
Fixed a bug in the transit spectrum calculation.
Verified the functionality of transits.
Changed the radiative transfer algorithm from one-stream to two-stream, using
the short form of Mark Marley's algorithm used by Ben Burningham. This
algorithm uses a tridiagonal matrix solver.
Reversed order of profiles as needed for this algorithm.
Fixed bugs that came up in getH, getT, and getP that resulted from the
reversal of the profiles.
Added taulayer array that is needed because the new algorithm uses layer
optical depth instead of column optical depth.
Added the "verbatim" switch to use the T-P profile directly from the input
file rather than interpolating.

v0.9.3
Added the MakeHaze program and Makefile to make haze cross section tables from
indices of refraction.
Added index of refraction tables to the Opacities directory.
Fixed bug in MakeHaze that caused it to not account for particle size and
replaced the applicale haze files.
Fixed the Gaussian quadrature integration over solid angle to achieve more
accurate results.
Wrote a bash script to create the tarball of all needed files.
Verified the functionality of an opaque cloud deck.
Added a check to ensure the cloud deck is below the top of the model
atmosphere.
Removed the inconsisted usage of NIRISS_Throughput.dat from Filter.py and
removed NIRISS_Throughput.dat from the prerequisites list.
Harmonized the import lists between APOLLO and Apollo Makespectrum.
Harmonized the read in files sections between APOLLO and Apollo Makespectrum.
Fixed bugs in the handling of defaults.
Verified functionality of default settings except for the number of parameters
and the atmosphere profile.
Verified functionality of the Howe & Burrows (2012) haze model.
Changed the naming system to match the number of parameters in the input file.
Fixed a bug that set the output file name incorrectly in APOLLO.

v0.9.2
Set defaults for End parameters.
Added an override switch for serial-run tests.
Verified stability of the Rad and RtoD radius parameters.
Verified functionality of the 5-parameter T-P profile (Line et al., 2012).
* Note that Tint~300 and log(kappa)~0.3 provide a good fit to the example
spectrum.
Fixed a bug in GetScaOpac where some of the molecules were mislabelled.
Implemented separate Eclipse and Resolved object modes accounting for
stellar irradiation of close-in planets.
Added a check to correct an odd number of walkers.

v0.9.1
Added a field to the input file format for the output directory.
Added a switch to the input file format for serial vs. parallel operation.
Changed switches from On/Off to True/False.
Added a field to the input file format for number of walkers.
Added fields to the input file format for the minimum and maximum pressure over
which to calculate the radiative transfer.
Changed the order of the fields in the input file format for intuitiveness.
Changed labelling of alkali metal files to reflect their sources.
Made the reading of the input file more robust in the order of the
non-parameter settings.
Added default values for all of the input settings and the Basic parameters.
Added a bundled example input file and observations file, currently GJ 570D,
also allowing the code to be called with no arguments.
Fixed a bug in the mean molecular weight calculation that omitted the
high-temperature molecules.
Fixed a bug that could cause a crash if extraneous lines are inserted in the
parameter list.
Documented dependencies.
Bundled the spectrum creator.

WISH LIST

Add additional cloud models. This will probably be the top priority from this
    sub-list.
Add additional parametric T-P profiles.
Add additional molecules. H2-minus seems fairly important, but I'm not sure
    what else would be.
Add a chemical equilibrium code. Note: suppress reactions that are slower than
    the rotation and/or circulation timescale.
Make sure the 2-stream algorithm works with clouds. (Split the cloud
    absorption and scattering opacities, etc.) Note 1: need to set the optical
    depth to infinity if it is below an opaque cloud deck. Note 2: this is what
    Arthur Adams is working on.
Make opacity tables that extend 1% beyond the official edges of their coverage.

CHANGES NEED FOR 2-PART MODELS

NOTE: A better band-handling system could simplify this.
Implement a parameter for partial cloud coverage. Think about how many
    parameters I'll need to fit a cloudy and a clear model to a spectrum.
    But need to get single-component models working first.
Implement combined transit and eclipse retrieval. This may involve dayside and
    nightside atmosphere models, but the other part is that it will need to
    compute separate transit and emission spectra for each model.
Implement phase curve observations. Basically the same as partial cloud
    coverage: dayside and nightside models with the "partial coverage"
    parameter being the phase angle.
Implement separate evening and morning terminator models for transiting planets.

IDEAS FOR VERSION 2.

Break up the functions in APOLLO/Makespectrum and Planet_layer/Planet_auto into
    separate files to be imported to clean things up.
Create a single method to handle all of the spectroscopic mode calculations.
Combine Planet_layer and Planet_auto into a single Planet class.
Combine APOLLO and Makespectrum into a single code with options for its various
    functions.
Have the code recognize the specific parameters listed for particular
    parametric T-P profiles and cloud models and implement default values for
    missing parameters.
Have APOLLO create a retrieved spectrum file and best fit plot directly after
    an MCMC run, possibly with a separate spectrum-computing function and/or
    using blobs.
Clean up the wavelength v. wavenumber mess.
Change all labeling of "haze" to "aerosol."

DISCARDED IDEAS FROM VERSION 1

Make the settings not depend on being on separate lines. Dropped for the same
    reason Python uses whitespace. I've seen what happens when it's not
    enforced (*cough* CoolTLUSTY *cough*), and it's an unreadable mess.
Allow extrapolation of temperature outside the table. Using the boundary values
    is probably pretty safe in pressure space. In temperature, I've got it
    covered for 75 K to 4000 K, and that will cover the actual region being
    probed for everything from Saturn to the top of the L-dwarf range. Plus,
    I don't particularly trust the temperature extrapolation.
Turn the APOLLO and Makespectrum codes into functions in a single script. Even
    though so much of the two codes overlap, it doesn't make sense with them
    being *functionally* so different, with one being an MCMC that must be run
    on a cluster, while the other is a spectrum generator that can be run
    locally.
Add the functionality to compute haze scattering directly from the indices of
    refraction. It's so slow that it doesn't make sense to add it to APOLLO
    directly, especially since I've cleaned up the interface for MakeHaze. A
    lot of the issues that motivated this idea tie into the issue of how to
    handle the opacity files more generally.
Have APOLLO create a retrieved spectrum file and plot directly. This doesn't
    work because the spectrum calculation is hidden inside the lnlike function,
    and emcee isn't designed to modify it.
Add a default T-P profile for "Parametric" atmospheres. This doesn't work
    because the T-P profile parameters are tied in with the parameter list
    length. By default, the code throws an error message instead.
Add a symbol for the average error bars on the spectral plots. Discarded in
    favor of plotting all of the error bars.
Make the code able to handle tables of different sizes in T-P space.
    Currently, it assumes the table for each species is the same size. (This is
    mainly if I want to implement things like Exomol and/or mix and match
    them.) This is a problem because for efficient memory handling and spectrum
    calculation, the tables need to be the same size and not interpolated.
Add an average error bar size symbol to the JWST mode plot. Discarded because
    the error bars are too small to be visible.
