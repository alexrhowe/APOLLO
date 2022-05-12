USER GUIDE TO THE APOLLO ATMOSPHERE RETRIEVAL CODE

APOLLO v0.12 Beta
12 March 2021

NOTE: this public beta is generally stable for the cases described in
Section 2, but development is ongoing. We are running a full suite of test and
are hoping to produce a reliable v1.0 by the end of the year.

Molecular cross section files are not included with this release due to memory
constraints. Please contact the developer at alex.r.howe@nasa.gov for more
information, or build your own based on the format described in Section 8.

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
	Plot Haze source code
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
OpenMPI 4.0.3

Warning: in our tests, we have fonud that Intel MPI does NOT work.

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

APOLLO version 0.11.5 is a test build. v0.10.4 was verified to be accurate and
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

Plot Haze does the same for haze particle cross sections.

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

Plot Haze is run from the command line in the same way. It plots the particle
cross sections from two specified haze tables (which may be absorption,
scattering, or asymmetry parameter).
Plot Haze takes 3 command line arguments:
1. First cross section file name.
2. Second cross section file name.
3. Particle size in microns (rounded down to the nearest size on the table).

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

Mass_Limits, 2 fields:
1. Minimum mass of the planet in Jupiter masses. Default: 0.5
2. Maximum mass of the planet in Jupiter masses. Default: 80.0
NOTE: These set an additional prior on the radius and gravity of the planet. The
   minimum mass is intended in part to prevent atmospheres too extended for the
   code to find a structural solution.

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

Output_Mode, 2 fields:
1. Output spectroscopic or photometric mode. Used only in "Spectrum" mode. If
   one of JWST's spectroscopic modes is specified, Makespectrum will create a
   spectrum with noise using a pipeline for that mode in addition to the
   full-resolution spectrum. If a file containing a filter throughput function
   is specified, Makespectrum will compute the weighted mean flux through the
   filter. Default: none
2. Exposure time in seconds. Default: 1000.
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

Subsequent lines list all of the paramters for each sample. By default, only
the last 10% of samples are printed, after the burn-in, to reduce file sizes.
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
