USER GUIDE TO THE APOLLO ATMOSPHERE RETRIEVAL CODE

APOLLO v0.9.13.2rc Release Candidate
25 November 2019

Alex R. Howe

Included files:
Apollo.tar
	README.txt
	maketar.sh
	APOLLO source code
	Makespectrum source code
	MakeHaze source code
	Plot Apollo source code
	Plot Converge source code
	JWST filter list
	JWST filter functions
	Example files
Apollo.opacs.tar.gz
	Molecular cross section tables
	Aerosol optical properties tables
	Aerosol optical cross section tables

CONTENTS
1.  Description
2.  Functionality
3.  Installation
4.  Running
5.  Input File Format: Settings
6.  Input File Format: Parameters
7.  Spectrum File Format
8.  Molecular Cross Section Table Format
9.  Haze Table Format
10. Output File Formats

1. DESCRIPTION

APOLLO is an atmosphere code for extrasolar planets designed to fit
forward models to both transit and emission spectra and retrieve atmospheric
properties. APOLLO uses an MCMC algorithm to find the best fit to a spectrum
for a given parameter set. The code is designed to be modular enough to have
many options that may be turned on and off, including the type of
observations, a flexible molecular composition, multiple cloud prescriptions,
multiple temperature-pressure profile prescriptions, multiple priors, and
continuum normalization. It uses cgs units unless otherwise specified.

APOLLO is written in C++ and Python 2.7 using the Cython extension. (Python 3
is in the works.) Anaconda is the recommended build of Python given that scipy
is required. The emcee and corner packages are needed in addition to Anaconda.
The original build environment includes:

Python 2.7.13
Anaconda 4.4.0 for Python 2
emcee 2.2.1
corner 2.0.1
Cython 0.25.2
GCC 4.9.4

NOTE: the version of emcee is especially important. Emcee version 3 makes
significant changes an is not fully backward-compatible. In particular, the
MPI Pool functionality is forked to another module, schwimmbad, in Emcee 3.

The build environment needed to run APOLLO on a computing cluster will vary
with the particular cluster. You may need to load Anaconda and the C++
compiler, and possibly create a Python virtual environment. Contact your local
support team for assistance.

2. FUNCTIONALITY

APOLLO version 0.9.10 is verified to be accurate and stable with the one-stream
radiative transfer algorithm for transit, secondary eclipse, and direct
imaging; for layered and 5-parameter T-P profiles; for cloud-free models, an
opaque cloud deck, a cloud slab of uniform particle size and density, and a
cloud of gaussian particle density and uniform particle size.

Later version testing is pending.

The two-stream radiative transfer algorithm is verified only for cloud-free and
opaque cloud deck eclipse spectra.

APOLLO is bundled with the following programs:

Apollo Makespectrum program computes spectra of forward models complete with
noise levels for a model JWST pipeline. Makespectrum uses the same input file
format as APOLLO, although it uses a slightly different set of parameters. It
generates spectrum files that can be used as input "observations" in APOLLO.
Makespectrum currently produces only JWST spectroscopic and photometric modes.
However, it does output the full-resolution spectrum to compute noise levels
for model pipelines for other observatories.

Plot Apollo, a standalone version of the plotting routines for APOLLO and
Makespectrum.

Plot Converge plots the motion of the walkers in each dimension to the screen.
This can be useful to troubleshoot if the walking is malfunctioning.

MakeHaze computes new haze opacity tables using Mie scattering from a table of
indices of refraction for the haze material. Note that MakeHaze is a C++
program, unlike the others.

3. INSTALLATION

To install APOLLO, unpack the Apollo.tar tarfile, and compile the code IN THE
SRC DIRECTORY using the Cython setup file:

python setup.py build_ext --inplace

This script contains all of the setup needed to run both APOLLO and
Makespectrum. APOLLO is also bundled with tables of molecular and haze
particle cross sections for NIR and MIR observations, which are also necessary
to run the code. However, these tables total about 5 GB, so due to their size,
they are available only via USB drive or FTP.

NOTE: You will need to set the C and C++ compiler names in setup.py to the
compilers used by your system.

By default, these tables are read from a directory named Opacities in the
parent directory to the code's work directory. This allows multiple builds of
the code to be run in different directories. However, the memory allocation of
a cluster is often such that you will need to store the opacity tables in an
archival or scratch directory. The filepath to this directory is specified with
the Opacities parameter in the input file.

Use the included Makefile to build the MakeHaze program.

4. RUNNING

APOLLO is normally called with a single command line argument, which is the
name of the input file. If no input file is specified, it will default to the
bundled example file, example.resolved.dat. If a second command line argument,
"Serial", is added, it will override parallel operation and initiate a
minimal-length serial run for testing purposes. APOLLO may also be run in
serial mode if parallel mode is turned off in the input file.

In serial model, APOLLO may be run as a regular Python program:

python Apollo.py <input file> Serial

The default mode of APOLLO is parallel processing. On Flux, this can be done
using a PBS script incorporating the following command:

mpirun python Apollo.py <input file> >> <scratch directory>/Apollo.log

Sending the log output to a scratch directory is often needed on a cluster to
reduce memory usage in the work directory. If it is not included, the log file
will be sent to the work directory. Other MPI arguments may also be needed,
depending on the cluster.

The example files include the parameters for GJ 570D for resolved spectra, and
for HD 189733b for transit and eclipse spectra.

Apollo Makespectrum is run as a regular Python program and does not need MPI:

python Apollo.makespectrum.py <input file>

Plot Apollo is run from the command line as Plot.Apollo.py. It can replicate
the output plots for both APOLLO and Makespectrum. Plot Apollo takes 4 command
line paramters:
1. Input file. This is a sample file or a spectrum file.
2. Object name. This will determine the output file names.
3. Data type. "Samples" for a samples file, which will create the APOLLO plots.
"Spectrum" for a spectrum file, which will create the Makespectrum plot.
4. Reference spectrum. Required for spectrum files only. The input (observed)
spectrum to plot against the spectral fit. Residuals will be calculated by
binning the spectral fit to the reference spectrum.

Plot Converge is run from the command line as Plot.Converge.py. It plots the
motion of all of the walkers over a single dimension in each frame and will
show if APOLLO has converged. It takes a samples file name as a command line
argument.

MakeHaze is run from the command line as ./MakeHaze. It takes 3 optional
command line arguments:
1. Opacities directory for both input and output files. Default: ../Opacities
2. Input file name, a table of real and complex indices of refraction versus
wavelength. Default: blank.r.dat, a table with all indices of refraction set
to 1, representing no refraction.
3. Output file name. Default: haze.table.dat

Note that MakeHaze.cpp should be able to compute >100 wavelengths x 161
particle sizes per minute. If it is taking significantly longer, there may
be a problem with the input file.

5. INPUT FILE FORMAT: SETTINGS

The input file for APOLLO and Makespectrum consists of two sections: a list of
settings for a particular configuration of the code, and a list of parameters
for the forward model to be fit. The settings may be listed in any order, but
they all must come before the parameters. The parameters must be grouped into
blocks with specific labels.

Defaults are set so that any parameters may be omitted from the input file. The
header line for the parameter section must be included, but the code will
execute without any other lines present.

The input files are designed so that they can be used interchangeably with
APOLLO and Makespectrum, so some variables that are used for one are dummies in
the other.

The SETTINGS block of the input file may contain any of the following settings.
(Settings not included in this list will be ignored.)

Mode, 1 field:
1. Type of spectrum. "Resolved" for resolved objects in emission. "Transit"
for transit spectra. "Eclipse" for secondary eclipse. Default: "Resolved"

Parallel, 1 field:
1. Mode of operation. "Parallel" for parallel operation using the MPI Pool.
"Serial" for serial operation from the command line.
NOTE: The "Serial" command line argument overrides this setting.

Object, 1 field:
1. The name of the target object. This field may be anything and is used only
to set the name of the output file. Default: "example"

Star, 3 fields:
1. Stellar Temperature in kelvin. Default: 5770 (Solar temperature)
2. Stellar radius in Solar radii. Default: 1.0
3. Semi-major axis of the planet in AU. Default: 1.0
NOTE: These parameters are not used in the "Resolved" mode.

Location, 3 fields:
1. Distance in parsecs. Default: 10.0
2. Right ascension in decimal degrees. Default: 0.0
3. Declination in decimal degrees. Default: 0.0

Data, 3 fields:
1. Spectrum file name. This file contains the observations to be fit. Default:
"example.obs.dat" (the bundled example file)
2. Binning factor observations. This is an integer factor by which to bin down
the observations. This may be useful because it is recommended for the forward
model to be 5-10 times higher resolution than the observations for the most
reliable results. However, if the spectrum has gaps or varies significantly in
resolution, it should be set to 1. Default: 1
3. Polynomial fitting flag. If the codeword "Polyfit" is included, the code
will normalize a subset of observations to a polynomial fit based on another
set of observations. The boundary between the two sets is set where the
wavenumber of the observations in the input spectrum increases.

Star_Spec, 1 field:
1. Stellar spectrum file name. Used only in Makespectrum. Defines a stellar
spectrum to use in the noise model in place of a blackbody. If no file is
specified, it will default to a blackbody spectrum.

Degrade, 1 field:
1. Degradation factor. This is an integer factor by which to degrade the
spectral resolution of the forward model. Degrading the spectrum speeds up
the code proportionally, although the lower resolution may reduce the
reliability of the code. Default: 3
NOTE: This setting is hard-coded to 1 in Makespectrum because this program is
intended to produce a high-resolution model.

Spectrum, 1 field:
1. Spectral range. Used only in Makespectrum. "NIR" for near-infrared
(0.6-5.0 microns), "MIR" for mid-infrared (5-28 microns). Default: "NIR"
NOTE: APOLLO will use the correct wavelength range automatically, but it is
implemented only for data sets solely in one of the two ranges.

N_Walkers, 1 field:
1. Number of walkers to be initialized by emcee. Default: 8*N, where N is the
number of parameters.
NOTE: if this number is set too low, it will reset to the minimum allowed value
of 2*N+2.
NOTE: N_Walkers must be even. If it is odd, the code will add 1.
NOTE: The "Serial" command line argument overrides this setting and sets it to
the minimum value.

N_Steps, 1 field:
1. The number of steps for the MCMC algorithm to take. This must be enough for
the walkers to converge and may require adjustment for any given object.
Default: 30000
NOTE: The "Serial" command line argument overrides this setting and sets nsteps
to 10.
NOTE: autocorrelation convergence checking is not yet implemented.

Pressure, 2 fields:
1. Minimum log-pressure in bars over which to compute the radiative transfer.
Default: -3.0
2. Maximum log-pressure in bars over which to compute the radiative transfer.
Default: 2.5

Vres, 1 field:
1. Vertical resolution: number of layers to use for the radiative transfer
calculation. Default: 71
NOTE: This applies to both the "Layers" and "Parametric" T-P profiles.

Streams, 1 field:
1. Which radiative transfer function to use. 1 for 1-stream, 2 for 2-stream.
Default: 1
NOTE: the 2-stream algorithm is much slower than the 1-stream algorithm, by as
much as 2 orders of magnitude. It is currently implemented only for cloudless
atmospheres and opaque cloud decks in emission.

Prior, 1 field:
1. Type of prior to be used. Currently supported priors are "Uniform" and
"Normal". Default: "Uniform"

Gray, 2 fields:
1. Boolean switch to use a gray atmosphere model. Default: False
2. Temperature of the gray atmosphere model. Default: 1500 K.

Output, 3 fields:
1. Output directory for the samples file. This may need to be a scratch or
archival directory on a cluster. This field is included only for the samples
file in APOLLO because of the large sizes of these files. All other outputs
automatically go to one of the local subdirectories included in Apollo.tar,
either "modelspectra" or "plots". Default: "samples"
2. Short filename flag. If the codeword "Short" is included, the code will
create output files with short names incorporating only the object name rather
than the default names, which incorporate other properties of the particular
calculation.
3. Print full output flag. If the codeword "Full" is included, the code will
output the full sample array to the samples file instead of the last 10% of the
array. This allows analysis of the burn-in phase of the run.
NOTE: The order of "Short" and "Full" may be reversed for APOLLO.
NOTE: The "Short" flag may be in any position for Makespectrum.

Opacities, 1 field:
1. Opacity files directory. This will often need to be a scratch or archival
directory on a cluster. Default: "../Opacities"

6. INPUT FILE FORMAT: PARAMETERS

The PARAMETERS block of the input file MUST begin with a header beginning with
the string "Parameter". This header shows the format of the parameters section:

"Parameter   Initial   Mu    Sigma    Min     Max"

The model parameters are listed in blocks with headers to separate them in the
code's handling. There are 5 blocks currently implemented, which may be in any
order, and any individual block may be omitted. However, if a block is
included, general all of the parameters for that block MUST be included. Each
parameter is listed with 6 fields:
1. Parameter name.
2. Initial guess.
3. Mean of the prior. (Only used for normal priors.)
4. Standard deviation of the prior. (Only used for normal priors.)
5. Minimum cutoff of the prior.
6. Maximum cutoff of the prior.

The implemented values of each block are listed below.

Basic
Bulk parameters of the planet. Currently supported parameters:
Rad: Radius in Earth radii
RtoD: Log of radius-to-distance ratio
RtoD2U: Log of squared radius-to-distance ratio
Default radius: 1 Jupiter radius in linear space
Log(g): log(g) in cm/s^2. Default: 4.1
Cloud_Base: log-pressure in bars at the cloud deck that serves as the base of
the atmosphere for the forward model. Default: 2.5
P_cl: alternate name for Cloud_Base.

NOTE: If the radius is not specified, the default is set to 1 Jupiter radius.
However, the code will automatically normalize the spectrum to the total flux
of the observations.
NOTE: APOLLO includes an internal prior that the mass must be less than 80
Jupiter masses. It is important to set the initial radius and gravity so that
the mass falls below this limit.
NOTE: Cloud_Base/P_cl allows the inclusion of an opaque cloud deck in addition
to another cloud model. (This also works for refraction limits.)

Gases
List of molecules included in the model in log-abundance. Currently supported
molecules:
h2
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

NOTE: The H2 (actually H2 and He combined) abundance is set by the difference
of 1 minus the other abundances.
NOTE: Alkali metals are automatically removed from MIR calculations because
they have negligible opacity in that range.

Atm, 2 fields:
1. Type of T-P profile for the atmosphere. "Layers" for an arbitrary number of
temperature points spaced evenly in log-pressure space. "Parametric" for a
5-parameter analytic T-P profile taken from Line et al. (2012).
2. Override vertical resolution flag. If the code word "Verbatim" is included,
the code will use a layered profile directly from the input file instead of the
standard (usually higher) vertical resolution.
NOTE: If the T-P profile is omitted entirely, the code will default to an
isothermal atmosphere with a temperature of 1500 K.

Parameters for the "Layers" profile:
T#: an arbitrary number of temperature points in kelvin. 15 is recommended.
gamma: a smoothing parameter that penalizes wide swings in the T-P profile in
the likelihood function. If included, the program will automatically implement
a smoothed T-P profile. Defined using an Inverse Gamma Function as in
Line et al. (2015).
NOTE: Although the peak in the distribution for gamma is set to 5.e-5, the
upper bound must be set to 1000 because this is the approximate point where the
roughness penalty approaches its asymptotic value.

NOTE: A single T# point will result in an isothermal atmosphere.
NOTE: The default cross section tables have boundaries of 75 K and 4000 K,
which are the default limits on the priors.

Parameters for the "Parametric" profile:
Tint: Effective temperature due to internal heat only.
log(kappaIR): Mean opacity in the infrared.
gamma1: gamma1 parameter from Line et al. (2012).
gamma2: gamma2 parameter from Line et al. (2012).
alpha: alpha parameter from Line et al. (2012).

Clouds, 2 fields:
1. Cloud model number.
2. Haze type (optional).

Currently supported cloud (aerosol) models:
0: No clouds (0 parameters).
1: Opaque cloud deck (1 parameter).
2: Slab cloud profile with uniform particle size (4 parameters).
3: Gaussian cloud profile with uniform particle size (4 parameters).

Model 1:
P_cl: log-pressure in bars at the cloud deck that serves as the base of
the atmosphere for the forward model.

NOTE: this overrides any setting of P_cl/Cloud_Base in the Basic block.

Model 2:
Haze_abund: log-number density of haze particles in the cloud in cm^-3.
Haze_size: log-radius of uniform haze particles in microns.
Haze_minP: log-pressure at the cloud top in bars.
Haze_thick: thickness of the slab cloud expressed as the difference in
log-pressure between the cloud top and cloud base.

Model 3:
Haze_Pabund: log-peak number density of haze particles in the cloud in cm^-3.
Haze_size: log-radius of uniform haze particles in microns.
Haze_meanP: log-pressure of the peak density of the cloud.
Haze_scale: standard deviation of the cloud density in dex of pressure.

Currently supported cloud types:
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

End:
Statistical parameters needed for an MCMC fit:
deltaL: wavelength offset in nm to correct for calibration errors of the 
observations. Default: 0.0
logf: log of error bar correction factor (multiplies the size of the error
bars of the input spectrum). Default: 1.0

NOTE: In APOLLO, the output will add one additional parameter to the end of this
list, which is the effective temperature computed for the model.

7. SPECTRUM FILE FORMAT

The spectrum file format has 6 columns and no header:
1. Wavenumber of short-wave end of wavelength bin in cm^-1
2. Wavenumber of long-wave end of wavelength bin in cm^-1
3. For input spectrum files, not used, usually identical to column 6. In
   Apollo Makespectrum output files, the model spectrum without noise added.
4. Lower error bar of flux
5. Upper error bar of flux
6. For emission spectra, observed flux in erg s^-1 cm^-3; for transit spectra,
   the observed planet-star radius ratio. In Apollo Makespectrum output files,
   only this column includes the computed noise.

The spectrum file must be ordered from highest wavenumber to lowest wavenumber
within each band because it is used to truncate the computed spectrum for
speed.

If the input file is give as one continuous band (with no gaps between bins),
the code will compute and fit the spectrum over the full width of the
wavelength range. If there are gaps between the bins, it will interpret the
spectrum as a set of wavelength bands or filters and fit only on those bands.

It is also possible to do a fit with continuum normalization with APOLLO. To do
this, add the codeword "Polyfit" to the "Data" line in the input file. Then,
list the bands to be fit in order from highest to lowest wavenumber, and then,
list the bands to normalizing on from highest to lowest wavenumber. The code
will interpret everything after a band increases in wavenumber as the set to be
normalized upon. APOLLO will average each of these bands and do a polynomial
fit of the input spectrum with the same degree as the number of normalization
bands. The bands to be fit will be divided by the polynomial fit to normalize
out any pseudocontinuum that would complicate the fit.

8. MOLECULAR CROSS SECTION TABLE FORMAT

The molecular cross section files created by MakeHaze consist of a grid of
630 x 21205 entries. The columns of the table correspond to a grid 18 pressure
points by 35 temperature points. The pressure points are logarithmically spaced
from 300 bar to 1 microbar, and the temperature points are logarithmically
spaced from 75 K to 3500 K, which may be extrapolated up to 4000 K.

Two sets of tables are currently packaged with APOLLO: one covering the range
of 0.6-5.0 microns and one covering the range of 5-28 microns.
The columns of the table correspond to wavelength, spaced logarithmically. As
packaged, the tables have a resolution of R = 10000*ln(10) = 4,343. However,
the code is currently set up to use only one set of tables at a time.

Note that APOLLO uses the boundary values of the cross section tables if the
temperature or pressure values go outside the table.

9. HAZE TABLE FORMAT

Header: 6 columns describing the table.
1. Number of frequency points (ideally the same as the cross section tables).
2. Minimum wavelength in microns.
3. Maximum wavelength in microns.
4. Number of particle sizes.
5. Minimum log-particle size in microns.
6. Maximum log-particle size in microns.

Each line of the table begins with the reference wavenlength in microns and is
followed by the the cross section per particle at each particle size.

NOTE: Atmosphere.cpp contains functions to compute haze cross sections directly
from the indices of refraction via Mie scattering, but these are not supported
in v0.9.13.1 because it is much slower. Use MakeHaze to compute new haze
opacity tables.

In addition to the pre-computed cross sections, the included MakeHaze.cpp
program can compute new cross section tables from a list of real and complex
indices of refraction of the haze material, which are included in the
Opacities directory. The index of refraction tables have three columns:
1. Wavelength in microns.
2. Real index of refraction.
3. Complex index of refraction.

10. OUTPUT FILE FORMATS

APOLLO prints a list of MCMC samples to a file in the specified output
directory. Be sure there is enough memory to store the output in that
directory, especially on a computing cluster. The output file is about 6 MB per
1000 steps, using the standard format of printing only the last 10% of samples.

The first line of the output file contains the number of walkers, number of
steps, and number of parameters. The second line is the parameter list as
specified in the input file, plus the extra derived parameter for effective
temperature. Subsequent lines list all of the paramters for each sample. By
default, only the last 10% of samples are printed, after the burn-in, to reduce
file sizes. The full sample array can be printed by adding the flag "Full" to
the "Output" line in the input file.

APOLLO also produces a new input file for the object with the best fit
parameters from the retrieval. This makes it much easier to feed it back into
Makespectrum to plot it. (This can't be done directly by APOLLO because of how
emcee handles things.) This file goes to the work directory and ends with
"retrieved.dat"

APOLLO also produces three plots, which go to the "plots" subdirectory:
1. A waterfall plot of the "basic" parameters such as radius and gravity, with
effective temperature added.
2. A waterfall plot of the molecular abundances.
3. A plot of the retrieved temperature-pressure profile with 1-sigma limits.

Makespectrum outputs a full-resolution spectrum file. Commented sections of the
code also create spectra in JWST's spectroscopic and photometric modes. For
photometric modes, the file name for the thruput function must be hard-coded.
These output files have the same format as the input spectrum and are sent to
the "modelspectra" subdirectory.

Makespectrum also plots the retrieved spectrum against the input spectrum with
residuals. This figure goes to the "plots" subdirectory.

Plot Apollo replicates all of the plots produced by APOLLO or Makespectrum,
depending on the type of input file.

*******************************************************************************

PRIORITY CHANGES

NOTE: The h-minus opacity PLUS the fencepost error in gamma PLUS the passing of
    the normalization to the Planet object could all be contributing to the
    code failing on L-dwarfs.

1. Write a script to compute the effective temperature with error bars.

2. Improve the observation binning for spectra with gaps.

3. Implement a default parametric T-P profile.

4. Extend the cross section tables to 4000 K (36 temps).

5. Need to find a way to set minP and maxP correctly in Plot Apollo.

6. Add an option to request specific spectroscopic and photometric modes in
    Makespectrum in the input file.

7. Make the code work for NIR and MIR simultaneously.

8. Look at the needs for other opacity tables: higher-resolution tables, wider
    frequency range, and/or binary format. Need to talk to people about the
    needs of the code to decide what to do about it.

CHANGE LOG

v0.9.13.2
Cleaned up maketar.sh to prevent it from including extra files.
Corrected the H2 scattering cross section in Makespectrum.
Fixed the alkali metal opacities, which were wrongly multiplied by 10,000.
Removed the effective temperature calculation because it was a shortcut that
didn't work. A script to compute it correctly from the luminosity is pending.
Verified that the fencepost error in gamma appears to have resolved the hot top
of the atmosphere problem.

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

TO DO LIST

*Interface*
Give each cloud model a name to create a switch for them, including 2-param?
    Almost definitely, but not a priority because I need to get the code up and
    running again to figure out where I need to go with the clouds.
Make a function to allow Makespectrum to run from the command line in the
    Python environment without rereading the cross section tables. Important,
    but long-term development working with the people doing the web interface.
     
*Code Mechanics*
Combine the Planet_layer and Planet_auto files in some respect, possibly
    moving the functions they share into a separate file. Viable, but would be
    a pain to do right. It might require overloading constructors. But if I
    add multiple parametric T-P profiles, I'll probably need to.
    Actually, take a closer look at how I handle the profiles now.
Add options to test for convergence on the fly as both the only and an
    additional stopping criterion to max_steps. I think emcee can do that, but
    I need to make sure it's finding the correct answers before I try it.
Make the code able to handle tables of different sizes in T-P space. The number
    of temperature and pressure points are hard-coded in the Planet
    constructors, but more importantly, it assumes the table for each species
    is the same size. This is mainly if I want to implement things like Exomol,
    but I've got enough problems figuring out what to do with the opacity
    tables already, so I can't say much about it until I figure those out.
Implement autocorrelation to check for convergence.

*Model Completeness*
Figure out what's up with the hot top of the atmosphere. Maybe add a constraint
    to force it to be close to the layer below? Or would that make it worse?
    Could be the fencepost error in gamma.
I'm not accounting for reflection, which requires single-scattering albedo.
    This might be part of Arthur Adams's program.

*Wish List*
Add additional cloud models. This will probably be the top priority from this
    sub-list.
Allow nitrogen as a filler along with hydrogen. This will be important if I
    ever get around to doing terrestrial planets.
Add additional molecules. H2-minus seems fairly important, but I'm not sure
    what else would be.
Add a chemical equilibrium code. Note: suppress reactions that are slower than
    the rotation and/or circulation timescale.
Make sure the 2-stream algorithm works with clouds. (Split the cloud
    absorption and scattering opacities, etc.) Note 1: need to set the optical
    depth to infinity if it is below an opaque cloud deck. Note 2: this is what
    Arthur Adams is working on.
Use blobs and alternate likelihood/probability functions to combine the
    functionality of APOLLO and Makespectrum. (Probably in Version 2 or
    something.)

*2-Part Models*
Implement a parameter for partial cloud coverage. Think about how many
    parameters I'll need to fit a cloudy and a clear model to a spectrum.
    But need to get single-component models working first.
Implement combined transit and eclipse retrieval. This may involve dayside and
    nightside atmosphere models, but the other part is that it will need to
    compute separate transit and emission spectra for each model.
Implement phase curve observations. Basically the same as partial cloud
    coverage: dayside and nightside models with the "partial coverage"
    parameter being the phase angle.

DISCARDED IDEAS

Make the settings not depend on being on separate lines. Dropped for the same
    reason Python uses whitespace. I've seen what happens when it's not
    enforced, and it's an unreadable mess.
Allow extrapolation of temperature outside the table. Using the boundary values
    is probably pretty safe in pressure space. In temperature, I've got it
    covered for 75 K to 3500 K, and that will cover the actual region being
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
