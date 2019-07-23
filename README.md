# APOLLO
Atmosphere retrieval code for exoplanets

USER GUIDE TO THE APOLLO ATMOSPHERE RETRIEVAL CODE

APOLLO v0.9.10 Release Candidate
18 June 2019

Alex R. Howe

Needed files:\
Apollo.tar\
	README.txt\
	maketar.sh\
	APOLLO source code\
	MakeHaze source code\
	JWST filter list\
	JWST filter functions\
	Example files\
Apollo.opacs.tar.gz\
	Molecular cross section tables\
	Aerosol optical properties tables\
	Aerosol optical cross section tables

DESCRIPTION

APOLLO is an atmosphere code for extrasolar planets designed to fit\
forward models to both transit and emission spectra and retrieve atmospheric\
properties. APOLLO uses an MCMC algorithm to find the best fit to a spectrum\
for a given parameter set. The code is designed to be modular enough to have\
many options that may be turned on and off, including the type of\
observations, a flexible molecular composition, multiple cloud prescriptions,\
multiple temperature-pressure profile prescriptions, multiple priors, and\
continuum normalization. It uses cgs units unless otherwise specified.

APOLLO is written in C++ and Python using the Cython extension. Anaconda is\
the recommended build of Python given that scipy is required. The emcee\
package is also needed. The original build environment includes:

Python 2.7.13\
Anaconda 4.4.0 for Python 2\
emcee 2.2.1\
Cython 0.25.2\
GCC 4.9.4

The build environment for APOLLO on the Flux cluster at U of M requires only:

module load python-anaconda2

FUNCTIONALITY

APOLLO version 0.9.11 is verified to be accurate and stable with the one-stream\
radiative transfer algorithm for transit, secondary eclipse, and direct\
imaging; for layered and 5-parameter T-P profiles; for cloud-free models, an\
opaque cloud deck, a cloud slab of uniform particle size and density, and a\
cloud of gaussian particle density and uniform particle size.

The two-stream radiative transfer algorithm is verified only for cloud-free and\
opaque cloud deck eclipse spectra.

APOLLO is bundled with the Apollo Makespectrum program, which computes spectra\
of forward models complete with noise levels for a model JWST pipeline.\
Makespectrum uses the same input file format as APOLLO, although it uses a\
slightly different set of parameters. It generates spectrum files that can be\
used as "Model"-type observations in APOLLO. Makespectrum currently produces\
only JWST spectroscopic modes. However, it does output the full-resolution\
spectrum to compute noise levels for photometric modes or model pipelines for\
other observatories.

APOLLO is also bundled with the MakeHaze program, which computes new haze\
opacity tables using Mie scattering from a table of indices of refraction for\
the haze material.

INSTALLATION

To install APOLLO, unpack the Apollo.tar tarfile, and compile the code using\
the Cython setup file:

python setup.py build_ext --inplace

This script contains all of the setup needed to run both APOLLO and\
Makespectrum. APOLLO is also bundled with tables of molecular and haze\
particle cross sections. You will need to unzip and unpack this bundle to run\
the code. The zipped version is 620 MB, so you will need to acquire it via USB\
drive or FTP.

As packaged, these tables are to be read from a directory named Opacities\
in the parent directory to the code directory. This allows multiple builds of\
the code to be run in different directories. If you use a different directory,\
you will need to change the file path in the constructors of Atmosphere.cpp and\
Atmosphere.haze.cpp, and in readopac in Planet_layer.cpp and Planet_auto.cpp.

Use the included Makefile to build the MakeHaze program.

RUNNING

APOLLO is called with a single command line argument, which is the name of\
the input file. If no input file is specified, it will default to the bundled\
example file, example.resolved.dat. If a second command line argument,\
"Serial", is added, it will override parallel operation and initiate a\
minimal-length serial run for testing purposes. APOLLO may also be run in\
serial mode if parallel mode is turned off in the input file.

In serial model, APOLLO may be run as a regular Python program:

python Apollo.py <input file>\
OR\
python Apollo.py <input file> Serial

The default mode of APOLLO is parallel processing. On Flux, this can be done\
using a PBS script incorporating the following command:

mpirun python Apollo.py <input file> >> /scratch/lsa_flux/<uniqname>/Apollo.log

Sending the log output to Scratch is recommended to reduce memory usage on\
Flux. If it is not included, the log file will be sent to the local directory.

The example files include the parameters for GJ 570D for resolved spectra, and\
for HD 189733b for transit and eclipse spectra.

INPUT FILE FORMAT

The input file consists of two sections: a list of settings for a particular\
configuration of the code, and a list of parameters for the forward model to\
be fit. The settings may be listed in any order, but they all must come before\
the parameters. The parameters must also be grouped into blocks with specific\
labels.

Defaults are set so that any parameters may be omitted from the input file.\
The header line for the parameter section must be included, but the code will\
execute with only that line present.

The settings block of the input file may contain any of the following settings.

Mode, 1 field:\
1. Type of spectrum. "Resolved" for resolved objects in emission. "Transit"\
for transit spectra. "Eclipse" for secondary eclipse. Default: "Resolved"

Parallel, 1 field:\
1. Mode of operation. "Parallel" for parallel operation using emcee's\
built-in pool function. "Serial" for serial operation from the command line.\
The "Serial" command line argument overrides this setting.

Object, 1 field:\
1. The name of the target object. This field may be anything and is used only\
to set the name of the output file. Default: "example"

Star, 3 fields:\
1. Stellar Temperature in kelvin. Default: 5770 (Solar)\
2. Stellar radius in Solar radii. Default: 1.0\
3. Semi-major axis of the planet in AU. Default: 1.0\
These parameters are not used in the "Resolved" mode.

Location, 3 fields:\
1. Distance in pc. Default: 10.0\
2. Right ascension in decimal degrees. Default: 0.0\
3. Declination in decimal degrees. Default: 0.0

Data, 2 fields:\
1. Spectrum file name. This file contains the observations to be fit. Default:\
"example.obs.dat" (the bundled example file)\
2. Spectrum file type. "Observations" for a real observed spectrum. "Model"\
for a forward model from the code. This is used only to set the output file\
name. Default: "Observations"

Degrade, 1 field:\
1. Degradation factor. This is an integer factor by which to degrade the\
spectral resolution of the forward model. Degrading the spectrum speeds up\
the code proportionally. If the degraded forward model is still a factor of\
5-10 times higher-resolution than the input spectrum, this will still\
produce reliable results. A degradation factor of 3 is recommended unless the\
input spectrum is especially high resolution (e.g. JWST observations).\
Default: 3. Hard-coded to 1 in Makespectrum.

N_Walkers, 1 field:\
1. Number of walkers to be initialized by emcee. Default: 8*N. The "Serial"\
command line argument overrides this setting and sets nwalkers to the\
minimum allowed value of 2*N+2.

N_Steps, 1 field:\
1. The number of steps for the MCMC algorithm to take. This must be enough\
for the walkers to converge and may require adjustment for any given object.\
Default: 30000. The "Serial" command line argument overrides this setting and\
sets nsteps to 10.

Pressure, 2 fields:\
1. Minimum log-pressure in bars over which to compute the radiative transfer.\
Default: -3.0\
2. Maximum log-pressure in bars over which to compute the radiative transfer.\
Default: 2.5

Vres, 1 field:\
1. Number of layers to use for the radiative transfer calculation. Default: 71

Streams, 1 field:\
1. Which radiative transfer function to use. 1 for 1-stream, 2 for 2-stream.\
Default: 1\
Note that the 2-stream algorithm is much slower than the 1-stream algorithm,\
by as much as 2 orders of magnitude. It is currently implemented only for\
cloudless atmospheres or opaque cloud decks in emission.

Prior, 1 field:\
1. Type of prior to be used. Currently supported priors are "Uniform" and\
"Normal". Default: "Uniform"

Gray, 1 field:\
1. Boolean switch to use a gray atmosphere model. Default: False.

Output, 1 field:\
1. Output directory. Default: "/nfs/turbo/lsa-arhowe". For operations on Flux,\
this should be set to "/nfs/turbo/lsa-<uniqname>/" to ensure there is\
sufficient long-term storage for the output file. Default "modelspectra/" in\
Makespectrum.

The parameters section of the input file must begin with a header beginning\
with the string "Parameter". This header shows the format of the parameter\
section:

"Parameter   Initial   Mu    Sigma    Min     Max"

The model parameters are listed in blocks with headers to separate them in the\
code's handling. There are 5 blocks currently implemented, which may be in any\
order. Any individual block may be omitted, but the header must be included to\
function. Each parameter is listed with 6 fields:\
1. Parameter name.\
2. Initial guess.\
3. Mean of the prior. (Only used for normal priors.)\
4. Standard deviation of the prior. (Only used for normal priors.)\
5. Minimum cutoff of the prior.\
6. Maximum cutoff of the prior.

The implemented values of each block are listed below.

Basic:\
Bulk parameters of the planet.\
Rad: Radius in Earth radii\
RtoD: Log of radius-to-distance ratio\
RtoD2U: Log of squared radius-to-distance ratio\
Default radius: 1 Jupiter radius in linear space\
Log(g): log(g) in cm/s^2. Default: 4.1\
Cloud_Base: log-pressure in bars at the cloud deck that serves as the base of\
the atmosphere for the forward model. Default: 2.5\
P_cl: alternate name for Cloud_Base.

If the radius is not specified, the default is set to 11.2 Earth radii.\
However, the code will automatically normalize the spectrum to the total flux\
of the observations.\
Also note that APOLLO includes an internal prior that the mass must be less\
than 80 Jupiter masses. It is important to set the initial radius and gravity\
so that the mass falls below this limit.

Gases:\
List of molecules included in the model in log-abundance. The h2\
(actually h2 and he combined) is set by the difference 1 minus the other\
abundances.\
Currently supported molecules:\
h2\
h2o\
ch4\
co\
co2\
nh3\
h2s\
Burrows_alk (Burrows Alkali Metals)\
Lupu_alk (Lupu Alkali Metals)\
crh\
feh\
tio\
vo

Atm, 2 fields:\
1. Type of T-P profile for the atmosphere. "Layers" for an arbitrary number of\
temperature points spaced evenly in log-pressure space. "Parametric" for a\
5-parameter analytic T-P profile.\
2. Optional. Type "Verbatim" to use a layered profile directly from the input\
file. If not included, the T-P profile will be interpolated to a number of\
points set by Vres (default: 71). The default if the T-P profile is omitted is\
an isothermal atmosphere at 1500 K.

T#: an arbitrary number of temperature points in kelvin. 15 is recommended.\
Use 1 point for an isothermal atmosphere. The default cross section tables have\
boundaries of 75 K and 4,000 K, which are the default limits.

gamma: a smoothing parameter that penalizes wide swings in the T-P profile in\
the likelihood function. If included, the program will automatically implement\
a smoothed T-P profile. Defined using an Inverse Gamma Function as in\
Line et al. (2015).

Clouds, 2 fields:\
1. Cloud model number.\
2. Haze type (optional).

Currently supported cloud models:\
0: No clouds (0 parameters).\
1: Opaque cloud deck (1 parameter).\
2: Slab cloud profile with uniform particle size (4 parameters).\
3: Gaussian cloud profile with uniform particle size (4 parameters).

Model 1\
P_cl: log-pressure in bars at the cloud deck that serves as the base of\
the atmosphere for the forward model.

Model 2\
Haze_abund: log-number density of haze particles in the cloud in cm^-3.\
Haze_size: log-radius of uniform haze particles in microns.\
Haze_minP: log-pressure at the cloud top in bars.\
Haze_thick: thickness of the slab cloud expressed as the difference in\
log-pressure between the cloud top and cloud base.

Model 3\
Haze_Pabund: log-peak number density of haze particles in the cloud in cm^-3.
Haze_size: log-radius of uniform haze particles in microns.\
Haze_meanP: log-pressure of the peak density of the cloud.\
Haze_scale: scale height of the cloud density in dex of pressure.

Currently supported haze types:\
H2SO4\
Polyacetylene\
Tholin\
Corundum\
Enstatite\
Forsterite\
Iron\
KCl\
Na2S\
NH3Ice\
Soot\
H2OCirrus\
H2OIce\
ZnS

End:\
Statistical parameters.\
deltaL: wavelength offset in nm to correct for calibration errors of\
observations. Default: 0.0\
logf: log of error bar correction factor (multiplies the size of the error\
bars of the input spectrum). Default: 1.0

In APOLLO, the output will add one additional parameter to the end of this\
list, which is the effective temperature computed for the model.

SPECTRUM FILE FORMAT

The spectrum file format has 6 columns and no header:\
1. Wavenumber of short-wave end of wavelength bin in cm^-1\
2. Wavenumber of long-wave end of wavelength bin in cm^-1\
3. Model flux in erg s^-1 cm^-3 (Identical to Column 6 for observations.)\
4. Lower error bar of flux\
5. Upper error bar of flux\
6. Observed flux in erg s^-1 cm^-3\
For Observations, Columns 3 and 6 are identical. For Models, Column 6 is\
Column 3 with noise added.\
The file must be ordered from highest wavenumber to lowest wavenumber because\
it is used to truncate the computed spectrum for speed.

In the input file is give as one continuous band (no gaps between bins), the\
code will compute and fit the spectrum over the full width of the wavelength\
range. If there are gaps between the bins, it will interpret the spectrum as a\
set of wavelength bands or filters and fit only on those bands.

It is also possible to do a fit with continuum normalization with APOLLO. To do\
this, list the bands to be fit in order from high to low wavenumber, and then\
list the bands to normalizing on from high to low wavenumber. The code will\
interpret everything after a band increases in wavenumber as the set to be\
normalized upon. APOLLO will average each of these bands and do a polynomial\
fit of the input spectrum with the same degree as the number of bands. The\
bands to be fit will be divided by the polynomial fit to normalize out any\
pseudocontinuum that would complicate the fit.

MOLECULAR CROSS SECTION TABLE FORMAT

The molecular cross section files consist of a grid of 630 x 21205 ASCII\
entries. The rows of the table correspond to a grid 18 pressure points by\
35 temperature points. The pressure points are logarithmically spaced from\
300 bar to 1 microbar. The temperature points are logarithmically spaced from\
75 K to 3500 K, which may be extrapolated up to 4000 K.

The columns of the table correspond to wavelength. They are spaced\
logarithmically in wavelength space from 0.6 microns to 5.0 microns.

Note that APOLLO uses the boundary values of the cross section tables if the\
temperature or pressure values go outside the table.

[Other table sizes: coming soon.]

HAZE TABLE FORMATS

Header: 6 columns describing the table.\
1. Number of frequency points (ideally the same as the cross section tables).\
2. Minimum wavelength in microns.\
3. Maximum wavelength in microns.\
4. Number of particle sizes.\
5. Minimum particle size in microns.\
6. Maximum particle size in microns.

Each line of the table begins with the reference wavenlength in microns and is\
followed by the [I think it's supposed to be the cross section per particle\
at each particle size, but it doesn't look right.]

[The original Atmosphere.cpp, now Apollo/Atmosphere.backup.cpp, computed the\
haze opacity by qe * pi*d^2/4 * n_haze, where qe is the extinction efficiency.]

Note that Atmosphere.cpp contains functions to compute haze cross sections\
directly from the indices of refraction via Mie scattering, but these are\
not supported in v0.9.11. Use MakeHaze to compute haze opacity tables.

In addition to the pre-computed cross sections, the included MakeHaze.cpp\
program can compute new cross section tables from a list of real and complex\
indices of refraction of the haze material, which are included in the\
Opacities directory. The index of refraction tables have three columns\
representing wavelength in microns, real index of refraction, and complex\
index of refraction. MakeHaze is compiled using the included Makefile. There\
is also a "blank.r.dat" table included as a dummy haze material\
with no refraction.

Note that MakeHaze.cpp should be able to compute >100 wavelengths x 161\
particle sizes per minute. If it is taking significantly longer, there may\
be a problem with the indices file.

OUTPUT FORMAT

APOLLO prints a list of MCMC samples to a file in /nfs/turbo/lsa-<uniqname>/\
This is necessary on Flux due to the size of the output file. If you want to\
run the code locally, you will need to change the path of "foutname" in the\
code.

The first line of the output file contains the number of walkers,\
number of steps, and number of parameters. Subsequent lines list all of the\
paramters of each sample in the order they appear in the input file, plus the\
extra parameter for effective temperature.

Makespectrum outputs full-resolution and JWST spectroscopy mode spectra to the\
modelspectra directory. These files have the same format as the input spectra.

*******************************************************************************

CHANGE LOG

v0.9.10\
Added a gaussian function to constants.cpp.\
Changed the number of cloud parameters to a numerical cloud model in the input\
file format.\
Implemented a gaussian cloud profile model.\
Fixed a bug in the haze opacity calculation in getTauProf.\
Completely rebuilt the polynomial normalization based on experiments on the\
Flux cluster.\
Fixed various mistakes along the way.\
Cleaned up the outputs throughout.\
Made the section header comments more informative.\
Made some of the error messages more informative.\
Rearranged a few lines to better harmonize APOLLO and Makespectrum.\
Changed the output of APOLLO to only output the last 10% of the length of the
chain to reduce file sizes. (The original is left, but commented out.)

v0.9.9\
Limited ability of haze calculation to go outside the size range due to the\
granularity of the table.\
Fixed a bug in the haze calculation that caused a crash in the cloud layer was\
narrower than one layer in the T-P profile.\
Changed the haze in the example files from tholin to enstatite.\
Fixed some format problems with the indicies of refraction files.\
Doubled the width and resolution in particles size of the haze opacity tables.\
Made changes to Atmosphere.cpp to prevent crashing due to rounding errors.\
Added more decimal places to haze opacity tables to prevent rounding errors.\
Added automatic flux normalization if no radius variable is included in the\
input file\
Added an option for polynomial continuum normalization.\
Updated the README file with the current input format.\
Added another tweak to Atmosphere.cpp to prevent crashing due to rounding\
errors in the interpolation.

v0.9.8\
Fixed handling of the gray atmosphere in Makespectrum.\
Repaired the functionality of a parametric T-P profile in transit.\
Removed a test output filename that I kept forgetting to remove.\
Fixed a bug that caused an error if the pressure limits were too wide.\
Updated example files without the Hazetype line.\
Added an example file for transits with a parametric T-P profile.\
Fixed a bug in getH that broke Planet_auto.cpp.\
Fixed the noise models for eclipse and transit.\
Harmonized minor differences between Planet_layer.cpp and Planet_auto.cpp.\
Fixed handling of gray boolean input parameter.\
Adjusted T-P profile to correct for variation in gravity with altitude.\
Fixed the Teff calculation with an opaque cloud deck included.\
Change the transit and eclipse example files to published data for HD 189733b.\
Added an eclipse, parametric profile example file.\
Fixed a sign error in the transit spectrum calculation.\
Rechecked all of the example files by plotting the results.

v0.9.7\
Changed the code to count the number of parameters instead of reading it from\
an input file entry.\
Fixed a bug where bounds were not set for teff0.\
Made the number of layers of the radiative transfer calculation variable.\
Allowed the layered T-P profile to be as short as 1 entry (isothermal).\
Added a check to make sure the parametric T-P profile is the right length.\
Added a default for a T-P profile being omitted.\
Made the code check for the smoothing parameter automatically.\
Updated example files with the new format.\
Added checks to handle missing fields in the input files.\
Allowed the parameter blocks to be in any order.\
Fixed some indentation errors.\
Fixed some bugs in the radius handling.\
Added an option to compute a gray atmosphere.\
Implemented the other available haze materials.\
Simplified the molecule and haze name handling.\
Removed the Hazetype line from the input file format.

v0.9.6\
Cleaned up the example files a bit.\
Made one-stream radiative transfer the default.\
Added a try-except to make sure the input file actually exists.\
Made a small change to the output filename format.\
Fixed a bug in setting the bounds for the cloud deck.\
Added a crude method in getTauProf in Planet_layer.cpp to compute T_eff.\
Added a check in case the probability function returns NaN.\
Fixed a bug that didn't send haze parameters to setParams.\
Moved the Atmosphere initialization into the Planet constructor to speed up\
execution.\
Implemented computing the effective temperature of the atmosphere for emission\
spectra and outputting it with the retrieved parameters.\
Checked all the example input files for errors.\
Cleaned up comments and extra print statements in the code.\
Fixed the cloudless atmosphere base in APOLLO.

v0.9.5\
Harmonized Apollo.pbs with the other example files.\
Changed the binning from cubic interpolation to linear because cubic crashes\
on flux.\
Automatically cut the computed spectrum to the length of the input spectrum to\
speed up execution.\
Improved the 1-stream algorithm to split a layer at a cloud deck boundary.\
Fixed the opaque cloud deck in emission and some stuff related to it.\
Added the needed header line to the haze cross section tables.\
Modified MakeHaze to add the header line to new tables.\
Fixed the 4-parameter haze model in emission and some stuff related to it.\
Fixed both clouds and haze in transit.\
Added cloud and haze example files to the tar ball.\
Updated parts of the README file that I'd been neglecting.

v0.9.4.1\
Turned the low optical depth correction in getFlux back on.\
Changed the getSpectrum function to include a switch between the 1-stream and\
2-stream radiative transfer algorithms.\
Cleaned up the mess the new getFlux was in.\
Fixed the getFluxOneStream function.\
Fixed a misnumbering of haze parameters in Planet_auto.\
Fixed name of APOLLO output file.\
Verified the functionality of the parametric T-P profile with the two-stream\
algorithm.\
Verified the two-stream algorithm for an isothermal atmosphere.\
Corrected an error in the Exponential Integral solver in constants.cpp.\
Fixed a duplicate use of mu in the one-stream radiative transfer.\
Verified both radiative transfer functions against the analytic solution for\
a gray atmosphere, and double-checked the normalization.

v0.9.4\
Fixed a bug in the transit spectrum calculation.\
Verified the functionality of transits.\
Changed the radiative transfer algorithm from one-stream to two-stream, using\
the short form of Mark Marley's algorithm used by Ben Burningham. This\
algorithm uses a tridiagonal matrix solver.\
Reversed order of profiles as needed for this algorithm.\
Fixed bugs that came up in getH, getT, and getP that resulted from the\
reversal of the profiles.\
Added taulayer array that is needed because the new algorithm uses layer\
optical depth instead of column optical depth.\
Added the "verbatim" switch to use the T-P profile directly from the input\
file rather than interpolating.

v0.9.3\
Added the MakeHaze program and Makefile to make haze cross section tables from\
indices of refraction.\
Added index of refraction tables to the Opacities directory.\
Fixed bug in MakeHaze that caused it to not account for particle size and\
replaced the applicale haze files.\
Fixed the Gaussian quadrature integration over solid angle to achieve more\
accurate results.\
Wrote a bash script to create the tarball of all needed files.\
Verified the functionality of an opaque cloud deck.\
Added a check to ensure the cloud deck is below the top of the model\
atmosphere.\
Removed the inconsisted usage of NIRISS_Throughput.dat from Filter.py and\
removed NIRISS_Throughput.dat from the prerequisites list.\
Harmonized the import lists between APOLLO and Apollo Makespectrum.\
Harmonized the read in files sections between APOLLO and Apollo Makespectrum.\
Fixed bugs in the handling of defaults.\
Verified functionality of default settings except for the number of parameters\
and the atmosphere profile.\
Verified functionality of the Howe & Burrows (2012) haze model.\
Changed the naming system to match the number of parameters in the input file.\
Fixed a bug that set the output file name incorrectly in APOLLO.

v0.9.2\
Set defaults for End parameters.\
Added an override switch for serial-run tests.\
Verified stability of the Rad and RtoD radius parameters.\
Verified functionality of the 5-parameter T-P profile (Line et al., 2012).\
* Note that Tint~300 and log(kappa)~0.3 provide a good fit to the example\
spectrum.\
Fixed a bug in GetScaOpac where some of the molecules were mislabelled.\
Implemented separate Eclipse and Resolved object modes accounting for\
stellar irradiation of close-in planets.\
Added a check to correct an odd number of walkers.

v0.9.1\
Added a field to the input file format for the output directory.\
Added a switch to the input file format for serial vs. parallel operation.\
Changed switches from On/Off to True/False.\
Added a field to the input file format for number of walkers.\
Added fields to the input file format for the minimum and maximum pressure over\
which to calculate the radiative transfer.\
Changed the order of the fields in the input file format for intuitiveness.\
Changed labelling of alkali metal files to reflect their sources.\
Made the reading of the input file more robust in the order of the\
non-parameter settings.\
Added default values for all of the input settings and the Basic parameters.\
Added a bundled example input file and observations file, currently GJ 570D,\
also allowing the code to be called with no arguments.\
Fixed a bug in the mean molecular weight calculation that omitted the\
high-temperature molecules.\
Fixed a bug that could cause a crash if extraneous lines are inserted in the\
parameter list.\
Documented dependencies.\
Bundled the spectrum creator.

TO DO

*Input File Structure*\
Give each cloud model a name to create a switch for them, including 2-param?\
Make the settings not depend on being on separate lines?\
Get rid of the "Model-Observations" distinction?

*Code Mechanics*\
Figure out what's going on with the haze cross section.\
Turn the APOLLO and Makespectrum codes into functions in a single script.\
Make a function to allow Makespectrum to run from the command line in the\
     Python environment without rereading the cross section tables.\
Make the code work past 5 microns.\
Add the functionality to compute haze scattering directly from the indices\
    of refraction.\
Combine the Planet_layer and Planet_auto files in some respect, possibly\
	moving the functions they share into a separate file.\
Allow extrapolation of temperature outside the table?\
Create source directory?\
Set up Atmosphere.cpp to read in index of refraction tables?\
Implement a way to test the required Num_Steps for different models?\
As written, the Teff calculation may crash if the cloud deck is at the top of\
the atmosphere.\
Get the noise model working for photometric filters.

*Model Completeness*\
Implement importing a stellar model.\
Extend cross section and indicies of refraction tables past 5 microns.\
Figure out what's up with the hot top of atmosphere. Maybe add a constraint to\
       force it to be close to the layer below? Or would that make it worse?\
I'm not accounting for reflection, which requires single-scattering albedo.\
Add additional Rayleigh coefficients.
       
*Wish List*\
Add additional molecules.\
Add additional cloud models.\
Allow nitrogen as a filler along with hydrogen.\
Implement different sizes and/or shapes of cross section tables.\
Implement a parameter for partial cloud coverage.\
Implement combined transit and eclipse retrieval.\
	  (Different day and night if needed.)\
Implement phase curve observations: 2-component model with 2 T-P profiles and\
	  fully-mixed abundances.\
Add a chemical equilibrium code.\
Suppress reactions in the equilibrium code that are slower than the rotation\
	 and/or circulation timescale.\
Make sure the 2-stream algorithm works with clouds.\
Split the cloud absorption and scattering opacities to make\
     the 2-stream radiative transfer algorithm work right with clouds.\
The 2-stream algorithm needs something to set the optical depth very high if\
    it is below a cloud deck.\
Make my output file processors clean enough to include with the package.
