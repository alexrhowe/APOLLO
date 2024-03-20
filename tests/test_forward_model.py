import numpy as np
import sys

sys.path.append("../newAPOLLO")

from make_forward_model_from_file import generate_model_spectrum

TEST_2M2236_FILEPATH = "/Users/arthur/Documents/Astronomy/2019/Retrieval/Code/Results/2M2236/spectral_range/2M2236.Piette.G395H.retrieved.ensemble.2024-01-22.dat"
TEST_MATCH_2M2236_SPECTRUM = "/Users/arthur/Documents/Astronomy/2019/Retrieval/Code/Results/2M2236/2M2236.Piette.G395H.cloud-free.2024-01-22.retrieved.Spectrum.binned.dat"


def absolute_difference_of_spectra(spectrum, reference_spectrum):
    return np.max(np.abs(spectrum - reference_spectrum))


def load_and_check_if_spectra_match(
    spectrum_filepath, reference_spectrum_filepath, absolute_tolerance=1e-15
):
    spectrum = generate_model_spectrum(spectrum_filepath)
    reference_spectrum = np.loadtxt(reference_spectrum_filepath)[:, 2]

    return (
        absolute_difference_of_spectra(spectrum, reference_spectrum)
        <= absolute_tolerance
    )


#############################################################
########################### TESTS ###########################
#############################################################
def test_2M2236_forward_model():
    assert load_and_check_if_spectra_match(
        TEST_2M2236_FILEPATH, TEST_MATCH_2M2236_SPECTRUM
    )
