import numpy as np


def calculate_CtoO_and_metallicity(list_of_gases, gas_logabundances):
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
        if ccompounds[i] in list_of_gases:
            j = list_of_gases.index(ccompounds[i])
            carbon = carbon + cmult[i] * (
                10 ** gas_logabundances[j]
            )  # -1 because of hydrogen
    for i in np.arange(0, len(ocompounds)):
        if ocompounds[i] in list_of_gases:
            j = list_of_gases.index(ocompounds[i])
            oxygen = oxygen + omult[i] * (10 ** gas_logabundances[j])
    for i in np.arange(0, len(zcompounds)):
        if zcompounds[i] in list_of_gases:
            j = list_of_gases.index(zcompounds[i])
            metals = metals + zmult[i] * (10 ** gas_logabundances[j])

    ctoo = carbon / oxygen
    fetoh = np.log10(metals / 0.0196)

    return np.array([ctoo, fetoh])
