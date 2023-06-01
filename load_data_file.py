from ast import literal_eval as guess_type
from dataclasses import dataclass, field
from collections import namedtuple
from glob import glob
import numpy as np
import pandas as pd
from src.ApolloFunctions import ConvBands, BinBands
from warnings import warn


@dataclass
class SpectralData:
    filepath:     str
    column_names: list[str] = [
        'wavelength_bin_start',
        'wavelength_bin_end', 
        'model_flux',
        'error_lower',
        'error_upper',
        'model_flux_with_noise'
        ],
    band_labels:  list[str] = field(default_factory=list)

    def load_from_file(self) -> pd.DataFrame:
        return pd.read_csv(
            self.filepath,
            delimiter='\s+',
            index_col=False,
            names=self.column_names
            )

    def load_from_file_with_bands(self) -> pd.DataFrame:
        unbanded_data = self.load_from_file()

        band_boundary_indices = np.where(
            np.asarray(unbanded_data['wavelength_bin_start'])[1:] !=
            np.asarray(unbanded_data['wavelength_bin_end'])[:-1]
            )[0] + 1

        iterated_band_labels = np.concatenate(
            [[band_label]*len(split_data)
              for band_label, split_data
              in zip(
                  self.band_labels,
                  np.split(unbanded_data.index, band_boundary_indices)
                  )
              ]
            )

        banded_data = unbanded_data.set_index(iterated_band_labels)
        self.as_dataframe = banded_data

        band_boundaries = np.sort(
            [
                [np.min(banded_data.loc[band_label].wavelength_bin_start),
                 np.max(banded_data.loc[band_label].wavelength_bin_end)]
                for band_label in self.band_labels
                ],
            axis=0
            )
        self.band_boundaries

        return {
            "band_boundaries": band_boundaries,
            "banded_data": banded_data
            }


    def convolve_and_bin(
            self,
            convolution_factor: float,
            binning_factor:     float
            ) -> pd.DataFrame:
        if not hasattr(self, 'banded_data'):
            raise AttributeError('Please load the spectrum from file first!')

        convert_to_list = lambda entry: [
            np.asarray(self.banded_data.loc[band_label])
            for band_label in self.band_labels
            ]

        flux_in_binning_format = [
            np.asarray(self.banded_data.loc[band_label].model_flux)
            for band_label in self.band_labels
            ]

        uncertainties_in_binning_format = [
            np.asarray(self.banded_data.loc[band_label].error_upper)
            for band_label in self.band_labels
            ]

        convolved_data = ConvBands(convert_to_list(),
                                   uncertainties_in_binning_format,
                                   convolution_factor)

        binned_data = BinBands(self.banded_data)


def load_data_from_file(
        data_filepath:                     str,
        column_names: list[str] = [
            'wavelength_bin_start',
            'wavelength_bin_end', 
            'model_flux',
            'error_lower',
            'error_upper',
            'model_flux_with_noise'
            ],
        band_labels: list[str] = field(default_factory=list)
        ) -> pd.DataFrame:
    unbanded_data = pd.read_csv(data_filepath,
                                delimiter='\s+',
                                index_col=False,
                                names=column_names)
    
    band_boundaries = np.where(
        np.asarray(unbanded_data['wavelength_bin_start'])[1:] !=
        np.asarray(unbanded_data['wavelength_bin_end'])[:-1]
        )[0] + 1

    iterated_band_labels = np.concatenate(
        [[band_label]*len(split_data)
          for band_label, split_data
          in zip(band_labels, np.split(unbanded_data.index, band_boundaries))
          ]
        )

    banded_data = unbanded_data.set_index(iterated_band_labels)
    return banded_data


def get_band_boundaries_from_data(data_dataframe_index):
    band_labels = np.unique(data_dataframe_index)

    band_boundaries = np.sort(
        [
            [np.min(data_dataframe.loc[band_label].wavelength_bin_start),
             np.max(data_dataframe.loc[band_label].wavelength_bin_end)]
            for band_label in band_labels
            ],
        axis=0
        )
    return band_boundaries


x = load_data_from_file('data/mock-L_JHK.dat', band_labels = ['J', 'H', 'K'])
# print(x.loc["J"].wavelength_bin_end)
print(get_band_boundaries_from_data(x))