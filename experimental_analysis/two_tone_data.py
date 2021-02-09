"""
Two tone data measured made on the diploe qubit for the paper
"""
from collections import defaultdict
import numpy as np


class TwoToneData:
    TRANSITION_12_FILENAME_SUFFIXES = ["m3", "m2", "m1", "1", "2", "3"]
    TRANSITION_23_FILENAME_SUFFIXES = ["m3b", "m2b", "m1b", "1b", "2b"]

    def __init__(self):
        self.data = defaultdict()

    def __getitem__(self, field: str):
        return self.data[field]

    @staticmethod
    def convert_mA_to_flux(
        mA_array: np.ndarray, offset: float, period: float
    ) -> np.ndarray:
        """
        __ Parameters __
        mA_array: array of mA values to convert
        offsset:  offset from 0
        period:   period in mA

        __ Description __
        Converts array measured in experimental mA to flux units
        """

        return (mA_array - offset) / period

    def load_data(self, base_file_name: str, offset: float, period: float):
        """Base path to the files"""

        self.data = defaultdict(list)

        for suffix in self.TRANSITION_12_FILENAME_SUFFIXES:
            loaded_array = self.import_and_sort_array_from_file(
                f"{base_file_name}_{suffix}.txt", 0
            )
            loaded_array[0] = self.convert_mA_to_flux(loaded_array[0], offset, period)
            self.data["flux_12"].extend(loaded_array[0])
            self.data["spectrum_12"].extend(loaded_array[1])

        for suffix in self.TRANSITION_23_FILENAME_SUFFIXES:
            loaded_array = self.import_and_sort_array_from_file(
                f"{base_file_name}_{suffix}.txt", 0
            )
            loaded_array[0] = self.convert_mA_to_flux(loaded_array[0], offset, period)
            self.data["flux_23"].extend(loaded_array[0])
            self.data["spectrum_23"].extend(loaded_array[1])

    def import_and_sort_array_from_file(
        self, file_name: str, column_to_sort_by: int
    ) -> np.array:
        """
        __ Description __
        Import the 2D array and sort by the i-th column
        """

        # 1 - import array
        loaded_array = np.loadtxt(file_name).transpose()

        # 2 - sort array by the i-th column
        sorted_index = np.argsort(loaded_array[column_to_sort_by])
        loaded_array[0] = np.array(loaded_array[0])[sorted_index]
        loaded_array[1] = np.array(loaded_array[1])[sorted_index]

        return loaded_array
