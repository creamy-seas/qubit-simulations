import numpy as np


class GenericConverter:
    """Class that runs generic conversions"""

    def __init__(self, quantum_constants):
        self.constants = quantum_constants

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

    def convert_energy_to_GHz(self, energy: float) -> float:
        return energy / (self.constants.h * 10 ** (9))
