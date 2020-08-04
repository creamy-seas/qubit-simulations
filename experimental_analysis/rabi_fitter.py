from typing import Tuple, List

import numpy as np

class RabiFitter:
    def load_data(self, filename: str, colX: int, colY: int) -> Tuple[List[float], List[float]]:
        """
        __ Parameters __
        [str] filename: where to load rabi data from
        [int] colX, colY: which columns to treat as the X (time) and Y(amplitude)

        __ Description __
        loads rabi oscillation data
        """

        temp_load = np.loadtxt(filename)
        rabi_x = temp_load[:, colX]
        rabi_y = temp_load[:, colY] * 10 ** 6
        
        return (rabi_x, rabi_y)

    def rabi_fit_function(self, x, A, tDec, t_p, phi, D):
        """
        __ Description __
        Fits Rabi oscillations of the format
        A e^(-t/tDec) cos(2pi*t/tP+phi) + D
        """

        return A * np.sin(2 * np.pi * x / t_p + phi) * np.exp(-x / tDec) + Drabi_fitter