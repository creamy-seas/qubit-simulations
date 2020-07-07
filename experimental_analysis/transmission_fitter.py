import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple


class TransmissionFitter:
    def __init__(self):
        self.transmission_x, self.transmission_y = (None, None)

    def load_data(
        self,
        data_file: str,
        x_column_index: int,
        y_column_index: int,
        convert_to_ghz: bool,
    ):
        """
        __ Parameters __
        x_column_index, y_column_index: column to treat as the X,Y coordinate
        convert_to_ghz: nromalise by 10^6 or not

        __ Description __
        Loads in the file and same
        """
        temp_load = np.loadtxt(data_file).transpose()
        self.transmission_x = temp_load[x_column_index, :]
        self.transmission_y = temp_load[y_column_index, :]

        if convert_to_ghz:
            self.transmission_x = self.transmission_x / 10 ** 9

    def filter(self, x_min: float, x_max: float, y_min: float, y_max: float):
        if self.transmission_x is None:
            raise RuntimeError("Load data first!")

        remove_index = []

        for i in range(0, len(self.transmission_x)):
            if (self.transmission_x[i] < x_min) or (self.transmission_x[i] > x_max):
                remove_index.append(i)

        for i in range(0, len(self.transmission_y)):
            if (self.transmission_y[i] < y_min) or (self.transmission_y[i] > y_max):
                remove_index.append(i)

        temp_x = []
        temp_y = []
        for i in range(0, len(self.transmission_y)):
            if i not in remove_index:
                temp_x.append(self.transmission_x[i])
                temp_y.append(self.transmission_y[i])

        self.transmission_y = temp_y
        self.transmission_x = temp_x

    def fit_function(self, x, Gamma1, Gamma2, Omega, offset) -> float:
        """
        __ Description __
        Fits the REAL part of general_data intesnity (which near the resonance will
        dominate over the IMAGINARY part (see major_project p.12))

        R[t]**2
        """

        x = x - offset
        return (
            1
            - Gamma1
            / (2 * Gamma2 * (1 + (x / Gamma2) ** 2 + Omega ** 2 / Gamma1 / Gamma2))
        ) ** 2

    def fit(self) -> Tuple[float, float, float, float]:
        """Perfrms fitting, returning parameters
        (Gamma1, Gamma2, Omega, Offset)
        """

        lowest_transmission = np.argmin(self.transmission_y)
        offset_at_lowest_transmission = self.transmission_x[lowest_transmission]

        (popt, _) = curve_fit(
            self.fit_function,
            self.transmission_x,
            self.transmission_y,
            bounds=(
                [0, 0, 0, offset_at_lowest_transmission - 0.1],
                [np.inf, np.inf, np.inf, offset_at_lowest_transmission + 0.1],
            ),
        )

        print(f"Gamma1   {popt[0]}")
        print(f"Gamma2:  {popt[1]}")
        print(f"Omega:   {popt[2]}")
        print(f"Offset:  {popt[3]}")

        return (popt[0], popt[1], popt[2], popt[3])
