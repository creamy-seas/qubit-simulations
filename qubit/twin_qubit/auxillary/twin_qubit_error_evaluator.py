import logging

import numpy as np
import pinject

# TODO: Work through this class
class ErrorEvaluator(object):
    """Class that evaluates error between simulation and experiment"""

    @pinject.copy_args_to_public_fields
    def __init__(self, twin_qubit_simulator, flux_list):
        pass

    def experimental_data_error(self):

        if not hasattr(self, "spectrum_simulation_12"):
            raise RuntimeError(
                "run the simulation first before calling 'experimental_data_error'"
            )
        if not hasattr(self, "spectrum_experimental_12"):
            raise RuntimeError(
                "import experimental data before calling 'experimental_data_error'"
            )

        error_cumulative = 0

        for i in range(0, len(self.flux_list_experimental_12)):
            data_index = np.where(self.flux_list == self.flux_list_experimental_12[i])[
                0
            ][0]

            # b - compute the difference
            error_cumulative = (
                error_cumulative
                + (
                    self.spectrum_simulation_12[data_index]
                    - self.spectrum_experimental_12[i]
                )[0]
                ** 2
            )

        # 3 - transition23
        for i in range(0, len(self.flux_list_experimental_23)):
            # a - find the simulation data_index corresponding to the experimental data point
            data_index = np.where(self.flux_list == self.flux_list_experimental_23[i])[
                0
            ][0]

            # b - compute the difference
            error_cumulative = (
                error_cumulative
                + (
                    self.spectrum_simulation_23[data_index]
                    - self.spectrum_experimental_23[i]
                )[0]
                ** 2
            )

        logging.debug("'experimental_data_error' finished")

        return error_cumulative
