class ErrorEvaluator(object):
    """Class that evaluates error between simulation and experiment

    """

    def __init__(self):
        pass

        def experimental_data_error(self):
        """
        __ Description __
        Error between the experimental data points and the simulation
        """
        try:
            logging.info(
                "==> 'experimental_data_error' comparing simulation to experiment"
            )

            if not hasattr(self, "spectrum_simulation_12"):
                raise RuntimeError(
                    "run the simulation first before calling 'experimental_data_error'"
                )
            if not hasattr(self, "spectrum_experimental_12"):
                raise RuntimeError(
                    "import experimental data before calling 'experimental_data_error'"
                )

            # 1 - prepare parameters
            error_cumulative = 0

            # 2 - transition12
            for i in range(0, len(self.flux_list_experimental_12)):
                # a - find the simulation entry corresponding to the experimental data point
                entry = np.where(self.flux_list == self.flux_list_experimental_12[i])[
                    0
                ][0]

                # b - compute the difference
                error_cumulative = (
                    error_cumulative
                    + (
                        self.spectrum_simulation_12[entry]
                        - self.spectrum_experimental_12[i]
                    )[0]
                    ** 2
                )

            # 3 - transition23
            for i in range(0, len(self.flux_list_experimental_23)):
                # a - find the simulation entry corresponding to the experimental data point
                entry = np.where(self.flux_list == self.flux_list_experimental_23[i])[
                    0
                ][0]

                # b - compute the difference
                error_cumulative = (
                    error_cumulative
                    + (
                        self.spectrum_simulation_23[entry]
                        - self.spectrum_experimental_23[i]
                    )[0]
                    ** 2
                )

            logging.info("'experimental_data_error' finished")

            return error_cumulative
        except TypeError as e:
            print(e)
