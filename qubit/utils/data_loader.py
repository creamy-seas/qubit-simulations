class DataLoader(object):
    """Load data in from file"""

    def __init__(self):
        super(DataLoader, self).__init__()

    def experimental_data_load(self, mpl_axes, set_flux_list):
        """
        __ Parameters __
        mpl_axes: axes to plot the raw data on
        set_flux_list: True to set experimental data for the flux list

        __ Description __
        1 - load up experimental data
        2 - perform conversion from mA to flux
        3 - plot data
        4 - store the magnetic field points of the experiment

        Plots and finds differences between simulation and experiment
        """
        logging.debug("'experimental_data_load' Importing data files")

        # 1 - data files to load
        base_data_name = "data/Qubit15_5JJ_Q2_"
        transition_12 = ["m3", "m2", "m1", "1", "2", "3"]
        transition_23 = ["m3b", "m2b", "m1b", "1b", "2b"]

        # 2 - for each data file
        self.flux_list_experimental_12 = []
        self.spectrum_experimental_12 = []
        self.flux_list_experimental_23 = []
        self.spectrum_experimental_23 = []

        for data_set in range(0, len(transition_12)):
            # a - generate file name
            temp_name = base_data_name + transition_12[data_set] + ".txt"

            # b - import, sort and convert to flux
            temp_data = self.import_and_sort_array(temp_name, 0)
            temp_data[0] = self.convert_mA_to_flux(temp_data[0], 0.125, 0.7)

            # c - plot data
            if self.plot_or_not:
                mpl_axes.plot(
                    temp_data[0],
                    temp_data[1],
                    marker="o",
                    color="#004BA8",
                    markeredgecolor="C2",
                    markersize=7,
                    alpha=0.95,
                    linestyle="",
                )

            # d - store imported data in 1 array
            self.flux_list_experimental_12.extend(temp_data[0])
            self.spectrum_experimental_12.extend(temp_data[1])

        for data_set in range(0, len(transition_23)):
            # a - generate file name
            temp_name = base_data_name + transition_23[data_set] + ".txt"

            # b - import, sort and convert to flux
            temp_data = self.import_and_sort_array(temp_name, 0)
            temp_data[0] = self.convert_mA_to_flux(temp_data[0], 0.125, 0.7)

            # c - plot data
            if self.plot_or_not:
                mpl_axes.plot(
                    temp_data[0],
                    temp_data[1],
                    marker="o",
                    color="C4",
                    markeredgecolor="#fb2c07",
                    markeredgewidth="0.4",
                    markersize=5,
                    alpha=0.95,
                    linestyle="",
                )

            # d - store imported fluxes
            self.flux_list_experimental_23.extend(temp_data[0])
            self.spectrum_experimental_23.extend(temp_data[1])

        # 3 - tidy up by sorting arrays
        temp_array = self.sort_array(
            np.array([self.flux_list_experimental_12, self.spectrum_experimental_12]), 0
        )
        self.flux_list_experimental_12 = temp_array[0]
        self.spectrum_experimental_12 = temp_array[1]

        temp_array = self.sort_array(
            np.array([self.flux_list_experimental_23, self.spectrum_experimental_23]), 0
        )
        self.flux_list_experimental_23 = temp_array[0]
        self.spectrum_experimental_23 = temp_array[1]

        logging.debug(
            "Imported %i flux points",
            len(
                list(self.flux_list_experimental_12)
                + list(self.flux_list_experimental_23)
            ),
        )

        # 4 - set the flux array is required (then simulations are only done to compare with
        # experimental points)
        if set_flux_list:
            logging.debug("Set experimental flux points for simulation")
            self.flux_list = np.array(
                list(self.flux_list_experimental_12)
                + list(self.flux_list_experimental_23)
            )
            self.flux_list.sort()

        if self.plot_or_not:
            plt.show()

        logging.debug("'experimental_data_load' finished")


            def import_and_sort_array(self, file_name, i):
        """
        __ Parameters __
        file_name: file to import, with path and extensions
        i: column to sort by

        __ Description __
        Import the 2D array and sort by the i-th column

        __ Return __
        return the sorted array
        """

        # 1 - import array
        array_to_return = np.loadtxt(file_name).transpose()

        # 2 - sort array by the i-th column
        sorted_index = np.argsort(array_to_return[i])
        array_to_return[0] = np.array(array_to_return[0])[sorted_index]
        array_to_return[1] = np.array(array_to_return[1])[sorted_index]

        return array_to_return

    def sort_array(self, array_to_sort, column_to_sort_by):
        """
        __ Parameters __
        array_to_sort: array to perform sorting for
        column_to_sort_by: which column to use for sorting

        __ Description __
        sorts the supplied array
        """
        # 1 - sort array by the i-th column
        sorted_index = np.argsort(array_to_sort[column_to_sort_by])
        for i in range(0, len(array_to_sort)):
            array_to_sort[i] = np.array(array_to_sort[i])[sorted_index]
        return array_to_sort



        def import_transmission_spectrum(self, file_name, plot_axes, plot_list):
        """Fiel should have the format

        # no_of_x_col by no_of_y_col
        xCol(field) yCol(freq) zCol(tranmission)

        __ Parameters __


        plot_list: list of the field points to plot

        *** Format is ***
        # no-of-xCol by no-of-yCol
        xCol(field) yCol(freq) zCol(transmission)

        __ Description __
        Import the transmission array, supplied as a commented file. The comment
        must specify the number of field points used
        """

        logging.debug("Importing transmission data file '{file_name}'")

        with open(file_name) as fin:
            first_line = fin.readline()

        no_field_points = int(first_line.split()[1])
        no_freq_points = int(first_line.split()[3])
        logging.debug(
            f"{no_field_points} field points and {no_freq_points} frequency points"
        )

        transmission_array = np.vsplit(np.loadtxt(file_name), no_field_points)
        for i in range(0, len(transmission_array)):
            transmission_array[i] = transmission_array[i].transpose()

        if self.plot_or_not:
            logging.debug("Plotting, cos why not")
            for i in plot_list:
                if (i < 0) or (i >= len(transmission_array)):
                    raise ValueError(
                        "{i} is outside the allowed field values [0:{no_field_points}]"
                    )

                plot_axes.plot(transmission_array[i][1], transmission_array[i][2])
        logging.debug("==> 'import_transmission_spectrum' finished")

    def save_plot(self, mpl_axes, filename, dpi=100):
        mpl_axes.set_facecolor("white")
        plt.savefig(f"output/{filename}.png", dpi=dpi)
