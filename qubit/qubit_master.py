import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes

from utils.logger import setup_logging


class MasterQubit:
    """
    """

    def __init__(self, plot_or_not, logging_level: int):
        """
        __ Parameters __
        [bool] plot_or_not:     whether to prepare and run plotting functions
        """

        setup_logging(logging_level)

        self.const_h = 6.64 * 10 ** (-34)
        self.const_eCharge = 1.6 * 10 ** (-19)
        self.const_Phi0 = self.const_h / (2 * self.const_eCharge)

        self.plot_or_not = plot_or_not

    def prepare_plot(self, nrows, ncols) -> Tuple[Figure, Axes]:

        plt.ioff()
        plt.close("all")
        fig = None
        mpl_axes = None

        if self.plot_or_not:
            logging.info("==> 'prepare_plot' is setting up figure and axes")

            # Interactive mode, to alow updating
            plt.ion()

            fig, mpl_axes = plt.subplots(nrows=nrows, ncols=ncols)
            try:
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(0, 30, 1280, 1600)
                fig.canvas.set_window_title("Run time data")

            except AttributeError:
                pass

        return fig, mpl_axes

    def convert_energy_to_GHz(self, energy):
        return energy / (self.const_h * 10 ** (9))

    def set_charging_energy(self, area_nm2):
        """
        __ Parameters __
        area: JJ area in nm^2

        __ Description __
        Evaluating charging energy and capacitance
            (2e)^2/2C
            C = e*e_0*A/d

        __ Sets __
        self.EC:                the charging energy
        self.capacitance:       capacitance for normalisation
        """

        THICKNESS_ALOX = 2 * 10 ** (-9)
        PERMITIVITY_ALOX = 10
        PERMITIVITY_VACUUM = 8.85 * 10 ** (-12)
        CAPACITANCE_AREA_OFFSET = 1.15 * 10 ** (-18)  # convert nm2 to m2

        # 2 - evaluate and return
        self.param_capacitance = (
            PERMITIVITY_ALOX
            * PERMITIVITY_VACUUM
            * CAPACITANCE_AREA_OFFSET
            * area_nm2
            / THICKNESS_ALOX
        )
        self.EC = (2 * self.const_eCharge) ** 2 / (2 * self.param_capacitance)
        self.EC = self.convert_energy_to_GHz(self.EC)

    def set_josephson_energy(self, squares_JJ):
        """
        __ Parameters __
        squares: width of the JJ in square units

        __ Description __
        Evaluate the Josephson energy, critical current, and resistance.
        The more squares of JJ, the bigger the resistance, and the larger it's energy

        We assume that it is a standard JJ: 20nm - AlOx - 30nm
        Count in 100x100nm^2 squares

        __ Sets __
        self.EJ:                josephson energy in GHz
        self.critical_current:  critical current in A
        """

        TRANSITION_AL = 1.2
        BOLTZMAN = 1.38 * 10 ** (-23)
        DELTA_AL = 1.764 * TRANSITION_AL * BOLTZMAN
        RESISTANCE_CONSTANT = self.const_h / (4 * self.const_eCharge ** 2)
        RESISTANCE_OF_100x100_AL = 18.4 * 10 ** 3

        self.param_resistance = RESISTANCE_OF_100x100_AL / squares_JJ
        self.param_critical_current = (
            np.pi * DELTA_AL / (2 * self.const_eCharge * self.param_resistance)
        )
        self.EJ = RESISTANCE_CONSTANT * DELTA_AL / (2 * self.param_resistance)
        self.EJ = self.convert_energy_to_GHz(self.EJ)

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

        logging.info("Importing transmission data file '{file_name}'")

        with open(file_name) as fin:
            first_line = fin.readline()

        no_field_points = int(first_line.split()[1])
        no_freq_points = int(first_line.split()[3])
        logging.info(
            f"{no_field_points} field points and {no_freq_points} frequency points"
        )

        transmission_array = np.vsplit(np.loadtxt(file_name), no_field_points)
        for i in range(0, len(transmission_array)):
            transmission_array[i] = transmission_array[i].transpose()

        if self.plot_or_not:
            logging.info("Plotting, cos why not")
            for i in plot_list:
                if (i < 0) or (i >= len(transmission_array)):
                    raise ValueError(
                        "{i} is outside the allowed field values [0:{no_field_points}]"
                    )

                plot_axes.plot(transmission_array[i][1], transmission_array[i][2])
        logging.info("==> 'import_transmission_spectrum' finished")

    def save_plot(self, mpl_axes, filename, dpi=100):
        mpl_axes.set_facecolor("white")
        plt.savefig(f"output/{filename}.png", dpi=dpi)
