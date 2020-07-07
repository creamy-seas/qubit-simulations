import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ilya_plot')


class qubit_master:
    """
    A set of common
    - quantum
    - plotting

    functions used by derived quantum classes
    """

    def __init__(self, plot_or_not, message_or_not):
        """
        __ Parameters __
        [bool] plot_or_not:     whether to prepare and run plotting functions
        [bool] message_or_not:  display information during run time

        __ Description __
        Initializes generic class that should be expanded for specific use-cases
        """

        """
        
        plot_or_not: whether to run the plotting functions
        """

        self.const_h = 6.64 * 10**(-34)
        self.const_eCharge = 1.6 * 10**(-19)
        self.const_Phi0 = self.const_h / (2 * self.const_eCharge)
        self.plot_or_not = plot_or_not
        self.message_or_not = message_or_not

    def prepare_plot(self, nrows, ncols):
        """
        __ Parameters __
        [int] nrows, ncols:     rows and columns of figures to prepare

        __ Description __
        Prepare class axes to plot on

        __ Return __
        [Figure] fig, [Axes] ax
        """

        # 1 - by default, nothing is plotting
        plt.ioff()
        plt.close("all")
        fig = None
        ax = None

        if(self.plot_or_not):
            if(self.message_or_not):
                print("==> 'prepare_plot' is setting up figure and axes")
            # 2 - interactive mode, to alow updating
            plt.ion()

            # 3 - define plots
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols)
            try:
                # 4  - adjust position on the screen
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(0, 30, 1280, 1600)
                fig.canvas.set_window_title('Run time data')

            except AttributeError:
                pass

            if(self.message_or_not):
                print("==> 'prepare_plot' finished")

        return fig, ax

    def convert_energy_to_GHz(self, energy):
        """
        __ Parameters __
        energy: in joules

        __ Description __
        Conversion from energy to GHz

        __ Return __
        Energy in GHz
        """
        return energy / (self.const_h * 10**(9))

    def energy_charging(self, area_nm2):
        """
        __ Parameters __
        area: JJ area in nm^2

        __ Description __
        Evaluating charging energy and capacitance
            (2e)^2/2C
            C = e*e_0*A/d

        __ Return __
        self.EC: the charging energy
        self.capacitance: capacitance for normalisation
        """

        # 1 - set constants
        thickness_alox = 2 * 10**(-9)
        permitivity_alox = 10
        permitivity_vacuum = 8.85 * 10**(-12)
        capacitance_area_offset = 1.15 * 10**(-18)  # convert nm2 to m2

        # 2 - evaluate and return
        self.param_capacitance = permitivity_alox * permitivity_vacuum * \
            capacitance_area_offset * area_nm2 / thickness_alox
        self.EC = (2 * self.const_eCharge)**2 / (2 * self.param_capacitance)
        self.EC = self.convert_energy_to_GHz(self.EC)

    def energy_josephson(self, squares_JJ):
        """
        __ Parameters __
        squares: width of the JJ in square units

        __ Description __
        Evaluate the Josephson energy, critical current, and resistance.
        The more squares of JJ, the bigger the resistance, and the larger it's energy

        We assume that it is a standard JJ: 20nm - AlOx - 30nm
        Count in 100x100nm^2 squares

        __ Return __
        self.EJ: josephson energy in GHz
        self.critical_current: critical current in A
        """

        # 1 - set constants
        transition_al = 1.2
        boltzman = 1.38 * 10**(-23)
        delta_al = 1.764 * transition_al * boltzman
        resistance_constant = self.const_h / (4 * self.const_eCharge**2)
        resistance_of100x100_Al = 18.4 * 10**3

        # 2 - evaluation
        self.param_resistance = resistance_of100x100_Al / squares_JJ
        self.param_critical_current = np.pi * delta_al / \
            (2 * self.const_eCharge * self.param_resistance)
        self.EJ = resistance_constant * delta_al / (2 * self.param_resistance)
        self.EJ = self.convert_energy_to_GHz(self.EJ)

    def data_transmission_slice(self, file_name, plot_axes, plot_list):
        """
        __ Parameters __
        [string] file_name:     file to load tranmission data from
        [Axes] plot_axes:       axes to perform plotting on

        plot_list: list of the field points to plot

        *** Format is ***
        # no-of-xCol by no-of-yCol
        xCol(field) yCol(freq) zCol(transmission)

        __ Description __
        Import the transmission array, supplied as a commented file. The comment
        must specify the number of field points used
        """

        print("==> 'data_transmission_slice' importing transmission data file")

        # 1 - check format
        with open(file_name) as fin:
            first_line = fin.readline()

        field_points = int(first_line.split()[1])
        freq_points = int(first_line.split()[3])

        # 2 - import and slice
        print("  > Importing file with\n\t%i field points\n\t%i frequency points" % (
            field_points, freq_points))

        transmission_array = np.vsplit(np.loadtxt(file_name), field_points)
        for i in range(0, len(transmission_array)):
            transmission_array[i] = transmission_array[i].transpose()

        # 3 - plot data
        print("  > Plotting data points")
        if(self.plot_or_not):
            for i in plot_list:
                if((i < 0) or (i >= len(transmission_array))):
                    self.raise_error(
                        "Invalid field points in 'field_points' parameter to data_transmission_slice function")
                plot_axes.plot(
                    transmission_array[i][1], transmission_array[i][2])
        print("==> 'data_transmission_slice' finished")

    def raise_error(self, display):
        output = "\n****************************************\n" + \
            display + "\n****************************************\n"
        print(output)
        raise ValueError(display)

    def save_plot(self, ax, filename, dpi=100):
        """
        __ Parameters __
        [Axes] ax:              to save
        [string] filename:      to save under
        [int] dpi:              save resolution

        __ Description __
        saves image with white background
        """

        ax.set_facecolor("white")

        plt.savefig("output/%s.png" % (filename), dpi=dpi)


if (__name__ == "__main__"):
    test = qubit_master(True)
    # test.experimental_data_load(test.ax, True)
    # test.simulate()
    # test.experimental_data_error()
    # a = list(np.arange(0, 17, 1))
    # test.data_transmission_slice(
    #     "qubit2_data/qubit2_m1_transmission.txt", test.ax, a)
