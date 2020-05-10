import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mp
from scipy.optimize import curve_fit

# TODO: Tidy this
class general_data:
    def transmission_load(self, colX, colY, convert_to_ghz):
        """
        __ Parameters __
        colX, colY: column to treat as the X,Y coordinate
        convert_to_ghz: nromalise by 10^6 or not
        convert_to_ghz: nromalise by 10^6 or not

        __ Description __
        Loads in the file to 'self.data_transmission'
        """
        # temp_load = np.loadtxt("data/extinction_22.txt").transpose()
        temp_load = np.loadtxt("data/extinction_smith27.txt").transpose()
        self.transmission_x = temp_load[colX, :]
        self.transmission_y = temp_load[colY, :]

        if convert_to_ghz:
            self.transmission_x = self.transmission_x / 10 ** 9

    def transmission_filter(self, x_min, x_max, y_min, y_max):
        if len(self.transmission_x) == 0:
            raise Exception("No data loaded!")

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

    def plot_2D(self, plot_axes, array_x, array_y, color):
        """
        __ Parameters __
        [array] array_x, array_y: data to plot
        [plt.Axes] plot_axes where to display data

        __ Description __
        plots data on the chosen axes
        """
        if self.plot_or_not:
            plot_axes.scatter(array_x, array_y, marker=".", color=color)
            plt.show()

    def transmission_plot(self, plot_axes, colour="C1"):
        """
        Plots loaded general_data data
        """
        if self.plot_or_not:
            plot_axes.scatter(
                self.transmission_x, self.transmission_y, marker=".", color=colour
            )
            plt.show()

    def __transmission_fit_function(self, x, Gamma1, Gamma2, Omega, offset):
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

    def transmission_fit(self, plot_axes, colour="C9"):
        # 1 - do fitting
        minTransmission = np.argmin(self.transmission_y)
        offset = self.transmission_x[minTransmission]

        popt, pcov = curve_fit(
            self.__transmission_fit_function,
            self.transmission_x,
            self.transmission_y,
            bounds=([0, 0, 0, offset - 0.1], [np.inf, np.inf, np.inf, offset + 0.1]),
        )

        print("  > Gamma1:\t%.4f" % (popt[0]))
        print("  > Gamma2:\t%.4f" % (popt[1]))
        print("  > Omega:\t%.4f" % (popt[2]))
        print("  > Offset:\t%.4f" % (popt[3]))

        # 2 - prepare arrays and plot
        temp_x = np.linspace(min(self.transmission_x), max(self.transmission_x), 500)
        temp_y = self.__transmission_fit_function(
            temp_x, popt[0], popt[1], popt[2], popt[3]
        )

        minTransmission = np.argmin(temp_y)
        xoffset = temp_x[minTransmission]
        yoffset = temp_y[minTransmission]

        if self.plot_or_not:
            plot_axes.plot(
                xoffset,
                yoffset,
                marker="o",
                color="#004BA8",
                markeredgecolor="C2",
                markersize=14,
                alpha=0.95,
                linestyle="",
            )
            # marker="o", color=colour, markeredgewidth=9, alpha=1)
            plot_axes.plot(temp_x, temp_y, color=colour)
            plot_axes.set_xlabel("$\omega_{21}/ 2 \pi$ (GHz)")
            plot_axes.set_ylabel("$|t|^2$")

            plt.show()

    def rabi_load(self, filename, colX, colY):
        """
        __ Parameters __
        [str] filename: where to load rabi data from
        [int] colX, colY: which columns to treat as the X (time) and Y(amplitude)

        __ Description __
        loads rabi oscillation data
        """

        temp_load = np.loadtxt(filename)
        self.rabi_x = temp_load[:, colX]
        self.rabi_y = temp_load[:, colY] * 10 ** 6

    def __rabi_fit_function(self, x, A, tDec, t_p, phi, D):
        """
        __ Description __
        Fits Rabi oscillations of the format
        A e^(-t/tDec) cos(2pi*t/tP+phi) + D
        """

        return A * np.sin(2 * np.pi * x / t_p + phi) * np.exp(-x / tDec) + D

    def rabi_fit(self, plot_axes, color):
        # 1 - do fitting

        popt, pcov = curve_fit(
            self.__rabi_fit_function,
            self.rabi_x,
            self.rabi_y,
            bounds=([0, 35, 5, -np.pi, -1e-5], [1, 45, 20, np.pi, 1e-5]),
        )

        print("  > Amplitude:\t%.3f" % (popt[0]))
        print("  > t_dec:\t%.3f" % (popt[1]))
        print("  > t_period:\t%.3f" % (popt[2]))
        print("  > phi_offset:\t%.3f" % (popt[3]))
        print("  > offset:\t%.3f" % (popt[4]))

        # 2 - prepare arrays and plot
        temp_x = np.linspace(min(self.rabi_x), max(self.rabi_x), 500)
        temp_y = self.__rabi_fit_function(
            temp_x, popt[0], popt[1], popt[2], popt[3], popt[4]
        )

        if self.plot_or_not:
            plot_axes.plot(temp_x, temp_y, color=color)
            plot_axes.set_xlabel("Pulse length, $\Delta t$ (ns)")
            plot_axes.set_ylabel("Real Amplitude (a.u)")

            plt.show()


if __name__ == "__main__":
    # ##########################################################
    # ################### Rabi  ################################
    ############################################################
    # utils.config_plot_size(0.2, 0.9, 0.15, 0.9)
    # mpl.rcParams["xtick.labelsize"] = 25
    # mpl.rcParams["ytick.labelsize"] = 25
    # mpl.rcParams["axes.labelsize"] = 25
    # test = general_data(True)
    # test.rabi_load("data/rabi_oscillation.txt", 1, 2)
    # test.plot_2D(test.ax, test.rabi_x, test.rabi_y, "C3")
    # test.rabi_fit(test.ax, "#7b68ee")
    # utils.save_ree(test.ax, "output/fig5_rabi", "svg")

    # ##########################################################
    # ################### Transmission #########################
    ############################################################
    utils.config_plot_size(0.25, 0.9, 0.28, 0.9)
    mpl.rcParams["text.color"] = (1, 1, 1)
    # mpl.rcParams["ytick.color"] = (0.53, 0.53, 1)
    # mpl.rcParams["axes.labelcolor"] = (0.53, 0.53, 1)
    mpl.rcParams["xtick.labelsize"] = 45
    mpl.rcParams["ytick.labelsize"] = 45
    mpl.rcParams["axes.labelsize"] = 45
    test = general_data(True)
    test.transmission_load(0, 1, True)
    test.transmission_plot(test.ax, colour="C2")
    test.transmission_filter(0.8, 10 ** 12, 0, 1.035)
    test.transmission_fit(test.ax, colour="C10")
    # test.ax.set_xticks([9.1, 9.2, 9.3, 9.4])
    # test.ax.set_ylim([5, 20])
    # test.ax.set_yticks([0.7, 0.8, 0.9, 1.0])
    utils.save_ree(test.ax, "output/fig2_transmission", "svg")
