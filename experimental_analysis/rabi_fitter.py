import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mp
from scipy.optimize import curve_fit

# TODO: Tidy this
class general_data:
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
