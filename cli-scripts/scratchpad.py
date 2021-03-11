import time
import matplotlib.pyplot as plt
import matplotlib as mpl

import qutils

plt.style.use("ilya_plot")


if __name__ == "__main__":

    ###########################################################################
    #                               Fig 2 insert                              #
    ###########################################################################
    # qutils.config_plot_size(0.2, 0.9, 0.15, 0.9)
    # mpl.rcParams["xtick.labelsize"] = 25
    # mpl.rcParams["ytick.labelsize"] = 25
    # mpl.rcParams["axes.labelsize"] = 25
    # EC = 13.5
    # EJ = 92
    # alpha = 1.023
    # assymetry = 1.011
    # test = twin(alpha, assymetry, 7, 5000, True, False)
    # test.experimental_data_load(test.ax, False)
    # test.prepare_operators()
    # test.override_parameters(EC, EJ, alpha, assymetry)
    # test.simulate([True, False])
    # test.plot_simulation(test.ax)
    # test.ax.set_xlim([0.4, 0.6])
    # test.ax.set_xticks([0.45, 0.5, 0.55])
    # test.ax.set_ylim([5, 20])
    # test.ax.set_xlabel("Normalized Magnetic Flux ($\Phi/\Phi_0$)")
    # test.ax.set_ylabel("$\omega/2\pi$ (GHz)")
    # qutils.save_ree(test.ax, "output/fig2_zoomed_spectrum", "svg")

    ###########################################################################
    #                                  Fig 3                                  #
    ###########################################################################
    qutils.config_plot_size(0.2, 0.9, 0.15, 0.9)
    mpl.rcParams["xtick.labelsize"] = 25
    mpl.rcParams["ytick.labelsize"] = 25
    mpl.rcParams["axes.labelsize"] = 25
    EC = 13.5
    EJ = 92
    alpha = 1.023
    assymetry = 1
    test = twin(alpha, assymetry, 7, 1000, True, False)
    test.plot_or_not = False
    test.experimental_data_load(None, False)
    test.prepare_operators()
    test.override_parameters(EC, EJ, alpha, assymetry)
    test.plot_or_not = True
    test.simulate([True, False])
    test.plot_simulation(test.ax)
    test.ax.set_xlim([0.365, 1 - 0.365])
    # test.ax.set_xticks([0.45, 0.5, 0.55])
    test.ax.set_ylim([3, 20])
    test.ax.set_xlabel("Normalized Magnetic Flux ($\Phi/\Phi_0$)")
    test.ax.set_ylabel("$\omega/2\pi$ (GHz)")
    qutils.save_ree(test.ax, "output/fig3_simulation", "svg")

    # ############################################################
    # ################### Simulation error analysis ##############
    # ############################################################
    # fig, ax = plt.subplots(nrows=2, ncols=2)
    # plt.ion()
    # plt.rcParams["agg.path.chunksize"] = 10000000
    # array_in = np.loadtxt("output/simulation_error_18apr2019.txt").transpose()
    # minValue = np.argmin(array_in[4])
    # for i in range(0, 5):
    #     print(array_in[i][minValue])

    # # 1 - filtering
    # list_0 = []
    # list_1 = []
    # list_2 = []
    # list_3 = []
    # list_4 = []

    # for i in range(0, len(array_in[0])):
    #     if (array_in[4][i] < 500):
    #         list_0.append(array_in[0][i])
    #         list_1.append(array_in[1][i])
    #         list_2.append(array_in[2][i])
    #         list_3.append(array_in[3][i])
    #         list_4.append(array_in[4][i])

    # array_out = np.array([list_0, list_1, list_2, list_3, list_4])

    # # 2 - pltting
    # # ax.plot(array_out[4])
    # ax[0][0].scatter(array_out[0], array_out[4], marker=",", s=1)
    # ax[0][0].set_title("EC")
    # ax[0][1].scatter(array_out[1], array_out[4], marker=",", s=1)
    # ax[0][1].set_title("EJ")
    # ax[1][0].scatter(array_out[2], array_out[4], marker=",", s=1)
    # ax[1][0].set_title("alpha")
    # ax[1][1].scatter(array_out[3], array_out[4], marker=",", s=1)
    # ax[1][1].set_title("assymetry")

    # # name = "output/simulation_error_17apr2019_comb"
    # # plt.savefig("%s.png" % name)
    # plt.show()

    # qutils.plot_column_data(ax, "output/simulation_error_16apr2019.txt")

    # ############################################################
    # ################### Dipole transition ######################
    # ############################################################
    # # 1 - twin qubit
    # EC = 13.5
    # EJ = 92
    # alpha = 1.023
    # assymetry = 1.011
    # twinclass = twin(alpha, assymetry, 7, no_points, True, False)
    # twinclass.prepare_operators()
    # twinclass.override_parameters(EC, EJ, alpha, assymetry)
    # twinclass.simulate([True, False])

    # # 2 - flux qubit
    # no_points = 10000
    # EC = 20
    # EJ = 30
    # alpha = 0.45
    # fluxclass = flux(alpha, 7, no_points, True, False)
    # fluxclass.prepare_operators()
    # fluxclass.override_parameters(EC, EJ, alpha)
    # fluxclass.simulate([True, False])

    # # 3 - plotting
    # # a - axes setup
    # fig, axTwin = plt.subplots(nrows=1, ncols=1)
    # mngr = plt.get_current_fig_manager()
    # mngr.window.setGeometry(0, 30, 1280, 1600)
    # fig.canvas.set_window_title('Run time data')
    # axFlux = fig.add_subplot(111, sharex=axTwin, frameon=False)
    # axFlux.yaxis.tick_right()
    # axFlux.yaxis.set_label_position("right")

    # # b - plot
    # fluxclass.plot_dipole_moment_voltage_beta(axFlux)
    # twinclass.plot_dipole_moment_voltage_beta(axTwin)

    # # c - coloring
    # # axTwin.spines['right'].set_color('C6')
    # # axTwin.yaxis.label.set_color('C6')
    # # axTwin.tick_params(axis='y', colors='C6')
    # # axFlux.spines['right'].set_color('C9')
    # # axFlux.yaxis.label.set_color('C9')
    # # axFlux.tick_params(axis='y', colors='C9')

    # # d - scaling
    # axTwin.set_ylim([0, 1e9])
    # axFlux.set_ylim([0, 1e9])

    # # e - saving
    # qutils.save_ree(axTwin, "output/fig6_dipoleBeta", "svg")
