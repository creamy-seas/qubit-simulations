import sys
sys.path.append('/Users/CCCP/creamy_seas/monday_starts_on_saturday/dipole/qubit_simulations/qubit_simulation')

import matplotlib.pyplot as plt
import honkler
from qubit_twin import twin
from qubit_flux import flux
import matplotlib as mpl
import numpy as np

plt.style.use('ilya_plot')

# ############################################################
# ################### Paper Plot Inset ####################
# ############################################################
# test = twin(1, 1, 3, 100, True, False)
# test.sparse_matrix_visualise()
# honkler.save_ree(test.ax, "output/fig4_matrix", "svg")

# ############################################################
# ################### Paper Plot Inset ####################
# ############################################################
# 1 - colouring and size
# mpl.rcParams["xtick.color"] = (0.53, 0.53, 1)
# mpl.rcParams["ytick.color"] = (0.53, 0.53, 1)
# mpl.rcParams["axes.labelcolor"] = (0.53, 0.53, 1)
# honkler.config_plot_size(0.2, 0.8, 0.2, 0.9)
# mpl.rcParams["xtick.labelsize"] = 16

# # 2 - load spectrum and simulation
# EC = 13.5
# EJ = 92
# alpha = 1.023
# assymetry = 1.011
# no_states = 7
# flux_points = 100

# test = twin(alpha, assymetry, no_states, flux_points, True, True, cut_interaction=1)
# test.experimental_data_load(test.ax, True)
# test.prepare_operators()
# test.simulate([False, False])
# test.plot_simulation(test.ax)


# # 3 - limits and labels
# # test.ax.set_xlim([0.4, 0.6])
# # test.ax.set_xticks([0.4, 0.45, 0.5, 0.55, 0.6])
# test.ax.set_ylim([5, 20])
# # test.ax.set_yticks([0, 10, 20])
# test.ax.set_xlabel("Magnetic Flux ($\Phi$)")
# test.ax.set_ylabel("$\omega/2\pi$ (GHz)")

# honkler.save_ree(test.ax, "output/fig2_spectrum", "svg")

############################################################
################### Dipole transition ######################
############################################################
mpl.rcParams["xtick.labelsize"] = 25
mpl.rcParams["ytick.labelsize"] = 25
mpl.rcParams["axes.labelsize"] = 25
no_points = 300
honkler.config_plot_size(0.2, 0.9, 0.15, 0.9)

# 1 - twin qubit ##############################################################
EC = 13.5
EJ = 92
alpha = 1.023
assymetry = 1
twinclass = twin(alpha, assymetry, 15, no_points, True, False, cut_interaction=1)
twinclass.prepare_operators()
twinclass.override_parameters(EC, EJ, alpha, assymetry)
twinclass.simulate([True, False])

twinclass.ax.plot(twinclass.flux_list,
                  (np.array(twinclass.dipole_moment_voltage['01'])[:, 0]**2 +
                   np.array(twinclass.dipole_moment_voltage['01'])[:, 1]**2)**(1 / 2) * 1000000,
                  label="1<->2", color='#004BA8')
twinclass.ax.plot(twinclass.flux_list,
                  (np.array(twinclass.dipole_moment_voltage['02'])[:, 0]**2 +
                   np.array(twinclass.dipole_moment_voltage['02'])[:, 1]**2)**(1 / 2) * 1000000,
                  label="1<->3", color='C6')
twinclass.ax.plot(twinclass.flux_list,
                  (np.array(twinclass.dipole_moment_voltage['12'])[:, 0]**2 +
                   np.array(twinclass.dipole_moment_voltage['12'])[:, 1]**2)**(1 / 2) * 1000000,
                  label="2<->3", color='C4')

twinclass.ax.set_xlabel("Normalized Magnetic Flux ($\Phi/\Phi_0$)")
twinclass.ax.set_ylabel("Matrix elements ($\mu$V)")

twinclass.ax.set_xlim([0.365, 0.635])
twinclass.ax.set_ylim([0, 50])
honkler.save_ree(twinclass.ax, "fig4_dipole", "svg")


# 2 - flux qubit ##############################################################
# no_points = 10000
# EC = 20
# EJ = 30
# alpha = 0.45
# fluxclass = flux(alpha, 7, no_points, True, False)
# fluxclass.prepare_operators()
# fluxclass.override_parameters(EC, EJ, alpha)
# fluxclass.simulate([True, False])

# 3 - plotting
# fig, axTwin = plt.subplots(nrows=1, ncols=1)
# mngr = plt.get_current_fig_manager()
# mngr.window.setGeometry(0, 30, 1280, 1600)
# fig.canvas.set_window_title('Run time data')
# axFlux = fig.add_subplot(111, sharex=axTwin, frameon=False)
# axFlux.yaxis.tick_right()
# axFlux.yaxis.set_label_position("right")

# b - plot
# twinclass.plot_dipole_moment_voltage_beta(axTwin)

# c - coloring
# axTwin.spines['right'].set_color('C6')
# axTwin.yaxis.label.set_color('C6')
# axTwin.tick_params(axis='y', colors='C6')
# axFlux.spines['right'].set_color('C9')
# axFlux.yaxis.label.set_color('C9')
# axFlux.tick_params(axis='y', colors='C9')
