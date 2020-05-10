import logging

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ilya_plot")

from common import plotter
from quantum_builder import QuantumBuilder

QUBIT_PARAMETERS = {
    "alpha": 1.023,
    "assymetry": 1.011,
    "jj_squares": 2,
    "states_per_island": 7,
}
FLUX_LIST = np.linspace(0.3, 0.7, 100)
EC = 13.5
EJ = 91

twin_qubit = QuantumBuilder.build_twin_qubit(
    QUBIT_PARAMETERS, FLUX_LIST, logging_level=logging.INFO
)
twin_qubit.override_parameters(EC=EC, EJ=EJ)
twin_qubit.run_simulation(evaluate_dipole_element=True)

fig = plt.figure(figsize=(6, 6))
ax = fig.subplots(nrows=1, ncols=1)
plotter.set_geometry(1000, 50, 800, 800)
# twin_qubit.plot_sparse_matrix(ax)
# ax[0].set_ylim([4, 20])
# twin_qubit.plot_transitions(ax[0])
# twin_qubit.plot_dipole_matrix_elements(ax)
ax.plot(twin_qubit.flux_list, twin_qubit.simulations["d32"])
# ax.set_ylim([0, 5 * 10 ** (-7)])
plt.show()
