import logging

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ilya_plot")

from quantum_builder import QuantumBuilder

QUBIT_PARAMETERS = {
    "alpha": 1.023,
    "assymetry": 1.011,
    "jj_squares": 2,
    "states_per_island": 17,
}
FLUX_LIST = np.linspace(0.3, 0.7, 100)

twin_qubit = QuantumBuilder.build_twin_qubit(
    QUBIT_PARAMETERS, FLUX_LIST, logging_level=logging.INFO
)
twin_qubit.override_parameters(EC=13.5, EJ=91)
twin_qubit.run_simulation(False)

fig = plt.figure(figsize=(6, 6))
ax = fig.subplots(nrows=1, ncols=1)
# twin_qubit.plot_sparse_matrix(ax)
ax.set_ylim([4, 20])
twin_qubit.plot_transitions(ax)
plt.show()
