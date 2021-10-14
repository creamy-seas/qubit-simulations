import logging

import matplotlib.pyplot as plt
from matplotlib.pyplot import Axes
import scipy.sparse as sp
import numpy as np

import pinject


class CQPSTwinQubitSparseMatrixVisualiser:
    """Visualise the sparse matrix"""

    @pinject.copy_args_to_public_fields
    def __init__(
        self, cqps_twin_qubit_hamiltonian_manager, cqps_twin_qubit_state_manager
    ):
        pass

    def generate_csr_matrix(self, hamiltonian_component: str) -> sp.csr_matrix:
        component = self.cqps_twin_qubit_hamiltonian_manager.hamiltonian_elements[
            hamiltonian_component
        ]

        return sp.coo_matrix(
            (np.ones(len(component["row"])), (component["row"], component["col"]))
        ).tocsr()

    def visualise_matrix(
        self,
        mpl_axes,
    ):
        MY_COLOURS = {
            "PINK": "#cc2375",
            "DeepPink4": "#8a154f",
            "DeepPink3": "#cc2375",
            "SeaGreen3": "#43cd80",
            "DodgerBlue4": "#104e8b",
            "DodgerBlue2": "#1c86ee",
            "DeepSkyBlue4": "#00688b",
            "DeepSkyBlue3": "#009acd",
            "DeepSkyBlue2": "#00bfff",
        }

        plt.rc("text.latex", preamble=r"\usepackage{braket}")

        diagonal_elements = self.generate_csr_matrix("inductance")
        loop_env_elements = self.generate_csr_matrix("loop-env-tunneling")
        loop_loop_elements = self.generate_csr_matrix("loop-loop-tunneling")

        mpl_axes.spy(diagonal_elements, color=MY_COLOURS["DeepSkyBlue3"], markersize=6)
        mpl_axes.spy(loop_env_elements, color=MY_COLOURS["SeaGreen3"], markersize=6)
        mpl_axes.spy(loop_loop_elements, color=MY_COLOURS["DeepPink3"], markersize=6)

        self.format_axes(mpl_axes)
        # self.add_ticks_to_axes(mpl_axes)

    def format_axes(self, mpl_axes: Axes):
        mpl_axes.grid(b=True, which="major", color="black")
        mpl_axes.grid(b=True, which="minor")

        states_total_number = self.cqps_twin_qubit_state_manager.states_total_number
        logging.info(f"📉 Using {states_total_number} for plotting")

        BORDER = 0.5
        mpl_axes.set_xlim([-BORDER, states_total_number - BORDER])
        mpl_axes.set_ylim([states_total_number - BORDER, -BORDER])

    def add_ticks_to_axes(self, mpl_axes: Axes):

        states_total_number = self.cqps_twin_qubit_state_manager.states_total_number

        mpl_axes.set_xticks(
            np.linspace(0, states_total_number - 1, states_total_number), minor=True
        )
        mpl_axes.set_yticks(
            np.linspace(0, states_total_number - 1, states_total_number), minor=True
        )

        mpl_axes.set_xticks([0, 13, 26])
        xticklabels = [""] * 3
        xticklabels[0] = r"$\left|-1, -1, -1 \right\rangle$"
        xticklabels[1] = r"$\left|0, 0, 0 \right\rangle$"
        xticklabels[2] = r"$\left|+1, +1, +1 \right\rangle$"

        mpl_axes.set_yticks([0, 6, 13, 19, 26])
        yticklabels = [""] * 5
        yticklabels[1] = r"$\left\langle -1, +1, -1 \right|$"
        yticklabels[2] = r"$\left\langle 0, 0, 0 \right|$"
        yticklabels[3] = r"$\left\langle +1, -1, 0 \right|$"
        yticklabels[4] = r"$\left\langle +1, +1, +1 \right|$"

        mpl_axes.set_xticklabels(xticklabels)
        mpl_axes.set_yticklabels(yticklabels)
