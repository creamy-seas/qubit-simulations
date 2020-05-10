from collections import defaultdict
from typing import Tuple
import logging

from pyprind import ProgBar
import numpy as np
from scipy.sparse.linalg import eigsh
import scipy.sparse as sp

import pinject


class TwinQubitSimulator:
    @pinject.copy_args_to_public_fields
    def __init__(
        self,
        twin_qubit_hamiltonian_manager,
        twin_qubit_state_manager,
        twin_qubit_constant_manager,
        flux_list,
        twin_qubit_operator_builder,
        quantum_constants,
    ):
        self.simulations = defaultdict(list)

    def evaluate_flux_in_loops(self, flux_number: float) -> Tuple[float, float]:
        phi_external = flux_number * 2 * np.pi
        phi_external_assymetric = (
            phi_external * self.twin_qubit_constant_manager.assymetry
        )

        return (phi_external, phi_external_assymetric)

    def simulate(self, evaluate_dipole_element=False):
        """
        __ Parameters __
        [bool] evaluate_dipole_element:   expectation values for the dipole

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        3 - plot out the spectrum
        """

        # ⦿ Reset for simulation
        self.simulations = defaultdict(list)
        self.simulations["eigvals"] = np.empty((0, 3))
        self.simulations["eigvecs"] = np.empty(
            (0, self.twin_qubit_state_manager.states_total_number)
        )
        self.twin_qubit_hamiltonian_manager.stage1_prepare_hamiltonian_skeleton()
        self.twin_qubit_hamiltonian_manager.stage2_prepare_constant_hamiltonian()
        (voltage_matrix, phi_matrix) = self.twin_qubit_operator_builder.build()

        logging.info("💻 Running simulation")
        progress_bar = ProgBar(len(self.flux_list), bar_char="█")
        for ext_flux_number in self.flux_list:

            (phi_external, phi_external_assymetric) = self.evaluate_flux_in_loops(
                ext_flux_number
            )

            self.twin_qubit_hamiltonian_manager.stage3_build_hamiltonian_for_simulation(
                phi_external, phi_external_assymetric
            )

            (eigvals, eigvecs) = eigsh(
                self.twin_qubit_hamiltonian_manager.hamiltonian_simulation,
                3,
                which="SA",
                tol=0,
            )
            (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(eigvals, eigvecs)

            self.append_results(eigvals, eigvecs)

            if evaluate_dipole_element:
                self.evaluate_dipole_element_and_append(eigvecs, voltage_matrix)

            progress_bar.update()

        logging.info("💻 Simulation completed")

    @staticmethod
    def sort_in_ascending_eigval_order(
        eigvals: np.array, eigvecs: np.array
    ) -> Tuple[np.array, np.array]:
        sort_idx = np.argsort(eigvals)
        eigvals = np.array(eigvals)[sort_idx]
        eigvecs = np.transpose(eigvecs)[sort_idx]

        return (eigvals, eigvecs)

    def append_results(self, eigvals: np.array, eigvecs: np.array):
        self.simulations["eigvals"] = np.vstack((self.simulations["eigvals"], eigvals))
        self.simulations["eigvecs"] = np.vstack((self.simulations["eigvecs"], eigvecs))
        self.simulations["1-2"].append(eigvals[1] - eigvals[0])
        self.simulations["2-3"].append(eigvals[2] - eigvals[1])

    def evaluate_dipole_element_and_append(
        self, eigvecs: np.array, voltage_matrix: sp.csr_matrix
    ):
        state1 = eigvecs[0]
        state2 = eigvecs[1]
        state3 = eigvecs[2]

        d21 = state2.dot(voltage_matrix.dot(state1))
        d32 = state3.dot(voltage_matrix.dot(state2))
        d13 = state1.dot(voltage_matrix.dot(state3))

        self.simulations["d21"].append(np.abs(d21))
        self.simulations["d32"].append(np.abs(d32))
        self.simulations["d13"].append(np.abs(d13))

        # self.simulations["d21-beta"] = self.simulations["d21-beta"].append(
        #     np.abs(d13 / (self.quantum_constants.phi0 * self.simulations["1-2"][-1]))
        # )
