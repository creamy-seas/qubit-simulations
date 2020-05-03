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
        self.simulation = defaultdict(list)

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
        self.simulation["eigvals"] = np.empty((0, 3))
        self.simulation["eigvecs"] = np.empty(
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
            (eigvals, eigvecs) = self.sort_in_ascending_order(eigvals, eigvecs)

            self.append_results(eigvals, eigvecs)

            if evaluate_dipole_element:
                self.evaluate_dipole_element_and_append(voltage_matrix)

            progress_bar.update()

        logging.info("💻 Simulation completed")

    @staticmethod
    def sort_in_ascending_order(
        eigvals: np.array, eigvecs: np.array
    ) -> Tuple[np.array, np.array]:
        sort_idx = np.argsort(eigvals)
        eigvals = np.array(eigvals)[sort_idx]
        eigvecs = np.transpose(eigvecs)[sort_idx]

        return (eigvals, eigvecs)

    def append_results(self, eigvals: np.array, eigvecs: np.array):
        self.simulation["eigvals"] = np.vstack((self.simulation["eigvals"], eigvals))
        self.simulation["eigvecs"] = np.vstack((self.simulation["eigvecs"], eigvecs))
        self.simulation["1-2"].append(eigvals[1] - eigvals[0])
        self.simulation["2-3"].append(eigvals[2] - eigvals[1])

    def evaluate_dipole_element_and_append(self, voltage_matrix: sp.csr_matrix):
        state_0 = self.simulation["eigvecs"][0]
        state_1 = self.simulation["eigvecs"][1]

        dipole_moment_voltage = state_0.dot(voltage_matrix.dot(state_1))

        dipole_moment_beta = dipole_moment_voltage / (
            self.quantum_constants.phi0 * self.simulation["1-2"][-1]
        )

        self.simulation["dipole-voltage"].append(np.abs(dipole_moment_voltage))
        self.simulation["dipole-beta"].append(np.abs(dipole_moment_beta))
