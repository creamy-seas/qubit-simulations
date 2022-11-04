from collections import defaultdict
from typing import Tuple
import logging
import itertools

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

    def simulate(
        self, number_of_levels_to_simulate: int, evaluate_dipole_element=False
    ):
        """
        __ Parameters __
        [bool] evaluate_dipole_element:   expectation values for the dipole

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        3 - plot out the spectrum
        """

        self.simulations = defaultdict(list)
        self.simulations["eigvals"] = np.empty((0, number_of_levels_to_simulate))
        self.simulations["eigvecs"] = np.empty(
            (0, self.twin_qubit_state_manager.states_total_number)
        )
        self.twin_qubit_hamiltonian_manager.stage2_prepare_constant_hamiltonian()
        (
            voltage_matrix,
            voltage_scaling_constant,
            _,
        ) = self.twin_qubit_operator_builder.build()

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
                number_of_levels_to_simulate,
                which="SA",
                tol=0,
            )
            (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(eigvals, eigvecs)

            self.append_results(eigvals, eigvecs)

            if evaluate_dipole_element:
                self.evaluate_dipole_element_and_append(
                    eigvecs, voltage_matrix, voltage_scaling_constant
                )

            progress_bar.update()

        logging.info("💻 Simulation completed")
        return self.simulations

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
        self,
        eigvecs: np.array,
        voltage_matrix: sp.csr_matrix,
        voltage_scaling_constant: float,
    ):
        for (i, j) in itertools.combinations(range(0, len(eigvecs)), 2):
            matrix_element = np.conjugate(eigvecs[i]).dot(
                voltage_matrix.dot(eigvecs[j])
            )
            self.simulations[f"d{i}-{j}"].append(
                np.abs(matrix_element) * voltage_scaling_constant
            )
