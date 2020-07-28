"""This simulator steps indendently across phi_r and phi_l
to evalute the transition energies into a grid
"""

from collections import defaultdict
from typing import Tuple, List, Dict
import logging

from pyprind import ProgBar
import numpy as np
from scipy.sparse.linalg import eigsh
import scipy.sparse as sp

import pinject


class TwinQubitSimulatorPhilPhir:
    @pinject.copy_args_to_public_fields
    def __init__(
        self,
        twin_qubit_hamiltonian_manager,
        twin_qubit_state_manager,
        twin_qubit_constant_manager,
        twin_qubit_operator_builder,
        quantum_constants,
    ):
        simulation_dictionary = defaultdict(list)

    def simulate(
        self,
        phi_l_list: List[float],
        phi_r_list: List[float],
        number_of_levels_to_simulate: int,
        phil_phir_coordinates_supplied: bool = False,
    ) -> Dict:
        """
        __ Parameters __
        [bool] phil_phir_coordinates_supplied:         whether the external flux arrays are passed
                                                       in according to the following mapping:
                        phi_l <-> phi_minus
                        phi_r <-> phi_plus

        in that case, a coordinate trasnformation will be required

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        """
        dim_l = len(phi_l_list)
        dim_r = len(phi_r_list)

        simulation_dictionary = {}
        simulation_dictionary["eigvals"] = np.empty(
            (dim_l, dim_r, number_of_levels_to_simulate)
        )
        simulation_dictionary["eigvecs"] = np.empty(
            (
                dim_l,
                dim_r,
                number_of_levels_to_simulate,
                self.twin_qubit_state_manager.states_total_number,
            )
        )
        simulation_dictionary["d0-1"] = np.empty((dim_l, dim_r, 1))
        simulation_dictionary["d0-2"] = np.empty((dim_l, dim_r, 1))

        self.twin_qubit_hamiltonian_manager.stage2_prepare_constant_hamiltonian()
        (voltage_matrix, phi_matrix) = self.twin_qubit_operator_builder.build()

        logging.info("💻 Running simulation")
        progress_bar = ProgBar(dim_l * dim_r, bar_char="█")
        for (phi_l_idx, phi_l) in enumerate(phi_l_list):
            for (phi_r_idx, phi_r) in enumerate(phi_r_list):
                if phil_phir_coordinates_supplied:
                    phi_l_adjusted = phi_r + phi_l
                    phi_r_adjusted = phi_r - phi_l
                else:
                    phi_l_adjusted = phi_l
                    phi_r_adjusted = phi_r

                self.twin_qubit_hamiltonian_manager.stage3_build_hamiltonian_for_simulation(
                    phi_l_adjusted, phi_r_adjusted
                )

                (eigvals, eigvecs) = eigsh(
                    self.twin_qubit_hamiltonian_manager.hamiltonian_simulation,
                    number_of_levels_to_simulate,
                    which="SA",
                    tol=0,
                )
                (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                    eigvals, eigvecs
                )

                simulation_dictionary = self.store_results(
                    simulation_dictionary, eigvals, eigvecs, phi_l_idx, phi_r_idx
                )

                # if evaluate_dipole_element:
                # simulation_dictionary = self.evaluate_dipole_element_and_append(eigvecs, voltage_matrix)

                progress_bar.update()

        logging.info("💻 Simulation completed")
        return simulation_dictionary

    @staticmethod
    def sort_in_ascending_eigval_order(
        eigvals: np.array, eigvecs: np.array
    ) -> Tuple[np.array, np.array]:
        sort_idx = np.argsort(eigvals)
        eigvals = np.array(eigvals)[sort_idx]
        eigvecs = np.transpose(eigvecs)[sort_idx]

        return (eigvals, eigvecs)

    def store_results(
        self,
        simulation_dictionary: Dict,
        eigvals: np.array,
        eigvecs: np.array,
        phi_l_idx: int,
        phi_r_idx: int,
    ) -> Dict:
        """Save the eigenvalues and eigenvector results to a dictionary
        """

        simulation_dictionary["eigvals"][phi_l_idx][phi_r_idx] = eigvals
        for (idx, vec) in enumerate(eigvecs):
            simulation_dictionary["eigvecs"][phi_l_idx][phi_r_idx][idx] = vec
        simulation_dictionary["0-1"][phi_l_idx][phi_r_idx] = eigvals[1] - eigvals[0]
        simulation_dictionary["1-2"][phi_l_idx][phi_r_idx] = eigvals[2] - eigvals[1]

        return simulation_dictionary

    def evaluate_dipole_element_and_append(
        self, eigvecs: np.array, voltage_matrix: sp.csr_matrix
    ) -> Dict:
        state1 = eigvecs[0]
        state2 = eigvecs[1]
        state3 = eigvecs[2]

        d01 = state2.dot(voltage_matrix.dot(state1))
        d32 = state3.dot(voltage_matrix.dot(state2))
        d13 = state1.dot(voltage_matrix.dot(state3))

        simulation_dictionary["d01"].append(np.abs(d01))
        simulation_dictionary["d32"].append(np.abs(d32))
        simulation_dictionary["d13"].append(np.abs(d13))

        return simulation_dictionary
        # simulation_dictionary["d01-beta"] = simulation_dictionary["d01-beta"].append(
        #     np.abs(d13 / (self.quantum_constants.phi0 * simulation_dictionary["0-1"][-1]))
        # )
