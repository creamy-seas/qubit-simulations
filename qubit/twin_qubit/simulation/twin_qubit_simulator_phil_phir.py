"""This simulator steps indendently across phi_2 and phi_1
to evalute the transition energies into a grid

Well, I learnt about a trap
scipy.linalg.eig and eigh are exact, while
scipy.sparse.linalg eigs and eigh are approximation - which are worse than the MatLab ones
"""
import itertools
from collections import defaultdict
from typing import Tuple, List, Dict
import logging

from pyprind import ProgBar
import numpy as np
import scipy
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
        phi_1_list: List[float],
        phi_2_list: List[float],
        number_of_levels_to_simulate: int,
        phil_phir_coordinates_supplied: bool = False,
        use_sparse_matrix: bool = True,
    ) -> Dict:
        """
        __ Parameters __
        [bool] phil_phir_coordinates_supplied:         if True:
                        phi_l --> phi_1
                        phi_r --> phi_2
                                                       if False
                        phi_plus  --> phi_1
                        phi_minus --> phi_2

        in that case, a coordinate trasnformation will be required
        [bool] use_sparse_matrix:                               whether to use the approximated faster evaluation for
                                                        eigenvalues and eigenvectors

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        """
        dim_1 = len(phi_1_list)
        dim_2 = len(phi_2_list)

        simulation_dictionary = defaultdict(lambda: np.empty((dim_1, dim_2, 1)))
        simulation_dictionary["eigvals"] = np.empty(
            (dim_1, dim_2, number_of_levels_to_simulate), dtype=float
        )
        simulation_dictionary["eigvecs"] = np.empty(
            (
                dim_1,
                dim_2,
                number_of_levels_to_simulate,
                self.twin_qubit_state_manager.states_total_number,
            ),
            dtype=np.cdouble,
        )

        self.twin_qubit_hamiltonian_manager.stage2_prepare_constant_hamiltonian()
        (voltage_matrix, phi_matrix) = self.twin_qubit_operator_builder.build()

        logging.info("💻 Running simulation")
        progress_bar = ProgBar(dim_1 * dim_2, bar_char="█")

        for (phi_1_idx, phi_1) in enumerate(phi_1_list):
            for (phi_2_idx, phi_2) in enumerate(phi_2_list):
                if phil_phir_coordinates_supplied == True:
                    phi_l = phi_1
                    phi_r = phi_2
                else:
                    phi_l = phi_1 + phi_2
                    phi_r = phi_1 - phi_2

                self.twin_qubit_hamiltonian_manager.stage3_build_hamiltonian_for_simulation(
                    phi_l, phi_r
                )

                ###############################################################
                #                    h stands for Hermetian                   #
                # Note - the two methods return vectors as row in once case and columns in the other
                ###############################################################
                if use_sparse_matrix:
                    # Approximation - may not be valid ########################
                    (eigvals, eigvecs) = scipy.sparse.linalg.eigsh(
                        self.twin_qubit_hamiltonian_manager.hamiltonian_simulation,
                        number_of_levels_to_simulate,
                        # which="SA",
                        which="SR",
                        tol=0,
                    )
                    (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                        eigvals, eigvecs
                    )
                else:
                    # Exact, but slow #########################################
                    (eigvals, eigvecs) = scipy.linalg.eigh(
                        self.twin_qubit_hamiltonian_manager.hamiltonian_simulation.todense(),
                    )
                    (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                        eigvals, eigvecs
                    )
                    eigvals = eigvals[:number_of_levels_to_simulate]
                    eigvecs = eigvecs[:number_of_levels_to_simulate, :]

                simulation_dictionary = self.store_results(
                    simulation_dictionary, eigvals, eigvecs, phi_1_idx, phi_2_idx
                )

                simulation_dictionary = self.evaluate_dipole_element_and_store(
                    simulation_dictionary, eigvecs, voltage_matrix, phi_1_idx, phi_2_idx
                )

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
        phi_1_idx: int,
        phi_2_idx: int,
    ) -> Dict:
        """Save the eigenvalues and eigenvector results to a dictionary"""

        simulation_dictionary["eigvals"][phi_1_idx][phi_2_idx] = np.real(eigvals)
        # for (idx, vec) in enumerate(eigvecs):
        #     simulation_dictionary["eigvecs"][phi_1_idx][phi_2_idx][idx] = vec
        # simulation_dictionary["0-1"][phi_1_idx][phi_2_idx] = eigvals[1] - eigvals[0]
        # simulation_dictionary["1-2"][phi_1_idx][phi_2_idx] = eigvals[2] - eigvals[1]

        return simulation_dictionary

    def evaluate_dipole_element_and_store(
        self,
        simulation_dictionary: Dict,
        eigvecs: np.array,
        voltage_matrix: sp.csr_matrix,
        phi_1_idx: int,
        phi_2_idx: int,
    ) -> Dict:

        for (i, j) in itertools.combinations(range(0, len(eigvecs)), 2):
            matrix_element = np.conjugate(eigvecs[i]).dot(
                voltage_matrix.dot(eigvecs[j])
            )
            simulation_dictionary[f"d{i}-{j}"][phi_1_idx][phi_2_idx] = np.abs(
                matrix_element
            )

        return simulation_dictionary
