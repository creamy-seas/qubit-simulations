"""
Class that runs simulation - Hamiltonian should be prepared
"""

import logging
from collections import defaultdict
from typing import List, Dict

from pyprind import ProgBar
import numpy as np
import scipy

import pinject


class CqpsQubitSimulator:
    @pinject.copy_args_to_public_fields
    def __init__(self, cqps_qubit_hamiltonian_manager, cqps_qubit_constant_manager):
        pass

    def simulate(
        self,
        flux_ext_list: List[float],
        number_of_levels_to_simulate: int,
        use_sparse_matrix: bool = True,
    ) -> Dict:
        """
        __ Parameters __
        [list] flux_ext_list:           Flux in units of Φ₀
        [bool] use_sparse_matrix:       whether to use the approximated faster evaluation for
                                        eigenvalues and eigenvectors

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        """
        simulation_dictionary = defaultdict(list)

        dim = len(flux_ext_list)

        simulation_dictionary = {}
        simulation_dictionary["eigvals"] = np.empty(
            (dim, number_of_levels_to_simulate), dtype=float
        )
        simulation_dictionary["eigvecs"] = np.empty(
            (
                dim,
                number_of_levels_to_simulate,
                self.cqps_qubit_constant_manager.number_of_states,
            ),
            dtype=np.cdouble,
        )

        logging.info("💻 Running simulation")
        if logging.root.level == logging.INFO:
            progress_bar = ProgBar(dim, bar_char="█")

        for (flux_ext_idx, flux_ext) in enumerate(flux_ext_list):

            hamiltonian = self.cqps_qubit_hamiltonian_manager.prepare_hamiltonian(
                self.cqps_qubit_constant_manager.EL,
                self.cqps_qubit_constant_manager.ES,
                flux_ext,
            )

            if use_sparse_matrix:
                # Approximation - may not be valid ########################
                (eigvals, eigvecs) = scipy.sparse.linalg.eigsh(
                    hamiltonian,
                    number_of_levels_to_simulate,
                    which="SA",
                    tol=0,
                )
                (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                    eigvals, eigvecs
                )
            else:
                # Exact, but slow #########################################
                (eigvals, eigvecs) = scipy.linalg.eigh(
                    hamiltonian.todense(),
                )
                (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                    eigvals, eigvecs
                )
                eigvals = eigvals[:number_of_levels_to_simulate]
                eigvecs = eigvecs[:number_of_levels_to_simulate, :]

            simulation_dictionary = self.store_results(
                simulation_dictionary, eigvals, eigvecs, flux_ext_idx
            )

            if logging.root.level == logging.INFO:
                progress_bar.update()

        logging.info("💻 Simulation completed")
        return simulation_dictionary

    @staticmethod
    def sort_in_ascending_eigval_order(
        eigvals: np.array, eigvecs: np.array
    ) -> (np.array, np.array):
        sort_idx = np.argsort(eigvals)
        eigvals = np.array(eigvals)[sort_idx]
        eigvecs = np.transpose(eigvecs)[sort_idx]

        return (eigvals, eigvecs)

    @staticmethod
    def store_results(
        simulation_dictionary: Dict,
        eigvals: np.array,
        eigvecs: np.array,
        flux_ext_idx: int,
    ) -> Dict:
        """Save the eigenvalues and eigenvector results to a dictionary"""

        simulation_dictionary["eigvals"][flux_ext_idx] = np.real(eigvals)
        for (idx, vec) in enumerate(eigvecs):
            simulation_dictionary["eigvecs"][flux_ext_idx][idx] = vec

        return simulation_dictionary
