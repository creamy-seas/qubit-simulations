"""
Class that runs simulation
"""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple

from pyprind import ProgBar
import numpy as np
import scipy
import scipy.sparse as sp

import pinject


class TransmonQubitSimulator:
    @pinject.copy_args_to_public_fields
    def __init__(
        self, transmon_qubit_hamiltonian_manager, transmon_qubit_constant_manager
    ):
        pass

    def simulate(
        self,
        N_ext_list: List[float],
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

        dim_1 = len(N_ext_list)
        dim_2 = len(flux_ext_list)

        simulation_dictionary = defaultdict(lambda: np.empty((dim_1, dim_2, 1)))
        simulation_dictionary["eigvals"] = np.empty(
            (dim_1, dim_2, number_of_levels_to_simulate), dtype=float
        )
        simulation_dictionary["eigvecs"] = np.empty(
            (
                dim_1,
                dim_2,
                number_of_levels_to_simulate,
                self.transmon_qubit_constant_manager.number_of_charge_states,
            ),
            dtype=np.cdouble,
        )

        logging.info("💻 Running simulation")
        if logging.root.level == logging.INFO:
            progress_bar = ProgBar(dim_1 * dim_2, bar_char="█")

        for (N_ext_idx, N_ext) in enumerate(N_ext_list):
            for (flux_ext_idx, flux_ext) in enumerate(flux_ext_list):

                EJ = self.transmon_qubit_constant_manager.evaluate_EJ_from_external_flux(
                    flux_ext
                )
                hamiltonian = self.transmon_qubit_hamiltonian_manager.prepare_hamiltonian(
                    self.transmon_qubit_constant_manager.EC, EJ, N_ext
                )

                if use_sparse_matrix:
                    # Approximation - may not be valid ########################
                    (eigvals, eigvecs) = scipy.sparse.linalg.eigsh(
                        hamiltonian, number_of_levels_to_simulate, which="SR", tol=0,
                    )
                    (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                        eigvals, eigvecs
                    )
                else:
                    # Exact, but slow #########################################
                    (eigvals, eigvecs) = scipy.linalg.eigh(hamiltonian.todense(),)
                    (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(
                        eigvals, eigvecs
                    )
                    eigvals = eigvals[:number_of_levels_to_simulate]
                    eigvecs = eigvecs[:number_of_levels_to_simulate, :]

                simulation_dictionary = self.store_results(
                    simulation_dictionary, eigvals, eigvecs, N_ext_idx, flux_ext_idx
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

    def store_results(
        self,
        simulation_dictionary: Dict,
        eigvals: np.array,
        eigvecs: np.array,
        N_ext_idx: int,
        flux_ext_idx: int,
    ) -> Dict:
        """Save the eigenvalues and eigenvector results to a dictionary"""

        simulation_dictionary["eigvals"][N_ext_idx][flux_ext_idx] = np.real(eigvals)
        for (idx, vec) in enumerate(eigvecs):
            simulation_dictionary["eigvecs"][N_ext_idx][flux_ext_idx][idx] = vec

        return simulation_dictionary
