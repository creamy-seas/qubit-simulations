"""
Class that runs simulation - Hamiltonian should be prepared
"""

import logging
from collections import defaultdict
from typing import List, Dict, Tuple

from pyprind import ProgBar
import numpy as np
import scipy

import pinject


class CqpsTwinQubitSimulator:
    @pinject.copy_args_to_public_fields
    def __init__(
        self, cqps_twin_qubit_hamiltonian_manager, cqps_twin_qubit_constant_manager
    ):
        pass

    def simulate(
        self,
        phi_dict: Dict,
        number_of_levels_to_simulate: int,
        use_sparse_matrix: bool = True,
    ) -> Dict:
        """
        __ Parameters __
        [list] phi_dict:                Information on the bias in the form of a dict
                If phi_dict['type'] = can be balanced, phi_l-phi_r or phi_+-phi_-
                If phi_dict['balanced'] = false, run simulation in each direction separately
        [bool] use_sparse_matrix:       whether to use the approximated faster evaluation for
                                        eigenvalues and eigenvectors

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        """

        # Prepare sweeping parameters
        (phi_l_list, phi_r_list) = self.unpack_phi_sweep_dictionary(phi_dict)
        dim_l = len(phi_l_list)
        dim_r = len(phi_r_list)

        simulation_dictionary = {
            "eigvals": np.empty(
                (dim_l, dim_r, number_of_levels_to_simulate), dtype=float
            ),
            "eigvecs": np.empty(
                (
                    dim_l,
                    dim_r,
                    number_of_levels_to_simulate,
                    self.cqps_twin_qubit_constant_manager.states_total_number,
                ),
                dtype=np.cdouble,
            ),
        }

        # Prepare constant part of the Hamiltonian
        self.cqps_twin_qubit_hamiltonian_manager.stage1_prepare_hamiltonian_skeleton()

        logging.info("💻 Running simulation")
        if logging.root.level == logging.INFO:
            progress_bar = ProgBar(dim_l * dim_r, bar_char="█")
        for (phi_l_idx, phi_l) in enumerate(phi_l_list):
            for (phi_r_idx, phi_r) in enumerate(phi_r_list):
                # In the case of balanced flux
                if phi_r is None:
                    phi_r = phi_l

                (eigvals, eigvecs) = self.simulation_kernel(
                    phi_l=phi_l,
                    phi_r=phi_r,
                    use_sparse_matrix=use_sparse_matrix,
                    number_of_levels_to_simulate=number_of_levels_to_simulate,
                )
                simulation_dictionary = self.store_results(
                    simulation_dictionary, eigvals, eigvecs, phi_l_idx, phi_r_idx
                )

                if logging.root.level == logging.INFO:
                    progress_bar.update()

        logging.info("💻 Simulation completed")
        return simulation_dictionary

    @staticmethod
    def unpack_phi_sweep_dictionary(
        phi_dict: Dict,
    ) -> Tuple[List[float], List[float]]:
        """
        Unpack the sweep parameters into a dict of phi_l and phi_r that can be used in simulatio
        """
        try:
            if "type" not in phi_dict:
                raise ValueError()

            if phi_dict["type"] == "balanced":
                phi_l_list = phi_dict["phi_list"]
                phi_r_list = [None]

            elif phi_dict["type"] == "phi_l-phi_r":
                phi_l_list = phi_dict["phi_l_list"]
                phi_r_list = phi_dict["phi_r_list"]

            elif phi_dict["type"] == "phi_+-phi_-":
                phi_l_list = []
                phi_r_list = []
                for (plus, minus) in zip(
                    phi_dict["phi_+_list"], phi_dict["phi_-_list"]
                ):
                    phi_l_list.append(plus + minus)
                    phi_r_list.append(plus - minus)

            else:
                raise ValueError()
        except KeyError as err:
            raise KeyError(
                f"Please supply {err} to the phi_dict to be used for {phi_dict['type']} simulation"
            )
        except ValueError:
            raise ValueError(
                """Please supply 'type' of flux sweep! Possible options:
                       - 'balanced'
                       - 'phi_l-phi_r'
                       - 'phi_+-phi_-'
                """
            )
        return (
            phi_l_list,
            phi_r_list,
        )

    def simulation_kernel(
        self,
        phi_l: float,
        phi_r: float,
        use_sparse_matrix: bool,
        number_of_levels_to_simulate: int,
    ) -> Tuple[np.array, np.array]:
        """
        Run simulation for supplied bias, returning the eigenvectors and eigenvalues,
        ordered by ascending eigenvalues
        """

        hamiltonian = self.cqps_twin_qubit_hamiltonian_manager.stage_2_build_hamiltonian_for_simulation(
            phi_l=phi_l, phi_r=phi_r
        )

        if use_sparse_matrix:
            # Approximation - may not be valid ########################
            (eigvals, eigvecs) = scipy.sparse.linalg.eigsh(
                hamiltonian,
                number_of_levels_to_simulate,
                which="SA",
                tol=0,
            )
            (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(eigvals, eigvecs)
        else:
            # Exact, but slow #########################################
            (eigvals, eigvecs) = scipy.linalg.eigh(
                hamiltonian.todense(),
            )
            (eigvals, eigvecs) = self.sort_in_ascending_eigval_order(eigvals, eigvecs)
            eigvals = eigvals[:number_of_levels_to_simulate]
            eigvecs = eigvecs[:number_of_levels_to_simulate, :]

        return (eigvals, eigvecs)

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
        phi_l_idx: int,
        phi_r_idx: int,
    ) -> Dict:
        """Save the eigenvalues and eigenvector results to a dictionary"""

        simulation_dictionary["eigvals"][phi_l_idx, phi_r_idx] = np.real(eigvals)
        for (idx, vec) in enumerate(eigvecs):
            simulation_dictionary["eigvecs"][phi_l_idx, phi_r_idx][idx] = vec

        return simulation_dictionary
