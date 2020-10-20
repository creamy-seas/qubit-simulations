"""
Matrix is built freshly each time - since small matrices
are used, this is the simpler and more direct approach
"""

from collections import defaultdict

import pinject
import numpy as np
import scipy.sparse as sp


class TransmonQubitHamiltonianManager:
    @pinject.copy_args_to_public_fields
    def __init__(self, transmon_qubit_constant_manager):
        pass

    def prepare_hamiltonian(self, EC: float, EJ: float, N_ext: float) -> "coo Matrix":

        hamiltonian_skeleton = defaultdict(lambda: defaultdict(list))
        number_of_charge_states = (
            self.transmon_qubit_constant_manager.number_of_charge_states
        )

        # Evalaute elements of the matrix #####################################
        for row in range(number_of_charge_states):
            hamiltonian_skeleton["charge"]["row"] += [row]
            hamiltonian_skeleton["charge"]["col"] += [row]
            hamiltonian_skeleton["charge"]["elm"] += [
                (
                    4 * EC
                    * (
                        self.transmon_qubit_constant_manager.convert_index_to_state(row)
                        - N_ext
                    )
                    ** 2
                )
            ]

            if row < number_of_charge_states - 1:
                col = row + 1
                hamiltonian_skeleton["flux"]["row"] += [row, col]
                hamiltonian_skeleton["flux"]["col"] += [col, row]
                hamiltonian_skeleton["flux"]["elm"] += [-EJ / 2, -EJ / 2]

        # Build the matrix ####################################################
        return sp.coo_matrix(
            (
                np.hstack(
                    (
                        hamiltonian_skeleton["charge"]["elm"],
                        hamiltonian_skeleton["flux"]["elm"],
                    )
                ),
                (
                    np.hstack(
                        (
                            hamiltonian_skeleton["charge"]["row"],
                            hamiltonian_skeleton["flux"]["row"],
                        )
                    ),
                    np.hstack(
                        (
                            hamiltonian_skeleton["charge"]["col"],
                            hamiltonian_skeleton["flux"]["col"],
                        )
                    ),
                ),
            )
        ).tocsr()
