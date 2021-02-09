"""
Matrix is built freshly each time - since small matrices
are used, this is the simpler and more direct approach

⋱         ⁞               ⁞               ⁞               ⋰
...        EL(N - f)       -ES/2           0               ...
...        -ES/2           EL(N+1 - f)     -ES/2           ...
...        0               -ES/2           EL(N+2 - f)     ...
⋰         ⁞               ⁞               ⁞                 ⋱
"""

from collections import defaultdict

import pinject
import numpy as np
import scipy.sparse as sp


class CqpsQubitHamiltonianManager:
    @pinject.copy_args_to_public_fields
    def __init__(self, cqps_qubit_constant_manager):
        pass

    def prepare_hamiltonian(self, EL: float, ES: float, f_ext: float) -> "coo Matrix":

        hamiltonian_skeleton = defaultdict(lambda: defaultdict(list))
        number_of_states = self.cqps_qubit_constant_manager.number_of_states

        # Evalaute elements of the matrix #####################################
        for row in range(number_of_states):
            hamiltonian_skeleton["charge"]["row"] += [row]
            hamiltonian_skeleton["charge"]["col"] += [row]
            hamiltonian_skeleton["charge"]["elm"] += [
                (
                    EL
                    * (
                        self.cqps_qubit_constant_manager.convert_index_to_state(row)
                        - f_ext
                    )
                    ** 2
                )
            ]

            if row < number_of_states - 1:
                col = row + 1
                hamiltonian_skeleton["flux"]["row"] += [row, col]
                hamiltonian_skeleton["flux"]["col"] += [col, row]
                hamiltonian_skeleton["flux"]["elm"] += [-ES / 2, -ES / 2]

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
