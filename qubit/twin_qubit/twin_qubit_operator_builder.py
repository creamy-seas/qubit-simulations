import logging
from typing import Tuple
from collections import defaultdict

import scipy.sparse as sp
import pinject
import numpy as np


class TwinQubitOperatorBuilder:
    """Builds the voltage and phase operators as np.csr_matrix"""

    @pinject.copy_args_to_public_fields
    def __init__(
        self, quantum_constants, twin_qubit_constant_manager, twin_qubit_state_manager
    ):
        pass

    def build(self) -> Tuple[sp.csr_matrix, sp.csr_matrix]:
        """Returns-> Builds voltage and phase operators """

        EC = self.twin_qubit_constant_manager.EC
        alpha = self.twin_qubit_constant_manager.alpha
        states_total_number = self.twin_qubit_state_manager.states_total_number
        states_per_island = self.twin_qubit_state_manager.states_per_island
        convert_island_state_to_index = (
            self.twin_qubit_state_manager.convert_island_state_to_index
        )
        convert_index_to_island_state = (
            self.twin_qubit_state_manager.convert_index_to_island_state
        )

        logging.debug(
            f"""🏭 Generating operators with
{'EC:':<10} {EC}
{'alpha:':<10} {alpha}
{'states_total_number:':<10}{states_total_number}"""
        )

        voltage_matrix_dict = defaultdict(list)
        phi_matrix_dict = defaultdict(list)
        voltage_constant = (
            self.quantum_constants.h
            * 10 ** 9
            * EC
            / (2 * self.quantum_constants.eCharge * (1 + alpha))
        )

        for x in range(0, states_total_number):

            (state_indx, state_cp) = convert_index_to_island_state(x)

            voltage_elm = voltage_constant * (np.dot(state_cp, [1, 2, 1]))
            voltage_matrix_dict["row-col"].append(x)
            voltage_matrix_dict["elm"].append(voltage_elm)

            # 4 - phase operator e^{i phi_20}
            if state_indx[1] < (states_per_island - 1):
                y = convert_island_state_to_index(state_indx + [0, 1, 0])
                phi_matrix_dict["row"].append(x)
                phi_matrix_dict["col"].append(y)
                phi_matrix_dict["elm"].append(1)

        voltage_matrix = sp.coo_matrix(
            (
                voltage_matrix_dict["elm"],
                (voltage_matrix_dict["row-col"], voltage_matrix_dict["row-col"]),
            )
        ).tocsr()
        phi_matrix = sp.coo_matrix(
            (phi_matrix_dict["elm"], (phi_matrix_dict["row"], phi_matrix_dict["col"]))
        ).tocsr()

        logging.debug("🏭 Completed voltage and phase operator building")

        return voltage_matrix, phi_matrix
