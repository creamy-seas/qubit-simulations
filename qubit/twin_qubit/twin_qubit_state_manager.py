import logging
from typing import Tuple, List

import numpy as np


class TwinQubitStateManager:
    """Class that converts between indexed and island representation of a quantum state"""

    def __init__(self, param_dictionary):
        self.states_per_island = param_dictionary["states_per_island"]
        self.states_total_number = self.states_per_island ** 3
        self.cp_offset_to_apply = (self.states_per_island - 1) // 2

        self.print_simulation_parameters()

    def print_simulation_parameters(self):
        """Correct common errors before program begins"""
        if (self.states_total_number % 2) == 0:
            raise ValueError(
                f"Total number of states is {self.states_total_number} -> it needs to be odd"
            )
        logging.info(
            f"""⚙ Quantum state manager is using:
{self.states_per_island:<5} states per island
{self.states_total_number:<5} total states
"""
        )

    def convert_index_to_island_state(self, index: int) -> Tuple[np.array, np.array]:
        """
        __ Parameters __
        [int] index:                    to convert into a one-to-one system state
                                        (0 <= index < states_total_number)

        __ Description __
        Takes an index and converts it to a unique cooper pair distribution
        across the 3 islands of the system - (n1, n2, n3).

        __ Returns __
        [int, int, int] state_indx:            3-island index distribution
        [int, int, int] state_cp:              3-island cp distribution, centered about (0,0,0)
        """

        if (index < 0) or (index >= self.states_total_number):
            raise ValueError(
                "Index {index} must be within 0 and {self.states_total_number} -> converting to 0"
            )

        # 1 - convert index to  represetnation on the 3 islands
        state_indx = np.zeros(3, dtype=int)
        for idx_i, i in enumerate(
            np.base_repr(index, base=self.states_per_island)[::-1]
        ):
            state_indx[idx_i] = int(i)

        # 2 - generate cp distribution, centered around (0, 0, 0)
        state_cp = state_indx - self.cp_offset_to_apply

        return state_indx, state_cp

    def convert_island_state_to_index(self, state_indx: List) -> int:
        """
        __ Parameters __
        [int, int, int] -> inique index from 0 to states_per_island^3
        """

        return int("".join(map(str, state_indx)), base=self.states_per_island)
