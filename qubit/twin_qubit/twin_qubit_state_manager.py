import logging
from typing import Tuple, List
import numpy as np


class TwinQubitStateManager:
    """Class that converts between indexed and island representation of a quantum state"""

    def __init__(self, param_dictionary):
        self.states_per_island = param_dictionary["states_per_island"]
        self.states_total_number = self.states_per_island ** 3
        self.cp_offset_to_apply = (self.states_per_island - 1) // 2

        self.verify_simulation_parameters()

    def verify_simulation_parameters(self):
        """Correct common errors before program begins"""
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
        across the 3 islands of the system - (n1, n2, n3), begging with the lowest number

        __ Returns __
        [int, int, int] state_numeric:            3-island index distribution
        [int, int, int] state_cp:              3-island cp distribution, centered about (0,0,0)
        """

        if (index < 0) or (index >= self.states_total_number):
            raise ValueError(
                f"Index {index} must be within 0 and {self.states_total_number} -> converting to 0"
            )

        # Convert index to  represetnation on the 3 islands
        # Order of bits is reversed: 011 -> 110 otherwise 0's will be trimmed
        state_numeric = np.zeros(3, dtype=int)
        for idx_i, i in enumerate(
            np.base_repr(index, base=self.states_per_island)[::-1]
        ):
            state_numeric[idx_i] = int(i, base=self.states_per_island)

        # Reverse the order: [5, 0, 0] -> [0, 0, 5] as we want the array to increase from rightmost entry
        state_numeric = state_numeric[::-1]

        state_cp = state_numeric - self.cp_offset_to_apply

        return state_numeric, state_cp

    def convert_numeric_state_to_index(self, state_numeric: List) -> int:
        """
        __ Parameters __
        [int, int, int] -> inique index from 0 to states_per_island^3
        """
        index = 0
        for (power, value) in enumerate(state_numeric[::-1]):
            index += value * self.states_per_island ** power

        return index
