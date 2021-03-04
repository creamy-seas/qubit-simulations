from typing import List, Tuple
import logging
import numpy as np

import pinject

from common.terminal_colour import TerminalColour


class CqpsTwinQubitConstantManager:
    # Proportional to resistance
    INDUCTANCE_PER_SQUARE = 1.5 * 10 ** (-9)

    @pinject.copy_args_to_public_fields
    def __init__(self, quantum_constants, generic_converter, param_dictionary):

        try:
            # Simulation parameters
            self.states_per_loop = param_dictionary["states_per_loop"]
            self.states_total_number = self.states_per_loop ** 2
            self.offset_to_apply = (self.states_per_loop - 1) // 2

            # Phase slip "amplitude" across the central and side constrictions
            self.ES = param_dictionary["ES"]
            self.ES_on_sides = param_dictionary["ES_on_sides"]

            # Energies associated with the branches of the inductive loops
            self.inductive_loop_squares_left = param_dictionary[
                "inductive_loop_squares_left"
            ]
            self.inductive_loop_squares_right = param_dictionary[
                "inductive_loop_squares_right"
            ]
            self.EL_left = self.generic_converter.convert_energy_to_GHz(
                self.quantum_constants.phi0 ** 2
                / (2 * self.inductive_loop_squares_left * self.INDUCTANCE_PER_SQUARE)
            )
            self.EL_right = self.generic_converter.convert_energy_to_GHz(
                self.quantum_constants.phi0 ** 2
                / (2 * self.inductive_loop_squares_right * self.INDUCTANCE_PER_SQUARE)
            )

            self.print_constants()

        except KeyError as err:
            logging.error(
                "Need to pass in parameter ({err}) to the parameter dictionary".format(
                    err=err
                )
            )
            raise RuntimeError(
                "Need to pass in parameter ({err}) to the parameter dictionary".format(
                    err=err
                )
            )

    def convert_index_to_state(self, index: int) -> Tuple[np.array, np.array]:
        """
        __ Parameters __
        [int] index:                    to convert into a one-to-one system state
                                        (0 <= index < states_total_number)

        __ Description __
        Takes an index and converts it to a unique flux state of the
        system (f1, f2), beggining with the lowest number

        __ Returns __
        [int, int] state_numeric:           2-loop index distribution
        [int, int] state_flux:              2-loop flux distribution, centered about (0,0)

        [0,0] -> [1,0] -> [2,0] -> [0,1]
        """
        if (index < 0) or (index >= self.states_total_number):
            raise ValueError(
                f"Index {index} must be within 0 and {self.states_total_number} -> converting to 0"
            )

        # Order of bits is reversed: 011 -> 110 otherwise 0's will be trimmed
        state_numeric = np.zeros(2, dtype=int)
        for (idx_i, i) in enumerate(
            np.base_repr(index, base=self.states_per_loop)[::-1]
        ):
            state_numeric[idx_i] = int(i, base=self.states_per_loop)

        state_flux = state_numeric - self.offset_to_apply

        return state_numeric, state_flux

    def convert_state_to_index(self, state: List[int]) -> int:
        """
        __ Parameters __
        [int, int] -> inique index from 0 to states_per_loop^2
        """
        index = 0
        for (power, value) in enumerate(state):
            index += value * self.states_per_loop ** power
        return index

    def override_parameters(
        self,
        EL_left: float = None,
        EL_right: float = None,
        ES: float = None,
        ES_on_sides: float = None,
    ):
        logging.info(
            f"{TerminalColour.CRED}{TerminalColour.BOLD}Overriding to the values{TerminalColour.ENDC}"
        )

        if ES_on_sides:
            self.ES_on_sides = ES_on_sides
        if ES:
            self.ES = ES
        if EL_right:
            self.EL_right = EL_right
        if EL_left:
            self.EL_left = EL_left

        self.print_constants()

    def print_constants(self):
        logging.info(
            f"""⚛ System setup with the following parameters (can choose to override them):

{TerminalColour.CBLUEBG}{TerminalColour.CYELLOW}{TerminalColour.BOLD}Energies{TerminalColour.ENDC}
{'EL_left:':<50}{self.EL_left:.2f} (GHz)
{'EL_right:':<50}{self.EL_right:.2f} (GHz)
{'ES:':<50}{self.ES:.2f} (GHz)
{'ES_on_sides:':<50}{self.ES_on_sides:.2f} (GHz)

{TerminalColour.CBLUEBG}{TerminalColour.CYELLOW}{TerminalColour.BOLD}Raw Parameters{TerminalColour.ENDC}
{'states_per_loop:':<50}{self.states_per_loop}
{'inductive_loop_squares_left:':<50}{self.inductive_loop_squares_left:.2f} (100x100nm²)
{'inductive_loop_squares_right:':<50}{self.inductive_loop_squares_right:.2f} (100x100nm²)
"""
        )
