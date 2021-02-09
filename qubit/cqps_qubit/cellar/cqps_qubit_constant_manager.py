"""

Assume that each JJ is a 100x100nm² square
"""

import logging
import sys

import numpy as np

import pinject

from common.terminal_colour import TerminalColour

pi = np.pi


class CqpsQubitConstantManager:
    # Proportional to resistance
    INDUCTANCE_PER_SQUARE = 1.5 * 10 ** (-9)

    @pinject.copy_args_to_public_fields
    def __init__(self, quantum_constants, generic_converter, param_dictionary):

        try:

            # Simulation parameters
            self.number_of_states = param_dictionary["number_of_states"]
            self.offset_to_apply = (self.number_of_states - 1) // 2

            self.ES = param_dictionary["ES"]

            # Physical parameters that can be used to derive the energies #########
            self.inductive_loop_squares = param_dictionary["inductive_loop_squares"]
            self.inductance = self.inductive_loop_squares * self.INDUCTANCE_PER_SQUARE
            self.EL = self.generic_converter.convert_energy_to_GHz(
                self.quantum_constants.phi0 ** 2 / (2 * self.inductance)
            )

            self.print_constants()

        except KeyError as err:
            logging.error(
                f"Need to pass in parameter ({err}) to the parameter dictionary"
            )
            sys.exit()

    def convert_index_to_state(self, index: int) -> int:
        if (index < 0) or (index >= self.number_of_states):
            raise ValueError(
                f"Index {index} must be within 0 and {self.number_of_states}"
            )

        return index - self.offset_to_apply

    def override_parameters(self, EL: float, ES: float):
        logging.info(
            f"""{TerminalColour.CRED}{TerminalColour.BOLD}Overriding to the following values:{TerminalColour.ENDC}
{'EL:':<10}{EL}
{'ES:':<10}{ES}
"""
        )

        self.ES = ES
        self.EL = EL

    def print_constants(self):
        logging.info(
            f"""⚛ System setup with the following parameters (can choose to override them):

{TerminalColour.CBLUEBG}{TerminalColour.CYELLOW}{TerminalColour.BOLD}Energies{TerminalColour.ENDC}
{'EL:':<50}{self.EL:.2f} (GHz)
{'ES:':<50}{self.ES:.2f} (GHz)

{TerminalColour.CBLUEBG}{TerminalColour.CYELLOW}{TerminalColour.BOLD}Raw Parameters{TerminalColour.ENDC}
{'number_of_states:':<50}{self.number_of_states}
{'inductive_loop_squares:':<50}{self.inductive_loop_squares:.2f} (100x100nm²)

{TerminalColour.CROMANBG}{TerminalColour.CWHITE}{TerminalColour.BOLD}Derived Parameters{TerminalColour.ENDC}
{'inductance:':<50}{self.inductance * 10**9:.2f} (nH)
"""
        )
