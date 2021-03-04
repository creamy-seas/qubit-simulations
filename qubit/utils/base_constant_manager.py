"""
Base class with useful functions and parmeters, that qubit classes can derive and use. Maybe -- still in workings.
"""
from abc import ABC
from typing import Dict
import logging
import sys

import pinject


class BaseConstantManager(ABC):
    @pinject.copy_args_to_public_fields
    def __init__(self, quantum_constants, generic_converter, param_dictionary: Dict):
        """
        Param dictionary should supply:
        - 'number_of_states'
        """
        try:
            # Simulation parameters
            self.number_of_states = param_dictionary["number_of_states"]
            self.offset_to_apply = (self.number_of_states - 1) // 2

        except KeyError as err:
            logging.error(
                f"Need to pass in parameter ({err}) to the parameter dictionary"
            )
            sys.exit()

    def convert_index_to_state(self, index: int) -> int:
        """
        Index used during iteration (0, 1, 2 ...) is converted to a state of the system
        """
        if (index < 0) or (index >= self.number_of_states):
            raise ValueError(
                f"Index {index} must be within 0 and {self.number_of_states}"
            )

        return index - self.offset_to_apply

    @ABC.abstractmethod
    def print_constants(self):
        """Method to display the derived settings"""
        pass
