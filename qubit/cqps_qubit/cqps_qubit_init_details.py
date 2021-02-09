"""
Required class to pass in parameters directly using the pinject framework
"""

from typing import Dict
import pinject


class CqpsQubitInitDetails(pinject.BindingSpec):
    def __init__(
        self,
        param_dictionary: Dict,
        logging_level: int,
    ):
        self.param_dictionary = param_dictionary
        self.logging_level = logging_level

    def configure(self, bind):
        bind("logging_level", to_instance=self.logging_level)
        bind("param_dictionary", to_instance=self.param_dictionary)
