from typing import Dict, List
import pinject


class TwinQubitInitDetails(pinject.BindingSpec):
    def __init__(self, param_dictionary: Dict, flux_list: List, logging_level: int):
        self.param_dictionary = param_dictionary
        self.logging_level = logging_level
        self.flux_list = flux_list

    def configure(self, bind):
        bind("param_dictionary", to_instance=self.param_dictionary)
        bind("logging_level", to_instance=self.logging_level)
        bind("flux_list", to_instance=self.flux_list)
