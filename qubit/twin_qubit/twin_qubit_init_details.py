from typing import Dict, List
import pinject


class TwinQubitInitDetails(pinject.BindingSpec):
    def __init__(
        self,
        param_dictionary: Dict,
        flux_list: List,
        logging_level: int,
        other_parameters: Dict,
    ):
        self.param_dictionary = param_dictionary
        self.logging_level = logging_level
        self.flux_list = flux_list
        self.other_parameters = other_parameters

    def configure(self, bind):
        bind("param_dictionary", to_instance=self.param_dictionary)
        bind("logging_level", to_instance=self.logging_level)
        bind("flux_list", to_instance=self.flux_list)
        bind("other_parameters", to_instance=self.other_parameters)
