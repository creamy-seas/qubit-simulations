from typing import Dict, List
import pinject


class InitDetails(pinject.BindingSpec):
    def __init__(
        self,
        param_dictionary: Dict,
        logging_level: int,
        other_parameters: Dict,
        flux_list: List = [],
    ):
        self.param_dictionary = param_dictionary
        self.logging_level = logging_level
        self.other_parameters = other_parameters
        self.flux_list = flux_list

    def configure(self, bind):
        bind("param_dictionary", to_instance=self.param_dictionary)
        bind("logging_level", to_instance=self.logging_level)
        bind("other_parameters", to_instance=self.other_parameters)
        bind("flux_list", to_instance=self.flux_list)
