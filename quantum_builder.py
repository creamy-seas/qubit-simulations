"""
Imports all of the classes that resolve all the dependencies of the project
"""

from typing import Dict, List

import pinject

# CQPS qubit ##################################################################
from qubit.cqps_qubit.cqps_qubit import CqpsQubit
from qubit.cqps_qubit.cellar.cqps_qubit_constant_manager import CqpsQubitConstantManager
from qubit.cqps_qubit.simulation.cqps_qubit_hamiltonian_manager import (
    CqpsQubitHamiltonianManager,
)
from qubit.cqps_qubit.simulation.cqps_qubit_simulator import CqpsQubitSimulator

# Twin Qubit ##################################################################
from qubit.twin_qubit.twin_qubit import TwinQubit
from qubit.twin_qubit.cellar.twin_qubit_constant_manager import TwinQubitConstantManager
from qubit.twin_qubit.cellar.twin_qubit_state_manager import TwinQubitStateManager
from qubit.twin_qubit.simulation.twin_qubit_hamiltonian_manager import (
    TwinQubitHamiltonianManager,
)
from qubit.twin_qubit.simulation.twin_qubit_operator_builder import (
    TwinQubitOperatorBuilder,
)
from qubit.twin_qubit.simulation.twin_qubit_simulator import TwinQubitSimulator
from qubit.twin_qubit.simulation.twin_qubit_local_fluctuation_simulator import (
    TwinQubitLocalFluctuationSimulator,
)
from qubit.twin_qubit.plotting.twin_qubit_plotter import TwinQubitPlotter
from qubit.twin_qubit.plotting.twin_qubit_sparse_matrix_visualiser import (
    TwinQubitSparseMatrixVisualiser,
)
from qubit.twin_qubit.simulation.twin_qubit_simulator_phil_phir import (
    TwinQubitSimulatorPhilPhir,
)

# Transmon Qubit ##############################################################
from qubit.transmon_qubit.transmon_qubit import TransmonQubit
from qubit.transmon_qubit.cellar.transmon_qubit_constant_manager import (
    TransmonQubitConstantManager,
)
from qubit.transmon_qubit.simulation.transmon_qubit_hamiltonian_manager import (
    TransmonQubitHamiltonianManager,
)
from qubit.transmon_qubit.simulation.transmon_qubit_simulator import (
    TransmonQubitSimulator,
)

# Common ######################################################################
from qubit.utils.quantum_constants import QuantumConstants
from qubit.utils.quantum_logger import QuantumLogger
from qubit.utils.generic_converter import GenericConverter
from qubit.utils.init_details import InitDetails


class QuantumBuilder:
    @classmethod
    def build_twin_qubit(
        cls,
        param_dictionary: Dict,
        flux_list: List,
        logging_level: int,
        other_parameters: Dict = {"empty": "empty"},
    ):

        # Load the details into the InitDetails class
        BINDING_SPECS = [
            InitDetails(
                param_dictionary=param_dictionary,
                logging_level=logging_level,
                other_parameters=other_parameters,
                flux_list=flux_list,
            )
        ]

        # Construct the graph, binding fields of the InitDetails class
        OBJ_GRAPH = pinject.new_object_graph(binding_specs=BINDING_SPECS)

        # Perform the build
        TWIN_QUBIT = OBJ_GRAPH.provide(TwinQubit)

        return TWIN_QUBIT

    @classmethod
    def build_transmon_qubit(
        cls,
        param_dictionary: Dict,
        logging_level: int,
    ):

        # Load the details into the InitDetails class
        BINDING_SPECS = [
            InitDetails(
                param_dictionary=param_dictionary,
                logging_level=logging_level,
                other_parameters=-1,
            )
        ]

        # Construct the graph, binding fields of the InitDetails class
        OBJ_GRAPH = pinject.new_object_graph(binding_specs=BINDING_SPECS)

        # Perform the build
        return OBJ_GRAPH.provide(TransmonQubit)

    @classmethod
    def build_cqps_qubit(
        cls,
        param_dictionary: Dict,
        logging_level: int,
    ):

        # Load the details into the InitDetails class
        BINDING_SPECS = [
            InitDetails(param_dictionary, logging_level, other_parameters=-1)
        ]

        # Construct the graph, binding fields of the InitDetails class
        OBJ_GRAPH = pinject.new_object_graph(binding_specs=BINDING_SPECS)

        # Perform the build
        return OBJ_GRAPH.provide(CqpsQubit)
