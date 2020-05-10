"""
Imports all of the classes that resolve all the dependencies of the project
"""

from typing import Dict, List

import pinject

from qubit.twin_qubit.twin_qubit import TwinQubit
from qubit.twin_qubit.twin_qubit_init_details import TwinQubitInitDetails
from qubit.twin_qubit.cellar.twin_qubit_constant_manager import TwinQubitConstantManager
from qubit.twin_qubit.cellar.twin_qubit_state_manager import TwinQubitStateManager
from qubit.twin_qubit.simulation.twin_qubit_hamiltonian_manager import (
    TwinQubitHamiltonianManager,
)
from qubit.twin_qubit.simulation.twin_qubit_operator_builder import (
    TwinQubitOperatorBuilder,
)
from qubit.twin_qubit.simulation.twin_qubit_simulator import TwinQubitSimulator
from qubit.twin_qubit.plotting.twin_qubit_plotter import TwinQubitPlotter
from qubit.twin_qubit.plotting.twin_qubit_sparse_matrix_visualiser import (
    TwinQubitSparseMatrixVisualiser,
)

from qubit.utils.quantum_constants import QuantumConstants
from qubit.utils.quantum_logger import QuantumLogger
from qubit.utils.generic_converter import GenericConverter


class QuantumBuilder:
    @classmethod
    def build_twin_qubit(
        cls, param_dictionary: Dict, flux_list: List, logging_level: int
    ):

        # Load the details into the InitDetails class
        BINDING_SPECS = [
            TwinQubitInitDetails(param_dictionary, flux_list, logging_level)
        ]

        # Construct the graph, binding fields of the InitDetails class
        OBJ_GRAPH = pinject.new_object_graph(binding_specs=BINDING_SPECS)

        # Perform the build
        TWIN_QUBIT = OBJ_GRAPH.provide(TwinQubit)

        return TWIN_QUBIT
