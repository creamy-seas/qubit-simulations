from typing import Tuple, List
from matplotlib.pyplot import Axes, Figure
import pinject


class TransmonQubit:
    """Quantum class to work with the twin, 5JJ qubit"""

    @pinject.copy_args_to_public_fields
    def __init__(
        self,
        quantum_logger,
        transmon_qubit_constant_manager,
        transmon_qubit_simulator,
        init_details,
        transmon_qubit_hamiltonian_manager,
    ):
        """
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        pass
