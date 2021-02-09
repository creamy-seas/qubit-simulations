"""
CQPS qubit
"""

from typing import Tuple, List, Dict
from matplotlib.pyplot import Axes, Figure
import pinject


class CqpsQubit:
    @pinject.copy_args_to_public_fields
    def __init__(
        self,
        quantum_logger,
        cqps_qubit_constant_manager,
        cqps_qubit_init_details,
        cqps_qubit_simulator,
        cqps_qubit_hamiltonian_manager,
    ):
        """
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        pass

    def override_parameters(self, EL=None, ES=None):
        self.cqps_qubit_constant_manager.override_parameters(EL, ES)

    def run_simulation(
        self,
        flux_ext_list: List[float],
        number_of_levels_to_simulate: int,
        use_sparse_matrix=True,
    ) -> Dict:
        """
        @brief      Run simulation of CQPS energies and eigenstates for each of the supplied fluxes

        @param      flux_ext_list: List of fluxes
                    use_sparse_matrix: True to use faster approximation

        @return     None
        """
        return self.cqps_qubit_simulator.simulate(
            flux_ext_list,
            number_of_levels_to_simulate,
            use_sparse_matrix,
        )
