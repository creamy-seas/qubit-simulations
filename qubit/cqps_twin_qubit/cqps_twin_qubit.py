"""
CQPS twin qubit
"""

from typing import List, Dict
import pinject


class CqpsTwinQubit:
    @pinject.copy_args_to_public_fields
    def __init__(
        self,
        init_details,
        quantum_logger,
        cqps_twin_qubit_constant_manager,
        cqps_twin_qubit_simulator,
        cqps_twin_qubit_hamiltonian_manager,
    ):
        """
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        pass

    def override_parameters(
        self,
        EL_right: float = None,
        EL_left: float = None,
        ES: float = None,
        ES_on_sides: float = None,
    ):
        self.cqps_twin_qubit_constant_manager.override_parameters(
            EL_left=EL_left,
            EL_right=EL_right,
            ES=ES,
            ES_on_sides=ES_on_sides,
        )

    def run_simulation(
        self,
        phi_dict: Dict,
        number_of_levels_to_simulate: int,
        use_sparse_matrix=True,
    ) -> Dict:
        """
        @brief      Run simulation of CQPS Twin energies and eigenstates for each of the supplied fluxes

        @param      phi_dict: Details on how to vary flux in simulation
                    use_sparse_matrix: True to use faster approximation

        @return     None
        """
        return self.cqps_twin_qubit_simulator.simulate(
            phi_dict=phi_dict,
            number_of_levels_to_simulate=number_of_levels_to_simulate,
            use_sparse_matrix=use_sparse_matrix,
        )
