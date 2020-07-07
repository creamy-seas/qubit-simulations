from typing import Tuple, List
from matplotlib.pyplot import Axes, Figure
import pinject


class TwinQubit:
    """Quantum class to work with the twin, 5JJ qubit"""

    @pinject.copy_args_to_internal_fields
    def __init__(
        self,
        quantum_logger,
        twin_qubit_constant_manager,
        twin_qubit_simulator,
        twin_qubit_plotter,
        twin_qubit_sparse_matrix_visualiser,
        twin_qubit_init_details,
        twin_qubit_local_fluctuation_simulator,
    ):
        """
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        self.simulations = None
        self.flux_list = twin_qubit_init_details.flux_list

    def override_parameters(self, EC=None, EJ=None, alpha=None, assymetry=None):
        self._twin_qubit_constant_manager.override_parameters(EC, EJ, alpha, assymetry)
        self._twin_qubit_constant_manager.print_constants()

    def run_simulation(self, evaluate_dipole_element=False):
        self._twin_qubit_simulator.simulate(evaluate_dipole_element)
        self.simulations = self._twin_qubit_simulator.simulations

    def run_fluctuation_simulations(
        self, mu: float, sigma: float, number_simulations: int, condition: str
    ):
        self._twin_qubit_local_fluctuation_simulator.run_simulation(
            mu, sigma, number_simulations, condition
        )
        self.fluctuation_simulations = (
            self._twin_qubit_local_fluctuation_simulator.simulations
        )

    def plot_transitions(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_transitions(mpl_axes)

    def plot_dipole_matrix_elements(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_dipole_matrix_elements(mpl_axes)

    def plot_dipole_moment_voltage(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_dipole_moment_voltage(mpl_axes)

    def plot_dipole_moment_voltage_beta(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_dipole_moment_voltage_beta(mpl_axes)

    def plot_sparse_matrix(self, mpl_axes: Axes):
        self._twin_qubit_sparse_matrix_visualiser.visualise_matrix(mpl_axes)
