from typing import Tuple, List
from matplotlib.pyplot import Axes, Figure
import pinject


class TwinQubit:
    """Quantum class to work with the twin, 5JJ qubit"""

    @pinject.copy_args_to_internal_fields
    def __init__(
        self,
        quantum_logger,
        _twin_qubit_constant_manager,
        _twin_qubit_simulator,
        _twin_qubit_plotter,
    ):
        """
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        pass

    def override_parameters(self, EC=None, EJ=None, alpha=None, assymetry=None):
        self._twin_qubit_constant_manager.override_parameters(EC, EJ, alpha, assymetry)
        self._twin_qubit_constant_manager.print_constants()

    def run_simulation(self, evaluate_dipole_element=False):
        self._twin_qubit_simulator.simulate(evaluate_dipole_element)

    def plot_transitions(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_transitions(mpl_axes)

    def plot_dipole_moment_voltage(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_dipole_moment_voltage(mpl_axes)

    def plot_dipole_moment_voltage_beta(self, mpl_axes: Axes):
        self._twin_qubit_plotter.plot_dipole_moment_voltage_beta(mpl_axes)

    def generate_axes(
        self, nrows: int = 1, ncols: int = 1
    ) -> Tuple[Figure, List[Axes]]:
        return self._twin_qubit_plotter.generate_axes(nrows, ncols)
