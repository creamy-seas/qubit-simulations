import logging
from collections import defaultdict
from typing import List, Tuple

from pyprind import ProgBar
import numpy as np
from scipy.sparse.linalg import eigsh
import pinject


class TwinQubitLocalFluctuationSimulator:
    CONDITION_SET = {"ANTICORRELATED", "ABSENT", "RANDOM"}

    @pinject.copy_args_to_public_fields
    def __init__(self, twin_qubit_simulator, twin_qubit_hamiltonian_manager):
        self.simulations = defaultdict(list)

    def generate_fluctuations(
        self, mu: float, sigma: float, number_simulations: int, condition: str
    ) -> Tuple[List, List]:

        assert (
            condition in self.CONDITION_SET
        ), f"'condition' must be one of {self.CONDITION_SET}"

        # Generate distribution
        fluctuation_list_1 = np.random.normal(mu, sigma, number_simulations)
        if condition == "ANTICORRELATED":
            fluctuation_list_2 = mu - (fluctuation_list_1 - mu)
        elif condition == "RANDOM":
            fluctuation_list_2 = np.random.normal(mu, sigma, number_simulations)
        else:
            fluctuation_list_1 = np.array(number_simulations * [0.5])
            fluctuation_list_2 = np.array(number_simulations * [0.5])

        logging.warning(
            f"""🐁 Generated fluctuations:
{'μ:':<30}{mu}
{'σ:':<30}{sigma}
{'number of simulations:':<30}{number_simulations}
{'condition:':<30}{condition}
{'Example pairs:':<30}{(fluctuation_list_1[-1], fluctuation_list_2[-1])}
{'':<30}{(fluctuation_list_1[-2], fluctuation_list_2[-2])}
"""
        )

        # Apply conversion
        fluctuation_list_1 *= 2 * np.pi
        fluctuation_list_2 *= 2 * np.pi

        return (fluctuation_list_1, fluctuation_list_2)

    def run_simulation(
        self, mu: float, sigma: float, number_simulations: int, condition: str
    ):

        (fluctuation_list_1, fluctuation_list_2) = self.generate_fluctuations(
            mu, sigma, number_simulations, condition
        )
        self.twin_qubit_hamiltonian_manager.stage2_prepare_constant_hamiltonian()
        self.simulations = defaultdict(list)

        progress_bar = ProgBar(number_simulations, bar_char="█")
        for (fluctuation_1, fluctuation_2) in zip(
            fluctuation_list_1, fluctuation_list_2
        ):

            self.twin_qubit_hamiltonian_manager.stage3_build_hamiltonian_for_simulation(
                fluctuation_1, fluctuation_2
            )

            (eigvals, eigvecs) = eigsh(
                self.twin_qubit_hamiltonian_manager.hamiltonian_simulation,
                3,
                which="SA",
                tol=0,
            )

            (
                eigvals,
                eigvecs,
            ) = self.twin_qubit_simulator.sort_in_ascending_eigval_order(
                eigvals, eigvecs
            )

            self.simulations["1-2-with-local-fluctuations"].append(
                eigvals[1] - eigvals[0]
            )

            progress_bar.update()
