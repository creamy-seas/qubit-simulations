import unittest
from unittest.mock import Mock
from unittest.mock import patch

import scipy.sparse as sp
import numpy as np

from qubit.cqps_twin_qubit.cqps_twin_qubit_hamiltonian_manager import (
    CqpsTwinQubitHamiltonianManager,
)
from qubit.cqps_twin_qubit.cqps_twin_qubit_constant_manager import (
    CqpsTwinQubitConstantManager,
)
from qubit.cqps_twin_qubit.cqps_twin_qubit_simulator import CqpsTwinQubitSimulator
from qubit.utils.generic_converter import GenericConverter
from qubit.utils.quantum_constants import QuantumConstants


class TestCqpsTwinQubitSimulator(unittest.TestCase):
    def test_unpack_parameters__balanced(self):
        (phi_l_list, phi_r_list,) = CqpsTwinQubitSimulator.unpack_phi_sweep_dictionary(
            {"type": "balanced", "trash_1": 2, "phi_list": [1, 2, 3]}
        )

        self.assertEqual(phi_l_list, [1, 2, 3])
        self.assertEqual(phi_r_list, [None])

    def test_unpack_parameters__phi_l_phi_r(self):
        (phi_l_list, phi_r_list,) = CqpsTwinQubitSimulator.unpack_phi_sweep_dictionary(
            {
                "type": "phi_l-phi_r",
                "trash_1": 2,
                "phi_r_list": [3, 3],
                "phi_l_list": [1],
            }
        )

        self.assertEqual(phi_l_list, [1])
        self.assertEqual(phi_r_list, [3, 3])

    def test_unpack_parameters__phi_plus_phi_minus(self):
        (phi_l_list, phi_r_list,) = CqpsTwinQubitSimulator.unpack_phi_sweep_dictionary(
            {
                "type": "phi_+-phi_-",
                "phi_+_list": [3, 3],
                "phi_-_list": [1, 2],
            }
        )

        self.assertEqual(phi_l_list, [4, 5])
        self.assertEqual(phi_r_list, [2, 1])

    def test_sort_in_ascending_order(self):

        # We simulate 4 energy levels (4 eigenvalues)
        # The sytem uses 2 states
        eigvals = np.array([3, 2, 1, 4])
        # Note that eigenvectors will need to be supplied in this transformed way i.e. amplitude of asingle vector at all it's eigevalues
        eigvecs = np.array(
            [
                [1, 2, 3, 4],
                [10, 20, 30, 40],
            ]
        )

        (eigvals, eigvecs) = CqpsTwinQubitSimulator.sort_in_ascending_eigval_order(
            eigvals, eigvecs
        )

        expected_eigvals = np.array([1, 2, 3, 4])
        expected_eigvecs = np.array(
            [
                [3, 30],
                [2, 20],
                [1, 10],
                [4, 40],
            ]
        )

        # Iterate eigenvalues
        for i in range(0, 4):
            self.assertEqual(eigvals[i], expected_eigvals[i])
            # For each eigenvalue, check eigenvector
            for j in range(0, 2):
                self.assertEqual(eigvecs[i][j], expected_eigvecs[i][j])

    def test_store_result(self):
        eigvals = np.array([1, 2, 3, 4])
        eigvecs = np.array(
            [
                [3, 30],
                [2, 20],
                [1, 10],
                [4, 40],
            ]
        )
        NO_LEVELS_SIMULATED = 4
        NO_STATES_USED = 2
        DIM_L = 2
        DIM_R = 3

        simulation_dictionary = {
            "eigvals": np.empty((DIM_L, DIM_R, NO_LEVELS_SIMULATED), dtype=float),
            "eigvecs": np.empty(
                (
                    DIM_L,
                    DIM_R,
                    NO_LEVELS_SIMULATED,
                    NO_STATES_USED,
                ),
                dtype=np.cdouble,
            ),
        }

        test_result = CqpsTwinQubitSimulator.store_results(
            simulation_dictionary, eigvals, eigvecs, phi_l_idx=1, phi_r_idx=0
        )

        for i in range(0, NO_LEVELS_SIMULATED):
            self.assertEqual(
                test_result["eigvals"][1][0][i],
                simulation_dictionary["eigvals"][1][0][i],
            )

            # need to test eignenvecs
