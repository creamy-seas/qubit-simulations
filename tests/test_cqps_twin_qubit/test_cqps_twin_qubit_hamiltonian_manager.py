import unittest
from unittest.mock import Mock
from unittest.mock import patch

import scipy.sparse as sp

from qubit.cqps_twin_qubit.cqps_twin_qubit_hamiltonian_manager import (
    CqpsTwinQubitHamiltonianManager,
)
from qubit.cqps_twin_qubit.cqps_twin_qubit_constant_manager import (
    CqpsTwinQubitConstantManager,
)
from qubit.qutils.generic_converter import GenericConverter
from qubit.qutils.quantum_constants import QuantumConstants


class TestCqpsTwinQubitHamiltonianManager(unittest.TestCase):
    def test(self):
        STATES_PER_LOOP = 3
        test_params = {
            "states_per_loop": STATES_PER_LOOP,
            "ES": 0,
            "ES_on_sides": 0,
            "inductive_loop_squares_left": 1,
            "inductive_loop_squares_right": 1,
        }
        mock_constant_manager = CqpsTwinQubitConstantManager(
            quantum_constants=QuantumConstants(),
            generic_converter=GenericConverter(QuantumConstants()),
            param_dictionary=test_params,
        )

        # Init the class and build Hamiltonian in stages
        EL_left = 2
        EL_right = 3
        ES = 5
        ES_on_sides = 7
        phi_l = 2
        phi_r = 1
        mock_constant_manager.override_parameters(
            EL_left=EL_left, EL_right=EL_right, ES=ES, ES_on_sides=ES_on_sides
        )
        sut = CqpsTwinQubitHamiltonianManager(mock_constant_manager)
        test_result = sut.stage_2_build_hamiltonian_for_simulation(
            phi_l=phi_l, phi_r=phi_r
        )

        inductance_val = [30, 20, 14, 21, 11, 5, 18, 8, 2]
        inductance_row_col = [0, 1, 2, 3, 4, 5, 6, 7, 8]

        loop_exchange_ES_val = [-2.5] * 8
        loop_exchange_ES_row = [1, 3, 2, 4, 6, 4, 7, 5]
        loop_exchange_ES_col = [3, 1, 4, 2, 4, 6, 5, 7]

        env_escape_val = [-3.5] * 24
        env_escape_col = [
            0,
            0,
            1,
            1,
            1,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            5,
            5,
            5,
            6,
            6,
            7,
            7,
            7,
            8,
            8,
        ]
        env_escape_row = [
            1,
            3,
            0,
            2,
            4,
            1,
            5,
            0,
            4,
            6,
            1,
            3,
            5,
            7,
            2,
            4,
            8,
            3,
            7,
            4,
            6,
            8,
            5,
            7,
        ]

        expected_result = sp.coo_matrix(
            (
                inductance_val + loop_exchange_ES_val + env_escape_val,
                (
                    inductance_row_col + loop_exchange_ES_row + env_escape_row,
                    inductance_row_col + loop_exchange_ES_col + env_escape_col,
                ),
            )
        ).tocsr()

        for i in range(0, 9):
            for j in range(0, 9):
                self.assertEqual(
                    test_result[i, j],
                    expected_result[i, j],
                    f"Failure in index {i},{j}. Expected {expected_result[i, j]} but got {test_result[i, j]}",
                )
