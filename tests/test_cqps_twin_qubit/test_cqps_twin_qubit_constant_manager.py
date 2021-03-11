import unittest
from unittest.mock import Mock
from unittest.mock import patch

import numpy as np
import scipy.sparse as sp

from tests.mock_templates import array_assertion
from qubit.cqps_twin_qubit.cqps_twin_qubit_constant_manager import (
    CqpsTwinQubitConstantManager,
)
from qubit.qutils.generic_converter import GenericConverter
from qubit.qutils.quantum_constants import QuantumConstants


class TestCqpsTwinQubitConstantManager(unittest.TestCase):
    def setUp(self):
        test_params = {
            "states_per_loop": 3,
            "ES": 1,
            "ES_on_sides": 2,
            "inductive_loop_squares_left": 5,
            "inductive_loop_squares_right": 4,
        }
        self.sut = CqpsTwinQubitConstantManager(
            quantum_constants=QuantumConstants(),
            generic_converter=GenericConverter(QuantumConstants()),
            param_dictionary=test_params,
        )

        self.sut.override_parameters(ES=2)
        self.assertEqual(self.sut.ES, 2)

    def test_bad_params(self):
        with self.assertRaises(RuntimeError):
            CqpsTwinQubitConstantManager(
                quantum_constants=QuantumConstants(),
                generic_converter=GenericConverter(QuantumConstants()),
                param_dictionary={
                    "states_per_loop": 3,
                    "ES": 1,
                    "ES_on_sides": 2,
                    "inductive_loop_squares_left": 5,
                },
            )

    def test_convert_index_to_state(self):

        # Out of bounds
        with self.assertRaises(ValueError):
            self.sut.convert_index_to_state(9)

        index_state, loop_state = self.sut.convert_index_to_state(0)
        array_assertion(self, np.array([0, 0]), index_state)
        array_assertion(self, np.array([-1, -1]), loop_state)

        index_state, loop_state = self.sut.convert_index_to_state(6)
        array_assertion(self, np.array([0, 2]), index_state)
        array_assertion(self, np.array([-1, 1]), loop_state)

    def test_convert_numeric_state_to_index(self):
        self.assertEqual(self.sut.convert_state_to_index([0, 1]), 3)
