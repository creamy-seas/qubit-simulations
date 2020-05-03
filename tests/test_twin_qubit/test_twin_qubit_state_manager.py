from typing import List

import unittest
from unittest.mock import patch
from unittest.mock import Mock
from unittest.mock import MagicMock
from unittest.mock import call

import numpy as np

from qubit.twin_qubit.twin_qubit_state_manager import TwinQubitStateManager


class TestLogRequest(unittest.TestCase):
    def setUp(self):
        self.param_dictionary = {"states_per_island": 17}
        self.sut = TwinQubitStateManager(self.param_dictionary)

    def tearDown(self):
        pass

    def test_init(self,):
        self.assertEqual(self.sut.states_total_number, 4913)
        self.assertEqual(self.sut.cp_offset_to_apply, 8)

    def array_assertion(self, expected_array: np.array, actual_array: np.array) -> int:
        self.assertEqual(
            list(expected_array),
            list(actual_array),
            f"""
Expected array: {expected_array}
Actual array: {actual_array}
            """,
        )

    def test_convert_index_to_island_state(self):
        index_state, cp_state = self.sut.convert_index_to_island_state(0)
        self.array_assertion(np.array([0, 0, 0]), index_state)
        self.array_assertion(np.array([-8, -8, -8]), cp_state)

        index_state, cp_state = self.sut.convert_index_to_island_state(1)
        self.array_assertion(np.array([0, 0, 1]), index_state)
        self.array_assertion(np.array([-8, -8, -7]), cp_state)

        index_state, cp_state = self.sut.convert_index_to_island_state(5)
        self.array_assertion(np.array([0, 0, 5]), index_state)
        self.array_assertion(np.array([-8, -8, -3]), cp_state)

        index_state, cp_state = self.sut.convert_index_to_island_state(12)
        self.array_assertion(np.array([0, 0, 12]), index_state)
        self.array_assertion(np.array([-8, -8, 4]), cp_state)

        index_state, cp_state = self.sut.convert_index_to_island_state(17)
        self.array_assertion(np.array([0, 1, 0]), index_state)
        self.array_assertion(np.array([-8, -7, -8]), cp_state)

        index_state, cp_state = self.sut.convert_index_to_island_state(4912)
        self.array_assertion(np.array([16, 16, 16]), index_state)
        self.array_assertion(np.array([8, 8, 8]), cp_state)

        with self.assertRaises(ValueError):
            self.sut.convert_index_to_island_state(4913)

        with self.assertRaises(ValueError):
            self.sut.convert_index_to_island_state(-1)

    def test_convert_island_state_to_index(self):
        self.generic_convert_island_state_to_index([16, 16, 16], 4912)
        self.generic_convert_island_state_to_index([0, 0, 0], 0)
        self.generic_convert_island_state_to_index([0, 0, 1], 1)
        self.generic_convert_island_state_to_index([0, 0, 2], 2)
        self.generic_convert_island_state_to_index([0, 1, 0], 17)
        self.generic_convert_island_state_to_index([1, 0, 0], 289)

    def generic_convert_island_state_to_index(
        self, index_state: List, expected_index: int
    ):

        result = self.sut.convert_numeric_state_to_index(index_state)

        self.assertEqual(
            result,
            expected_index,
            f"""
Expected {index_state} -> {expected_index}
Got      {index_state} -> {result}""",
        )
