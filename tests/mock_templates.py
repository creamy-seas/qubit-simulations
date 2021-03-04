import numpy as np


def array_assertion(self, expected_array: np.array, actual_array: np.array) -> int:
    self.assertEqual(
        list(expected_array),
        list(actual_array),
        f"""
Expected array: {expected_array}
Actual array: {actual_array}
            """,
    )
