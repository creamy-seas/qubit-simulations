"""
Hamiltonian for twin cqps system. Eigenstates enumerate the fluxes in the left and RIGHT loops: |nl, nr>

This will first produce a skeleton, that can will be scaled by the energies RIGHT before scaling by energies
"""
from collections import defaultdict, OrderedDict
from typing import Tuple
import logging

import pinject
import numpy as np
import scipy.sparse as sp

# Defined, so that the first element in numpy array is for the left
# [LEFT, RIGHT]
LEFT = 0
RIGHT = 1


class CqpsTwinQubitHamiltonianManager:
    @pinject.copy_args_to_public_fields
    def __init__(self, cqps_twin_qubit_constant_manager):
        pass

    def stage1_prepare_hamiltonian_skeleton(self):
        """
        __ Description __
        Generates Hamiltonian that will be scaled in stage2 and stage3 when bias in the left
        and RIGHT loops, (phi_l, phi_r) are applied.

        |---------------------+-----------------------------------------|
        |                     | Scaling                                 |
        |---------------------+-----------------------------------------|
        | diagonal            | (nl - phi_l)^2/2Ll + (fr - phi_r)^2/2Ll |
        | loop-loop-tunneling | - ES_center/2                                  |
        | loop-env-tunneling  | - ES_right/2                         |
        |---------------------+-----------------------------------------|
        """
        # Reset parameters for simulation
        self.hamiltonian_elements = OrderedDict()
        self.hamiltonian_elements["inductance"] = defaultdict(list)
        self.hamiltonian_elements["loop-loop-tunneling"] = defaultdict(list)
        self.hamiltonian_elements["left-loop-env-tunneling"] = defaultdict(list)
        self.hamiltonian_elements["right-loop-env-tunneling"] = defaultdict(list)

        self.hamiltonian_skeleton = defaultdict(list)

        # Defining tools used here ############################################
        states_per_loop = self.cqps_twin_qubit_constant_manager.states_per_loop
        states_total_number = self.cqps_twin_qubit_constant_manager.states_total_number
        convert_state_to_index = (
            self.cqps_twin_qubit_constant_manager.convert_state_to_index
        )

        self.cqps_twin_qubit_constant_manager.print_constants()
        logging.debug(
            f"""🏗 Creating skeleton
{'States total number:':<30}{states_total_number}
{'States per island:':<30}{states_per_loop}"""
        )

        # Run through the matrix and generate the row, col, val elements
        for matrix_col in range(0, states_total_number):
            (
                state_numeric,
                state_flux,
            ) = self.cqps_twin_qubit_constant_manager.convert_index_to_state(matrix_col)

            # Evaluation of the energy will be done in stage 2 using the (nl, nr) values
            self.hamiltonian_elements["inductance"]["row"].append(matrix_col)
            self.hamiltonian_elements["inductance"]["col"].append(matrix_col)
            self.hamiltonian_elements["inductance"]["(nl,nr)"].append(
                (state_flux[LEFT], state_flux[RIGHT])
            )

            # Offdiagonal due to tunnelling
            # - check that index is within matrix boundaries
            # - create a row-col offset from diagonal
            # - fill out the symmetrical entries
            if state_numeric[LEFT] < (states_per_loop - 1):

                # Tunneling into/out of left loop
                matrix_row = convert_state_to_index(state_numeric + [1, 0])
                self.hamiltonian_elements["left-loop-env-tunneling"]["row"] += [
                    matrix_row,
                    matrix_col,
                ]
                self.hamiltonian_elements["left-loop-env-tunneling"]["col"] += [
                    matrix_col,
                    matrix_row,
                ]
                self.hamiltonian_elements["left-loop-env-tunneling"]["val"] += [
                    -self.cqps_twin_qubit_constant_manager.ES_left / 2,
                    -self.cqps_twin_qubit_constant_manager.ES_left / 2,
                ]

                # flux exchange |nl, nr> <-> |nl-1, nl+1> between left and RIGHT loops
                if state_numeric[RIGHT] > 0:
                    matrix_row = convert_state_to_index(state_numeric + [1, -1])
                    self.hamiltonian_elements["loop-loop-tunneling"]["row"] += [
                        matrix_col,
                        matrix_row,
                    ]
                    self.hamiltonian_elements["loop-loop-tunneling"]["col"] += [
                        matrix_row,
                        matrix_col,
                    ]
                    self.hamiltonian_elements["loop-loop-tunneling"]["val"] += [
                        -self.cqps_twin_qubit_constant_manager.ES_center / 2,
                        -self.cqps_twin_qubit_constant_manager.ES_center / 2,
                    ]

            if state_numeric[RIGHT] < (states_per_loop - 1):
                # Tunneling out from right tloop
                matrix_row = convert_state_to_index(state_numeric + [0, 1])
                self.hamiltonian_elements["right-loop-env-tunneling"]["row"] += [
                    matrix_row,
                    matrix_col,
                ]
                self.hamiltonian_elements["right-loop-env-tunneling"]["col"] += [
                    matrix_col,
                    matrix_row,
                ]
                self.hamiltonian_elements["right-loop-env-tunneling"]["val"] += [
                    -self.cqps_twin_qubit_constant_manager.ES_right / 2,
                    -self.cqps_twin_qubit_constant_manager.ES_right / 2,
                ]

        # Collect up all rows and columns and put them in the skeleton
        # Order matters - that is why an OrderedDict is used
        for component in self.hamiltonian_elements.keys():
            self.hamiltonian_skeleton["row"] += self.hamiltonian_elements[component][
                "row"
            ]
            self.hamiltonian_skeleton["col"] += self.hamiltonian_elements[component][
                "col"
            ]

        print_string = "🏗 Skeleton Hamiltonian Information\n"
        for component in self.hamiltonian_elements.keys():
            print_string += f"{component:<20}"
            for index in ["row", "col", "val"]:
                print_string += (
                    f"{index:<5}{len(self.hamiltonian_elements[component][index]):<7}"
                )
            print_string += "\n"
        logging.info(print_string)

    def stage_2_build_hamiltonian_for_simulation(
        self, phi_l: float, phi_r: float
    ) -> sp.coo_matrix:
        """
        phi_l and phi_r are scaled to be from 0-1
        """

        def inductance_energy_evaluator(nl_nr_tuple: Tuple[int, int]):
            return (
                nl_nr_tuple[LEFT] - phi_l
            ) ** 2 * self.cqps_twin_qubit_constant_manager.EL_left + (
                nl_nr_tuple[RIGHT] - phi_r
            ) ** 2 * self.cqps_twin_qubit_constant_manager.EL_right

        return sp.coo_matrix(
            (
                np.hstack(
                    (
                        [
                            inductance_energy_evaluator(i)
                            for i in self.hamiltonian_elements["inductance"]["(nl,nr)"]
                        ],
                        self.hamiltonian_elements["loop-loop-tunneling"]["val"],
                        self.hamiltonian_elements["left-loop-env-tunneling"]["val"],
                        self.hamiltonian_elements["right-loop-env-tunneling"]["val"],
                    )
                ),
                (
                    self.hamiltonian_skeleton["row"],
                    self.hamiltonian_skeleton["col"],
                ),
            )
        ).tocsr()
