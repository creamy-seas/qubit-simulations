from functools import partial
from functools import reduce
import operator
from collections import defaultdict
import logging

import pinject
import numpy as np
import scipy.sparse as sp


class TwinQubitHamiltonianManager:
    """Builds and stores a custom Hamiltonian for the twin qubit"""

    @pinject.copy_args_to_public_fields
    def __init__(
        self, twin_qubit_state_manager, twin_qubit_constant_manager,
    ):
        self.hamiltonian_components = [
            "charge",
            "phi1",
            "phi2",
            "phi3",
            "+phi21",
            "-phi21",
            "+phi32",
            "-phi32",
        ]
        self.hamiltonian_skeleton = defaultdict(lambda: defaultdict(list))
        self.hamiltonian_constant = {}
        self.hamiltonian_flux_dependent = {}
        self.hamiltonian_simulation = None

    def stage1_prepare_hamiltonian_skeleton(self):
        """
        __ Description __
        Generates building blocks of the Hamiltonian that will be scaled with energies in stage2 and stage3
        |------------+------------+------------+---+-------------------------------|
        |            |            |            |   | Scaling of elm                |
        |------------+------------+------------+---+-------------------------------|
        | charge-elm | charge-row | charge-col |   | EC                            |
        | phi1-elm   | phi1-row   | phi1-col   |   | -EJ/2                         |
        | phi2-elm   | phi2-row   | phi2-col   |   | -alpha*EJ/2                   |
        | phi3-elm   | phi3-row   | phi3-col   |   | -EJ/2                         |
        | +phi21-elm | +phi21-row | +phi21-col |   | -EJ/2 * e^(+i * phi_ext)      |
        | -phi21-elm | -phi21-row | -phi21-col |   | -EJ/2 * e^(-i * phi_ext)      |
        | +phi32-elm | +phi32-row | +phi32-col |   | -EJ/2 * e^(-i * phi_ext * nu) |
        | -phi32-elm | -phi32-row | -phi32-col |   | -EJ/2 * e^(+i * phi_ext * nu) |
        |------------+------------+------------+---+-------------------------------|
        """

        capacitance_mat_inv = self.twin_qubit_constant_manager.capacitance_mat_inv
        states_per_island = self.twin_qubit_state_manager.states_per_island
        states_total_number = self.twin_qubit_state_manager.states_total_number
        convert_index_to_island_state = (
            self.twin_qubit_state_manager.convert_index_to_island_state
        )
        convert_island_state_to_index = (
            self.twin_qubit_state_manager.convert_island_state_to_index
        )

        logging.debug(
            f"""🏗 Creating skeleton
{'States total number:':<30}{states_total_number}
{'States per island:':<30}{states_per_island}
{'Capacitance Matrix Inverse:'}
{capacitance_mat_inv}"""
        )

        for x in range(0, states_total_number):

            (state_indx, state_cp) = convert_index_to_island_state(x)

            # Evaluation using nT * C-1 * n
            normalised_charging_energy = np.dot(
                state_cp, np.dot(capacitance_mat_inv, state_cp),
            )
            self.hamiltonian_skeleton["charge"]["row"].append(x)
            self.hamiltonian_skeleton["charge"]["col"].append(x)
            self.hamiltonian_skeleton["charge"]["elm"].append(
                normalised_charging_energy
            )

            # Offdiagonal JJ energy
            # - check that index is within matrix boundaries
            # - offset y coordinate from diagonal
            # - fill out the symmetrical entries
            if state_indx[0] < (states_per_island - 1):
                # cos(phi_1) => |n1><n1+1| and |n1+1><n1|
                y = convert_island_state_to_index(state_indx + [1, 0, 0])
                self.hamiltonian_skeleton["phi1"]["row"] += [x, y]
                self.hamiltonian_skeleton["phi1"]["col"] += [y, x]

                # cos(phi_2 - phi_1 - phi_ext), with cp exchange between 2 <-> 1
                if state_indx[1] > 0:
                    y = convert_island_state_to_index(state_indx + [1, -1, 0])
                    self.hamiltonian_skeleton["-phi21"]["row"] += [y]
                    self.hamiltonian_skeleton["-phi21"]["col"] += [x]
                    self.hamiltonian_skeleton["+phi21"]["row"] += [x]
                    self.hamiltonian_skeleton["+phi21"]["col"] += [y]

            if state_indx[1] < (states_per_island - 1):
                # alpha * cos(phi_2) => |n2><n2+1| and |n2+1><n2|
                y = convert_island_state_to_index(state_indx + [0, 1, 0])
                self.hamiltonian_skeleton["phi2"]["row"] += [x, y]
                self.hamiltonian_skeleton["phi2"]["col"] += [y, x]

            if state_indx[2] < (states_per_island - 1):
                # cos(phi_3) => |n3><n3+1| and |n3+1><n3|
                y = convert_island_state_to_index(state_indx + [0, 0, 1])
                self.hamiltonian_skeleton["phi3"]["row"] += [x, y]
                self.hamiltonian_skeleton["phi3"]["col"] += [y, x]

                # cos(phi_2 - phi_3 + nphi_ext), with cp exchange between 2 <-> 3
                if state_indx[1] > 0:
                    y = convert_island_state_to_index(state_indx + [0, -1, 1])
                    self.hamiltonian_skeleton["-phi32"]["row"] += [x]
                    self.hamiltonian_skeleton["-phi32"]["col"] += [y]
                    self.hamiltonian_skeleton["+phi32"]["row"] += [y]
                    self.hamiltonian_skeleton["+phi32"]["col"] += [x]

        # For kinetic charge energy, the elements have been evaluated
        self.hamiltonian_skeleton["charge"]["elm"] = np.array(
            self.hamiltonian_skeleton["charge"]["elm"]
        )
        # For potential JJ energy, we set the elements to 1
        for component in self.hamiltonian_components[1:]:
            self.hamiltonian_skeleton[component]["elm"] = np.ones(
                len(self.hamiltonian_skeleton[component]["row"])
            )

        self.print_skeleton_information()

    def print_skeleton_information(self):
        warn_string = "🏗 Skeleton Hamiltonian Information\n"
        for component in self.hamiltonian_components:
            warn_string += f"{component:<20}"
            for index in ["row", "col", "elm"]:
                warn_string += (
                    f"{index:<5}{len(self.hamiltonian_skeleton[component][index]):<7}"
                )
            warn_string += "\n"

        logging.info(warn_string)

    def stage2_prepare_constant_hamiltonian(self):
        """
        __ Description __
        Scale the skeleton by EC, EJ, alpha and stack the row and column indicies
        |------------+------------+------------+---+-------------------------------|
        |            |            |            |   | Scaling of elm                |
        |------------+------------+------------+---+-------------------------------|
        | charge-elm | charge-row | charge-col |   | EC                            |
        | phi1-elm   | phi1-row   | phi1-col   |   | -EJ/2                         |
        | phi2-elm   | phi2-row   | phi2-col   |   | -alpha * EJ/2                 |
        | phi3-elm   | phi3-row   | phi3-col   |   | -EJ/2                         |
        | ✘          | +phi21-row | +phi21-col |   |            ✘		   |
        | ✘          | -phi21-row | -phi21-col |   |            ✘		   |
        | ✘          | +phi32-row | +phi32-col |   |            ✘		   |
        | ✘          | -phi32-row | -phi32-col |   |            ✘		   |
        |------------+------------+------------+---+-------------------------------|
        | ['elm']    | ['row']    | ['col']    |   |                               |
        |------------+------------+------------+---+-------------------------------|
        """

        EC = self.twin_qubit_constant_manager.EC
        EJ = self.twin_qubit_constant_manager.EJ
        alpha = self.twin_qubit_constant_manager.alpha
        logging.info(
            f"""🏗 Scaling constant parts of Hamiltonian with
{'EC:':<10}{EC}
{'EJ:':<10}{EJ}
{'alpha:':<10}{alpha}"""
        )

        # Rows and Columns are stacked
        self.hamiltonian_constant["row"] = []
        self.hamiltonian_constant["col"] = []
        for component in self.hamiltonian_components:
            self.hamiltonian_constant["row"] += self.hamiltonian_skeleton[component][
                "row"
            ]
            self.hamiltonian_constant["col"] += self.hamiltonian_skeleton[component][
                "col"
            ]

        # Constant elements are scaled by EC, -EJ/2, -alpha * EJ/2
        self.hamiltonian_constant["elm"] = np.hstack(
            (
                EC * self.hamiltonian_skeleton["charge"]["elm"],
                -EJ / 2 * self.hamiltonian_skeleton["phi1"]["elm"],
                -alpha * EJ / 2 * self.hamiltonian_skeleton["phi2"]["elm"],
                -EJ / 2 * self.hamiltonian_skeleton["phi3"]["elm"],
            )
        )

        # Check number of elements produced (kind of redundant)
        check = reduce(
            operator.sub,
            [
                len(self.hamiltonian_constant["elm"]),
                len(self.hamiltonian_skeleton["charge"]["elm"]),
                len(self.hamiltonian_skeleton["phi1"]["elm"]),
                len(self.hamiltonian_skeleton["phi2"]["elm"]),
                len(self.hamiltonian_skeleton["phi3"]["elm"]),
            ],
        )
        if check:
            raise ValueError(
                f"{check} excess elements being produced for {len(self.hamiltonian_constant['elm'])}"
            )

        logging.info(
            f"""🏗 Constructed matrix 'hamiltonian_constant' with
{len(self.hamiltonian_constant['elm']):<10} elements (rest will be evaluated during stage3)
{len(self.hamiltonian_constant['row']):<10} rows
{len(self.hamiltonian_constant['col']):<10} cols
"""
        )

    def stage3_build_hamiltonian_for_simulation(
        self, phi_external, phi_external_assymetric
    ):
        """
        hamiltonian_constant['elm'] + hamiltonian_flux_dependent['elm']
        hamiltonian_constant['row']
        hamiltonian_constant['col']

        and construct a sparse matrix

        __ Parameters __
        [float] phi_external, phi_external_assymetric:  fluxes to the two loops

        __ Description __
        Scale the flux-dependent part of the Hamiltonian changes as the field is being swept.
        |--------+------------+------------+---+-------------------------------|
        | +phi21 |            |            |   | -EJ/2 * e^(+i * phi_ext)      |
        | -phi21 |            |            |   | -EJ/2 * e^(-i * phi_ext)      |
        | +phi32 |            |            |   | -EJ/2 * e^(-i * phi_ext * nu) |
        | -phi32 |            |            |   | -EJ/2 * e^(+i * phi_ext * nu) |
        |--------+------------+------------+---+-------------------------------|
        """

        EJ = self.twin_qubit_constant_manager.EJ

        self.hamiltonian_flux_dependent["elm"] = np.hstack(
            (
                -EJ
                / 2
                * np.exp(1j * phi_external)
                * self.hamiltonian_skeleton["+phi21"]["elm"],
                -EJ
                / 2
                * np.exp(-1j * phi_external)
                * self.hamiltonian_skeleton["+phi21"]["elm"],
                -EJ
                / 2
                * np.exp(-1j * phi_external_assymetric)
                * self.hamiltonian_skeleton["+phi32"]["elm"],
                -EJ
                / 2
                * np.exp(1j * phi_external_assymetric)
                * self.hamiltonian_skeleton["-phi32"]["elm"],
            )
        )

        logging.debug(
            f"Making {len(self.hamiltonian_flux_dependent['elm'])} flux-dependent elements"
        )

        self.hamiltonian_simulation = sp.coo_matrix(
            (
                np.hstack(
                    (
                        self.hamiltonian_constant["elm"],
                        self.hamiltonian_flux_dependent["elm"],
                    )
                ),
                (self.hamiltonian_constant["row"], self.hamiltonian_constant["col"]),
            )
        ).tocsr()
