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
            "+phi1-phi2+phiExt",
            "-phi1+phi2-phiExt",
            "-phi2+phi3-nphiExt",
            "+phi2-phi3+nphiExt",
        ]
        self.hamiltonian_skeleton = defaultdict(lambda: defaultdict(list))
        self.hamiltonian_constant = {}
        self.hamiltonian_flux_dependent = {}

    def stage1_prepare_hamiltonian_skeleton(self):
        """
        __ Description __
        Generates building blocks of the Hamiltonian that will be scaled with energies in stage2 and stage3
        |--------------------+-------------------------------|
        |                    | Scaling                       |
        |--------------------+-------------------------------|
        | charge             | EC                            |
        | phi1               | -EJ/2                         |
        | phi2               | -alpha*EJ/2                   |
        | phi3               | -EJ/2                         |
        | +phi1-phi2+phiExt  | -EJ/2 * e^(+i * phi_ext)      |
        | -phi1+phi2-phiExt  | -EJ/2 * e^(-i * phi_ext)      |
        | -phi2+phi3-nphiExt | -EJ/2 * e^(-i * phi_ext * nu) |
        | +phi2-phi3+nphiExt | -EJ/2 * e^(+i * phi_ext * nu) |
        |--------------------+-------------------------------|
        """

        capacitance_mat_inv = self.twin_qubit_constant_manager.capacitance_mat_inv
        states_per_island = self.twin_qubit_state_manager.states_per_island
        states_total_number = self.twin_qubit_state_manager.states_total_number
        convert_index_to_island_state = (
            self.twin_qubit_state_manager.convert_index_to_island_state
        )
        convert_numeric_state_to_index = (
            self.twin_qubit_state_manager.convert_numeric_state_to_index
        )
        state1, state2, state3 = 0, 1, 2

        logging.debug(
            f"""🏗 Creating skeleton
{'States total number:':<30}{states_total_number}
{'States per island:':<30}{states_per_island}
{'Capacitance Matrix Inverse:'}
{capacitance_mat_inv}"""
        )

        # mcol for matrixColumn
        for mcol in range(0, states_total_number):

            (mcol_state_numeric, mcol_state_cp) = convert_index_to_island_state(mcol)

            # Evaluation of diagonal charging energy using nT * C-1 * n
            normalised_charging_energy = np.dot(
                mcol_state_cp, np.dot(capacitance_mat_inv, mcol_state_cp),
            )
            self.hamiltonian_skeleton["charge"]["row"].append(mcol)
            self.hamiltonian_skeleton["charge"]["col"].append(mcol)
            self.hamiltonian_skeleton["charge"]["elm"].append(
                normalised_charging_energy
            )

            # Offdiagonal JJ energy
            # - check that index is within matrix boundaries
            # - create a row-col offset from diagonal
            # - fill out the symmetrical entries
            if mcol_state_numeric[state1] < (states_per_island - 1):
                # cos(phi_1) => |n1><n1+1| and |n1+1><n1|
                mrow = convert_numeric_state_to_index(mcol_state_numeric + [1, 0, 0])
                self.hamiltonian_skeleton["phi1"]["row"] += [mrow, mcol]
                self.hamiltonian_skeleton["phi1"]["col"] += [mcol, mrow]

                # cos(phi_2 - phi_1 - phi_ext), with cp exchange between 2 <-> 1
                if mcol_state_numeric[state2] > 0:
                    mrow = convert_numeric_state_to_index(
                        mcol_state_numeric + [1, -1, 0]
                    )
                    self.hamiltonian_skeleton["-phi1+phi2-phiExt"]["row"] += [mcol]
                    self.hamiltonian_skeleton["-phi1+phi2-phiExt"]["col"] += [mrow]
                    self.hamiltonian_skeleton["+phi1-phi2+phiExt"]["row"] += [mrow]
                    self.hamiltonian_skeleton["+phi1-phi2+phiExt"]["col"] += [mcol]

            if mcol_state_numeric[1] < (states_per_island - 1):
                # alpha * cos(phi_2) => |n2><n2+1| and |n2+1><n2|
                mrow = convert_numeric_state_to_index(mcol_state_numeric + [0, 1, 0])
                self.hamiltonian_skeleton["phi2"]["row"] += [mrow, mcol]
                self.hamiltonian_skeleton["phi2"]["col"] += [mcol, mrow]

            if mcol_state_numeric[state3] < (states_per_island - 1):
                # cos(phi_3)
                mrow = convert_numeric_state_to_index(mcol_state_numeric + [0, 0, 1])

                self.hamiltonian_skeleton["phi3"]["row"] += [mrow, mcol]
                self.hamiltonian_skeleton["phi3"]["col"] += [mcol, mrow]

                # cos(phi_2 - phi_3 + nphi_ext), with cp exchange between 2 <-> 3
                if mcol_state_numeric[1] > 0:
                    mrow = convert_numeric_state_to_index(
                        mcol_state_numeric + [0, -1, 1]
                    )
                    self.hamiltonian_skeleton["+phi2-phi3+nphiExt"]["row"] += [mrow]
                    self.hamiltonian_skeleton["+phi2-phi3+nphiExt"]["col"] += [mcol]
                    self.hamiltonian_skeleton["-phi2+phi3-nphiExt"]["row"] += [mcol]
                    self.hamiltonian_skeleton["-phi2+phi3-nphiExt"]["col"] += [mrow]

        self.convert_skeleton_charge_elements_to_array()
        self.set_skeleton_jj_elements_as_ones()

        self.print_skeleton_information()

    def convert_skeleton_charge_elements_to_array(self):
        self.hamiltonian_skeleton["charge"]["elm"] = np.array(
            self.hamiltonian_skeleton["charge"]["elm"]
        )

    def set_skeleton_jj_elements_as_ones(self):
        for component in self.hamiltonian_components[1:]:
            self.hamiltonian_skeleton[component]["elm"] = np.ones(
                len(self.hamiltonian_skeleton[component]["row"])
            )

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
        | ✘          | +phi1-phi2+phiExt-row | +phi1-phi2+phiExt-col |   |            ✘		   |
        | ✘          | -phi1+phi2-phiExt-row | -phi1+phi2-phiExt-col |   |            ✘		   |
        | ✘          | -phi2+phi3-nphiExt-row | -phi2+phi3-nphiExt-col |   |            ✘		   |
        | ✘          | +phi2-phi3+nphiExt-row | +phi2-phi3+nphiExt-col |   |            ✘		   |
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

        # Rows and Columns are stacked - order is critical
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
        | +phi1-phi2+phiExt |            |            |   | -EJ/2 * e^(+i * phi_ext)      |
        | -phi1+phi2-phiExt |            |            |   | -EJ/2 * e^(-i * phi_ext)      |
        | -phi2+phi3-nphiExt |            |            |   | -EJ/2 * e^(-i * phi_ext * nu) |
        | +phi2-phi3+nphiExt |            |            |   | -EJ/2 * e^(+i * phi_ext * nu) |
        |--------+------------+------------+---+-------------------------------|
        """

        EJ = self.twin_qubit_constant_manager.EJ

        self.hamiltonian_flux_dependent["elm"] = np.hstack(
            (
                -EJ
                / 2
                * np.exp(-1j * phi_external)
                * self.hamiltonian_skeleton["+phi1-phi2+phiExt"]["elm"],
                -EJ
                / 2
                * np.exp(+1j * phi_external)
                * self.hamiltonian_skeleton["-phi1+phi2-phiExt"]["elm"],
                -EJ
                / 2
                * np.exp(-1j * phi_external_assymetric)
                * self.hamiltonian_skeleton["-phi2+phi3-nphiExt"]["elm"],
                -EJ
                / 2
                * np.exp(+1j * phi_external_assymetric)
                * self.hamiltonian_skeleton["+phi2-phi3+nphiExt"]["elm"],
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
