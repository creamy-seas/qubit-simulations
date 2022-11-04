from typing import Dict, List, Tuple
import logging

import pinject
import numpy as np


class TwinQubitConstantManager:

    @pinject.copy_args_to_public_fields
    def __init__(self, quantum_constants, generic_converter, param_dictionary):

        # 1 - Set the constants
        self.alpha = param_dictionary["alpha"]
        self.assymetry = param_dictionary["assymetry"]
        self.jj_squares = param_dictionary["jj_squares"]
        self.jj_overlap_area = 200 * 200 * self.jj_squares
        self.delta = 0

        # 2 - set derived constatns
        (self.EC, self.jj_capacitance) = self.evaluate_ec_and_capacitance(
            self.jj_overlap_area
        )

        (
            self.EJ,
            self.jj_critical_current,
            self.jj_resistance,
        ) = self.evaluate_ej_criticalcurrent_and_resistance(self.jj_squares)

        (
            self.capacitance_mat,
            self.capacitance_mat_inv,
        ) = self.evaluate_capacitance_matrix(self.alpha, self.delta)

    def evaluate_capacitance_matrix(
        self, alpha: float, delta: float
    ) -> Tuple[np.array, np.array]:

        capacitance_matrix = np.matrix(
            [[2, -1, 0], [-1, 2 + alpha + delta, -1], [0, -1, 2]]
        )

        logging.debug(
            f"""⚛ Created capacitance matrix and it's inverse
C={np.array(capacitance_matrix)}
--------------------
C_inv={np.array(capacitance_matrix.I)}
--------------------
C*C_inv = {np.array(capacitance_matrix).dot(np.array(capacitance_matrix.I))}
"""
        )
        return (np.array(capacitance_matrix), np.array(capacitance_matrix.I))

    def evaluate_ec_and_capacitance(self, area_nm2: float) -> Tuple[float, float]:
        """
        __ Parameters __
        area: JJ area in nm^2

        __ Description __
        Evaluating charging energy and capacitance
        EC:             (2e)^2/2C
        capacitance:    C = e*e_0*A/d   for normalisation

        __ Returns __
        (EC, capacitance)
        """

        THICKNESS_ALOX = 2 * 10 ** (-9)
        PERMITIVITY_ALOX = 10
        PERMITIVITY_VACUUM = 8.85 * 10 ** (-12)
        CAPACITANCE_AREA_OFFSET = 1.15 * 10 ** (-18)  # convert nm2 to m2

        jj_capacitance = (
            PERMITIVITY_ALOX
            * PERMITIVITY_VACUUM
            * CAPACITANCE_AREA_OFFSET
            * area_nm2
            / THICKNESS_ALOX
        )
        EC = (2 * self.quantum_constants.eCharge) ** 2 / (2 * jj_capacitance)
        EC = self.generic_converter.convert_energy_to_GHz(EC)

        logging.debug(f"⚛ Evaluated EC: {EC}, jj_capacitance: {jj_capacitance}")
        return (EC, jj_capacitance)

    def evaluate_ej_criticalcurrent_and_resistance(self, squares_JJ):
        """
        __ Parameters __
        squares:        width of the JJ in square units

        __ Description __
        The more squares of JJ, the bigger the resistance, and the larger it's energy
        We assume that it is a standard JJ: 20nm - AlOx - 30nm
        Count in 100x100nm^2 squares

        __ Returns __
        (EJ, jj_critical_current, jj_resistance)
        """
        TRANSITION_AL = 1.2
        BOLTZMAN = 1.38 * 10 ** (-23)
        DELTA_AL = 1.764 * TRANSITION_AL * BOLTZMAN
        RESISTANCE_CONSTANT = self.quantum_constants.h / (
            4 * self.quantum_constants.eCharge ** 2
        )
        RESISTANCE_OF_100x100_AL = 18.4 * 10 ** 3

        jj_resistance = RESISTANCE_OF_100x100_AL / squares_JJ
        jj_critical_current = (
            np.pi * DELTA_AL / (2 * self.quantum_constants.eCharge * jj_resistance)
        )
        EJ = RESISTANCE_CONSTANT * DELTA_AL / (2 * jj_resistance)
        EJ = self.generic_converter.convert_energy_to_GHz(EJ)

        logging.debug(
            f"⚛ Evaluated EJ: {EJ}, jj_critical_current: {jj_critical_current}, jj_resistance: {jj_resistance}"
        )
        return (EJ, jj_critical_current, jj_resistance)

    def override_parameters(self, EC=None, EJ=None, alpha=None, assymetry=None):
        """Override parameters thare are not None"""

        if EC:
            logging.debug(f"⚛ Overriding {'EC:':<10}{self.EC:.4f} -> {EC:.4f}")
            self.EC = EC
        if EJ:
            logging.debug(f"⚛ Overriding {'EJ:':<10}{self.EJ:.4f} -> {EJ:.4f}")
            self.EJ = EJ
        if alpha:
            logging.debug(f"⚛ Overriding {'alpha:':<10}{self.alpha:.4f} -> {alpha:.4f}")
            self.alpha = alpha
        if assymetry:
            logging.debug(
                f"⚛ Overriding {'assymetry:':<10}{self.assymetry:.4f} -> {assymetry:.4f}"
            )
            self.assymetry = assymetry

    def print_constants(self):
        logging.info(
            f"""⚛ Constant Manager using parameters:
{'EC:':<30}{self.EC}
{'EJ:':<30}{self.EJ}
{'alpha:':<30}{self.alpha}
{'assymetry:':<30}{self.assymetry}
{'jj_critical_current:':<30}{self.jj_critical_current}
{'jj_resistance:':<30}{self.jj_resistance}
{'jj_capacitance:':<30}{self.jj_capacitance}
{'capacitance_matrix:':<30}
{self.capacitance_mat}"""
        )
