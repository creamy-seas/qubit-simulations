"""
Class that stores all the constants and their conversions

Assume that each JJ is a 100x100nm² square
"""

import logging
import sys

import numpy as np

import pinject

pi = np.pi


class TransmonQubitConstantManager:
    @pinject.copy_args_to_public_fields
    def __init__(self, quantum_constants, generic_converter, param_dictionary):

        try:
            # State ###############################################################
            self.number_of_charge_states = param_dictionary["number_of_charge_states"]
            self.offset_to_apply = (self.number_of_charge_states - 1) // 2

            # Physical parameters that can be used to derive the energies #########
            self.jj_squares = param_dictionary["jj_squares"]

            self.C_jj = self.evaluate_C_jj(self.jj_squares)
            self.C_transmon = param_dictionary["C_transmon"]
            self.C_gate = param_dictionary["C_gate"]
            self.C_sigma = self.C_gate + 2 * self.C_jj + self.C_transmon

            self.EC = self.evaluate_EC(self.C_sigma)
            (self.EJ0, self.critical_current, self.resistance_jj) = self.evaluate_EJ0_and_critical_current(
                self.jj_squares
            )

            self.print_constants()

        except KeyError as err:
            logging.error(
                f"Need to pass in parameter ({err}) to the parameter dictionary"
            )
            sys.exit()

    def evaluate_C_jj(self, jj_squares: int) -> float:
        PERMITIVITY_ALOX = 10
        PERMITIVITY_VACUUM = 8.85 * 10 ** (-12)
        CAPACITANCE_AREA_FACTOR = 1.15
        THICKNESS_ALOX = 2 * 10 ** (-9)

        area_nm2 = jj_squares * 10 ** (-14)

        return (
            PERMITIVITY_ALOX
            * PERMITIVITY_VACUUM
            * CAPACITANCE_AREA_FACTOR
            * area_nm2
            / THICKNESS_ALOX
        )

    def evaluate_EC(self, system_capacitance: float) -> float:
        """EC = e^2 / (2 * C_sgima)"""

        return self.generic_converter.convert_energy_to_GHz(
            self.quantum_constants.eCharge ** 2 / (2 * system_capacitance)
        )

    def evaluate_EJ0_and_critical_current(self, jj_squares: int) -> (float, float):
        """
        I_critical = pi * Delta(0) / (2 * e * Rsquare * NSquares)

        EJ0 = h/(2e)^2 * Delta(0)/2 * 1 / (Rsquare * NSquares)
        """

        TRANSITION_ALUMINIUM = 1.3
        DELTA_ALUMINIUM = 1.73 * self.quantum_constants.kb * TRANSITION_ALUMINIUM
        ALUMINIUM_OXIDE_SHEET_RESISTANCE = 1.84 * 10 ** (3)

        resistance_jj = ALUMINIUM_OXIDE_SHEET_RESISTANCE / jj_squares

        I_critical = (
            pi * DELTA_ALUMINIUM / (2 * self.quantum_constants.eCharge * resistance_jj)
        )
        EJO = self.generic_converter.convert_energy_to_GHz(
            (
                self.quantum_constants.h
                * DELTA_ALUMINIUM
                / (8 * self.quantum_constants.eCharge ** 2 * resistance_jj)
            )
        )

        return (EJO, I_critical, resistance_jj)

    def evaluate_EJ_from_external_flux(self, flux_ext: float) -> float:
        """Flux in units of Φ0"""

        return self.EJ0 * 2 * np.cos(pi * flux_ext)

    def convert_index_to_state(self, index: int) -> int:
        if (index < 0) or (index >= self.number_of_charge_states):
            raise ValueError(
                f"Index {index} must be within 0 and {self.number_of_charge_states}"
            )

        return index - self.offset_to_apply

    def convert_external_flux_to_josephson_energy(self, ext_flux: float) -> float:
        """An external flux will result in a effective Josephson energy
        as a result of the superconducting loop in the device:

                EJ = EJ0 * 2 * |cos(pi * ext_flux / Phi0)|
        """
        pass

    def override_parameters(self, EC: float, EJ0: float):
        logging.info(
            f"""Overriding to the following values:
{'EC:':<10}{EC}
{'EJ0:':<10}{EJ0}
"""
        )

        self.EC = EC
        self.EJ0 = EJ0

    def print_constants(self):
        logging.info(
            f"""⚛ System setup with the following parameters (can choose to override them):

-----Energies-----
{'EC:':<50}{self.EC:.2f} (GHz)
{'EJ0:':<50}{self.EJ0:.2f} (GHz)

-----Raw Parameters-----
{'number_of_charge_states:':<50}{self.number_of_charge_states}
{'jj_squares:':<50}{self.jj_squares}
{'C_transmon:':<50}{self.C_transmon * 10**15:.2f} (fF)
{'C_gate:':<50}{self.C_gate * 10**15:.2f} (fF)

-----Derived Parameters-----
{'C_jj:':<50}{self.C_jj * 10**15:.2f} (fF)
{'C_𝛴 = C_transmon + C_gate + 2 * C_jj:':<50}{self.C_sigma * 10**15:.2f} (fF)
{'Critical Current:':<50}{self.critical_current * 10**6:.2f} (µA)
{'JJ Resistance:':<50}{self.resistance_jj / 10**3:.2f} (k𝛀)
"""
        )
