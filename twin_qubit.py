from collections import defaultdict
import logging
from typing import List, Dict

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt


from .master import MasterQubit
import utils


class TwinQubit(MasterQubit):
    """
    Quantum class to work with the twin, 5JJ qubit.
    """

    def __init__(
        self,
        param_dictionary: Dict,
        flux_list: List,
        plot_or_not=False,
        logging_level=False,
    ):
        """
        __ Parameters __
        [dict] param_dictionary:        {[float] 'alpha':               middle junction ration to outside ones
                                         [float] 'assymetry':           of the flux loops
                                         [int] 'states_total_number':     for the simulation
                                         [1D-int] flux_list             flux list to evalute for in units of Phi0}
        __ Description __
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        MasterQubit.__init__(self, plot_or_not, logging_level)

        # 1 - store user supplied parameters
        self.EC = param_dictionary["EC"]
        self.EJ = param_dictionary["EJ"]
        self.alpha = param_dictionary["alpha"]
        self.assymetry = param_dictionary["assymetry"]
        self.states_per_island = int(param_dictionary["states_per_island"])
        self.states_total_number = self.states_per_island ** 3
        self.flux_list = flux_list
        # --------------------
        self.delta = 0

        # 2 - very common parameters, which should only be set once
        self.prepare_correction()
        self.cp_offset_to_apply = (self.states_per_island - 1) // 2

        # 3 - simulation preparation
        self.prepare_structure()
        self.prepare_hamiltonian_stage1()

        # 4 - plot preparation
        self.fig, self.ax = self.prepare_plot(1, 1)

    def override_parameters(
        self, EC=None, EJ=None, alpha=None, assymetry=None, flux_list=None
    ):
        """
        __ Parameters __
        EC, EJ, alpha, assymetry: parameters to set

        __ Description __
        For new simulations, set the new system parameters
        """

        logging.info("⚙ Overriding parameters:")
        if EC:
            self.EC = EC
            logging.info("> {'EC':<10} = {EC:.4f}")
        if EJ:
            self.EJ = EJ
            logging.info("> {'EJ':<10} = {EJ:.4f}")
        if alpha:
            self.alpha = alpha
            logging.info("> {'alpha':<10} = {alpha:.4f}")
        if assymetry:
            self.assymetry = assymetry
            logging.info("> {'assymetry':<10} = {assymtery:.4f}")
        if flux_list:
            self.flux_list = flux_list
            logging.info("> Updating flux list")

    def prepare_correction(self):
        """
        Correct common errors before program begins
        """
        if (self.states_total_number % 2) == 0:
            raise ValueError(
                f"Total number of states is {self.states_total_number} -> it needs to be odd"
            )

    def prepare_structure(self):
        """
        Prepares the structure parameters
        """
        logging.info("'prepare_structure' creating energies and capacitances")

        # 1 - set jj dimensions
        self.param_jj_squares = 2
        self.param_jj_overlap_area = 200 * 200 * self.param_jj_squares

        # 2 - set the energies EC and EJ
        self.set_charging_energy(self.param_jj_overlap_area)
        self.set_josephson_energy(self.param_jj_squares)

        # 3 - set capacitance matrix
        self.capacitance_normalised = np.matrix(
            [[2, -1, 0], [-1, 2 + self.alpha + self.delta, -1], [0, -1, 2]]
        )

    def prepare_operators(self):
        """
        __ Description __
        Prepare operators for evaluating expectation values once the eigenstates are found

        __ Return __
        self.op_V
        self.op_Phi
        """

        op_V_dict = defaultdict(list)
        op_Phi_dict = defaultdict(list)
        op_V_constant = (
            self.const_h
            * 10 ** 9
            * self.EC
            / (2 * self.const_eCharge * (1 + self.alpha))
        )

        # additional operator components
        op_V_
        self.op_V_row = []
        self.op_V_col = []
        self.op_V_elm = []
        self.op_Phi_row = []
        self.op_Phi_col = []
        self.op_Phi_elements = []

        for x in range(0, self.states_total_number):

            # 2 - convert the index to island occupation
            self.convert_index_to_island_state(x)

            # 3 - voltage operator (convert EC to Joules)
            voltage_element = (
                self.const_h
                * 10 ** 9
                * self.EC
                / (2 * self.const_eCharge * (1 + self.alpha))
                * (
                    self.state_cp_distribution[0, 0]
                    + 2 * self.state_cp_distribution[0, 1]
                    + self.state_cp_distribution[0, 2]
                )
            )
            self.op_V_row.append(x)
            self.op_V_col.append(x)
            self.op_V_elm.append(voltage_element)

            # 4 - phase operator phi_20
            if self.state_numerical_distribution[1] < (self.states_total_number - 1):
                # island 2 (element 1)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [0, 1, 0]
                )

                self.op_Phi_row.append(x)
                self.op_Phi_col.append(y)
                self.op_Phi_elements.append(1)

        # 5 - finalise operators
        self.op_V = sp.coo_matrix(
            (self.op_V_elm, (self.op_V_row, self.op_V_col))
        ).tocsr()
        self.op_Phi = sp.coo_matrix(
            (self.op_Phi_elements, (self.op_Phi_row, self.op_Phi_col))
        ).tocsr()

    def prepare_hamiltonian_stage1(self):
        """
        Generate elements for a normalised Hamiltonian.

        During the simulation 'prepare_hamiltonian()' must be run, to change energies
        """
        logging.info("'prepare_normalised_hamiltonian' creating Hamiltonian entries")

        # 0 - Individual Hamiltonian components.
        # constant
        self.op_H_charging_elm = []
        self.op_H_charging_row = []
        self.op_H_charging_col = []
        self.op_H_diag_row = []
        self.op_H_diag_col = []
        self.op_H_diagA_row = []
        self.op_H_diagA_col = []
        # changing
        self.op_H_phi_row = []
        self.op_H_phi_col = []
        self.op_H_phiAss_row = []
        self.op_H_phiAss_col = []

        # 1 - generate matrix elements, by going across all the states
        for x in range(0, self.states_total_number):

            # 2 - convert the index to island occupation
            self.convert_index_to_island_state(x)

            # 3 - diagonal charging energy
            charging_energy = (
                self.state_cp_distribution
                * (self.capacitance_normalised.I)
                * (self.state_cp_distribution.transpose())
            )

            self.op_H_charging_row.append(x)
            self.op_H_charging_col.append(x)
            self.op_H_charging_elm.append(charging_energy[0, 0])

            # 4 - offdiagonal JJ energy - check that within matrix boundaries
            # offset y coordinate from diagonal, and fill out the symmetrical entries
            # ! 'diag_elm' array is not created, since all elements will be the same
            if self.state_numerical_distribution[0] < (self.states_total_number - 1):
                # cos(phi_1)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [1, 0, 0]
                )
                self.op_H_diag_row.extend([x, y])
                self.op_H_diag_col.extend([y, x])

                # cos(phi_2 - phi_1 - phi_ext), with cp exchange between 2 <-> 1
                if self.state_numerical_distribution[1] > 0:
                    y = self.convert_numerical_state_to_index(
                        self.state_numerical_distribution + [1, -1, 0]
                    )
                    self.op_H_phi_row.extend([x, y])
                    self.op_H_phi_col.extend([y, x])

            if self.state_numerical_distribution[1] < (self.states_total_number - 1):
                # alpha * cos(phi_2)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [0, 1, 0]
                )
                self.op_H_diagA_row.extend([x, y])
                self.op_H_diagA_col.extend([y, x])

            if self.state_numerical_distribution[2] < (self.states_total_number - 1):
                # cos(phi_3)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [0, 0, 1]
                )
                self.op_H_diag_row.extend([x, y])
                self.op_H_diag_col.extend([y, x])

                # cos(phi_2 - phi_3 + nphi_ext), with cp exchange between 2 <-> 3
                if self.state_numerical_distribution[1] > 0:
                    y = self.convert_numerical_state_to_index(
                        self.state_numerical_distribution + [0, -1, 1]
                    )
                    self.op_H_phiAss_row.extend([x, y])
                    self.op_H_phiAss_col.extend([y, x])

        logging.warning(
            "Unchaning part of Hamiltonian has %i entries",
            len(self.op_H_charging_row + self.op_H_diag_row + self.op_H_diagA_row),
        )
        logging.warning(
            "Flux-dependent part of Hamiltonian has %i entries",
            len(self.op_H_phiAss_row + self.op_H_phiAss_row),
        )
        logging.info("'prepare_normalised_hamiltonian' finished")

    def prepare_hamiltonian_stage2(self):
        """
        __ Description __
        Using supplied EC, EJ, alpha, scale the normalised list created in 'prepare_normalised_hamiltonian' and
        combine them into 1.

        The lists are used in 'simulate' functions
        """

        logging.info(
            "==> 'prepare_hamiltonian' with EC=%.4f\tEJ=%.4f\talpha=%.4f\tass=%.4f",
            self.EC,
            self.EJ,
            self.alpha,
            self.assymetry,
        )

        # 1 - main part of the Hamiltonian, which remains unchanged during simulation
        temp_charging_elm = self.EC * np.array(self.op_H_charging_elm)
        temp_diag_elm = -self.EJ / 2 * np.ones(len(self.op_H_diag_row))
        temp_diagA_elm = -self.alpha * self.EJ / 2 * np.ones(len(self.op_H_diagA_row))

        self.op_H_elm = (
            list(temp_charging_elm) + list(temp_diag_elm) + list(temp_diagA_elm)
        )
        self.op_H_row = (
            self.op_H_charging_row + self.op_H_diag_row + self.op_H_diagA_row
        )
        self.op_H_col = (
            self.op_H_charging_col + self.op_H_diag_col + self.op_H_diagA_col
        )

        if (len(self.op_H_row) != len(self.op_H_elm)) or (
            len(self.op_H_col) != len(self.op_H_elm)
        ):
            self.raise_error(
                "Hamiltonin has %i rows, %i columns non zero coordinatres and only %i elements"
                % (len(self.op_H_row), len(self.op_H_col), len(self.op_H_elm))
            )

        # 2 - exchange elements - used by 'prepare_hamiltonian_stage3' at different fluxes
        no_exchange = int(len(self.op_H_phi_row) / 2)
        no_exchangeAss = int(len(self.op_H_phiAss_row) / 2)
        if no_exchange != no_exchangeAss:
            self.raise_error(
                "Number of 1 <-> 2 exchanges (%i) different from 2 <-> 3 exchanges (%i)"
                % (no_exchange, no_exchangeAss)
            )
        self.op_H_SUPPORT_exchange_elm = np.ones(no_exchange)

        logging.info("'prepare_hamiltonian' finished")

    def prepare_hamiltonian_stage3(self, phi_external, phi_externalAss):
        """
        __ Parameters __
        phi_external, phi_externalAss: fluxes to the two loops

        __ Description __
        The exchange parts of the Hamiltonian depend on the external fluxes, which changes for
        each simulation as the field is being swept. This method fills out a list
        of exchange values that needs to be added to the Hamiltonian list:

                    op_H_elm.extend(elm_list)
                    op_H_row.extend(self.op_H_phi_row)
                    op_H_col.extend(self.op_H_phi_col)

        __ Return __
        List that should extend the Hamiltonian in each simulation run
        """
        temp_phi_p = (
            -self.EJ / 2 * np.exp(1j * phi_external)
        ) * self.op_H_SUPPORT_exchange_elm
        temp_phi_n = (
            -self.EJ / 2 * np.exp(-1j * phi_external)
        ) * self.op_H_SUPPORT_exchange_elm
        temp_phiAss_n = (
            -self.EJ / 2 * np.exp(-1j * phi_externalAss)
        ) * self.op_H_SUPPORT_exchange_elm
        temp_phiAss_p = (
            -self.EJ / 2 * np.exp(1j * phi_externalAss)
        ) * self.op_H_SUPPORT_exchange_elm

        # 2 - interleave the lists, alternating between +ve and -ve fluxes, for the (x,y) and (y,x)
        # pairings in the Hamiltonian
        elm_list = list(np.ravel(np.column_stack((temp_phi_p, temp_phi_n))))
        elm_list = elm_list + list(
            np.ravel(np.column_stack((temp_phiAss_n, temp_phiAss_p)))
        )

        return elm_list

    def simulate(self, simulate_expectations):
        """
        __ Parameters __
        [array] simulate_expectations of the form [voltage=False, phi=False]

        __ Description __
        Method performs the eigenvalue simulations:
        1 - finish building the Hamiltonian, with the flux-dependent elements
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        3 - plot out the spectrum

        """
        simulate_voltage = simulate_expectations[0]
        simulate_phi = simulate_expectations[1]

        # 0 - prepare hamiltonian for this simulation
        self.prepare_hamiltonian_stage2()

        logging.info("'simulate' running")

        self.spectrum_eigvals = []
        self.spectrum_eigvecs = []
        self.spectrum_simulation_12 = []
        self.spectrum_simulation_23 = []
        self.dipole_moment_voltage = defaultdict(list)
        self.dipole_moment_voltage_beta = []
        self.dipole_moment_phi = []

        # copy main parts of Hamiltonian
        op_H_row = self.op_H_row.copy()
        op_H_col = self.op_H_col.copy()
        op_H_row.extend(self.op_H_phi_row)
        op_H_row.extend(self.op_H_phiAss_row)
        op_H_col.extend(self.op_H_phi_col)
        op_H_col.extend(self.op_H_phiAss_col)

        # 1 - simulation performed for each value of the external flux
        for ext_flux_number in range(0, len(self.flux_list)):

            # 2 - extraction of the onset phase, due to external flux
            phi_external = (self.flux_list[ext_flux_number]) * 2 * np.pi
            phi_externalAss = phi_external * self.assymetry

            ####################
            # 3 - FINISH MAIN HAMILTONIAN for this particular bias
            ####################
            # a - copy constant terms
            op_H_elm = self.op_H_elm.copy()

            # b - add on the phase dependent elements
            op_H_elm.extend(
                self.prepare_hamiltonian_stage3(phi_external, phi_externalAss)
            )

            if (len(op_H_row) != len(op_H_elm)) or (len(op_H_col) != len(op_H_elm)):
                self.raise_error(
                    "Hamiltonin lists have %i rows, %i columns  and only %i elements"
                    % (len(op_H_row), len(op_H_col), len(op_H_elm))
                )

            # 2 - construct sparse matrix
            self.op_H = sp.coo_matrix((op_H_elm, (op_H_row, op_H_col))).tocsr()

            # 3 - evaluate 3 lowest eigenenergies and eigenvectors, |1> |2> |3>
            eigvals, eigvecs = eigsh(self.op_H, 3, which="SA", tol=0)
            # sort in ascending order
            sort_index = np.argsort(eigvals)
            eigvals = np.array(eigvals)[sort_index]
            eigvecs = np.transpose(np.transpose(eigvecs)[sort_index])

            self.spectrum_eigvals.append(eigvals)
            self.spectrum_simulation_12.append([eigvals[1] - eigvals[0]])
            self.spectrum_simulation_23.append([eigvals[2] - eigvals[1]])

            # 4 - dipole transition
            if simulate_voltage:
                # evalute voltage in a sandwhich between the ground and excited states
                state_0 = eigvecs[:, 0]
                state_1 = eigvecs[:, 1]
                state_2 = eigvecs[:, 2]

                matrix01 = state_0.dot(self.op_V.dot(state_1))
                matrix02 = state_0.dot(self.op_V.dot(state_2))
                matrix12 = state_1.dot(self.op_V.dot(state_2))

                temp_dipole_moment_voltage_beta = (
                    matrix01
                    / self.const_Phi0
                    / self.spectrum_simulation_12[ext_flux_number]
                )

                self.dipole_moment_voltage["01"].append([matrix01.real, matrix01.imag])
                self.dipole_moment_voltage["02"].append([matrix02.real, matrix02.imag])
                self.dipole_moment_voltage["12"].append([matrix12.real, matrix12.imag])
                self.dipole_moment_voltage_beta.append(
                    [
                        temp_dipole_moment_voltage_beta.real,
                        temp_dipole_moment_voltage_beta.imag,
                    ]
                )

            self.track_progress(ext_flux_number, len(self.flux_list), 20, False)

        # 4 - finalise arrays
        self.dipole_moment_voltage_beta = np.array(self.dipole_moment_voltage_beta)
        self.spectrum_eigvals = np.array(self.spectrum_eigvals)
        self.spectrum_simulation_12 = np.array(self.spectrum_simulation_12)
        self.spectrum_simulation_23 = np.array(self.spectrum_simulation_23)

        logging.info("'simulate' finished")

    def sparse_matrix_visualise(self):
        """
        __ Parameters __

        __ Description __
        Visualise what the sparse matrix is going to look like with colours
        """
        plt.rc("text.latex", preamble=r"\usepackage{braket}")

        fig, ax = plt.subplots(nrows=1, ncols=1)
        plt.ion()

        # 1 - charing elements
        sparse_charge = sp.coo_matrix(
            (
                np.ones(len(self.op_H_charging_row)),
                (self.op_H_charging_row, self.op_H_charging_col),
            )
        ).tocsr()
        # 2 - diagonal elements
        diag_row = self.op_H_diag_row + self.op_H_diagA_row
        diag_col = self.op_H_diag_col + self.op_H_diagA_col
        sparse_diag = sp.coo_matrix(
            (np.ones(len(diag_row)), (diag_row, diag_col))
        ).tocsr()

        # 3 - phase dependent element
        phi_row = self.op_H_phi_row + self.op_H_phiAss_row
        phi_col = self.op_H_phi_col + self.op_H_phiAss_col
        sparse_phi = sp.coo_matrix((np.ones(len(phi_row)), (phi_row, phi_col))).tocsr()
        # 4 - plot spectrum
        ax.grid(b=True, which="major", color="black")
        ax.grid(b=True, which="minor")
        ax.spy(sparse_charge, color="C4", markersize=6)
        ax.spy(sparse_diag, color="C2", markersize=6)
        ax.spy(sparse_phi, markersize=6)
        ax.set_xlim([-0.5, self.states_total_number - 0.5])
        ax.set_ylim([self.states_total_number - 0.5, -0.5])

        # 5 - add the ticks
        total_ticks = self.states_total_number
        ax.set_xticks(
            np.linspace(0, self.states_total_number - 1, total_ticks), minor=True
        )
        ax.set_yticks(
            np.linspace(0, self.states_total_number - 1, total_ticks), minor=True
        )

        ax.set_xticks([0, 13, 26])
        xticklabels = [""] * 3
        xticklabels[0] = r"$\left|-1, -1, -1 \right\rangle$"
        xticklabels[1] = r"$\left|0, 0, 0 \right\rangle$"
        xticklabels[2] = r"$\left|+1, +1, +1 \right\rangle$"

        ax.set_yticks([0, 6, 13, 19, 26])
        yticklabels = [""] * 5
        yticklabels[1] = r"$\left\langle -1, +1, -1 \right|$"
        yticklabels[2] = r"$\left\langle 0, 0, 0 \right|$"
        yticklabels[3] = r"$\left\langle +1, -1, 0 \right|$"
        yticklabels[4] = r"$\left\langle +1, +1, +1 \right|$"

        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)
        utils.save_ree(ax, "output/fig4", "svg")
        plt.show()

    def track_progress(self, current_number, total_number, increment, heavy):
        """
        __ Parameters __
        current_number: current iteration
        total_number: full iteration
        increment: what % increments to plot in
        heavy: if false, a simple ouput is shown

        __ Description __
        print out current progress, with a given increment
        """
        no_stars = 50

        # 1 - find completed runs
        completion = current_number / total_number * 100

        # 2 - generate array is we hit increment
        if (int(completion * 1000) % int(increment * 1000)) == 0:
            if heavy:
                current_stars = int(completion / 100 * no_stars)
                stars = ["*"] * current_stars
                space = ["-"] * (no_stars - current_stars)
                stars.extend(space)
                output = "".join(stars)
                print("[%s][%i/%i]" % (output, current_number, total_number))
            else:
                print("  > [%i/%i]" % (current_number, total_number))

    def sparse_matrix_plot(self, sparse_matrix_to_plot):
        """
        __ Parameters __
        sparse_matrix_to_plot: matrix to show

        __ Description __
        show the sparse matrix via a plot
        """
        if self.plot_or_not:
            plt.spy(sparse_matrix_to_plot, markersize=3)
            plt.show()

    def plot_simulation(self, mpl_axes):
        """
        Plot the eigenvalues and transition spectrum
        """
        # 1 - prepare plot
        logging.info("Plotting results")

        if self.plot_or_not:
            mpl_axes.plot(
                self.flux_list,
                self.spectrum_simulation_12,
                label="1<->2",
                color="#004BA8",
                linewidth=1,
            )
            mpl_axes.plot(
                self.flux_list,
                self.spectrum_simulation_23,
                label="2<->3",
                color="C4",
                linewidth=1,
            )
            mpl_axes.set_ylim(0, 20)
            mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
            mpl_axes.set_ylabel("$\omega/2\pi$ (GHz)")

            plt.show()

    def plot_dipole_moment_voltage(self, mpl_axes):
        """
        __ Parameters __
        mpl_axes: where to output result to

         __ Description __
        Plot the absolute value of the dipole moment
        """

        # 1 - prepare plot
        logging.info("Plotting results")

        if self.plot_or_not:
            mpl_axes.plot(
                self.flux_list,
                (
                    self.dipole_moment_voltage[:, 0] ** 2
                    + self.dipole_moment_voltage[:, 1] ** 2
                )
                ** (1 / 2),
                label="1<->2",
                color="C6",
            )

            mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
            mpl_axes.set_ylabel(
                r"$\left|\left|\langle 1|\hat{V}_{20}|2 \rangle\right|\right|$ ($\mu$V)"
            )

    def plot_dipole_moment_voltage_beta(self, mpl_axes):
        """
        __ Parameters __
        mpl_axes: where to output result to

        __ Functionality __
        Plot beta evaluated from the dipole moment
        """

        # 1 - prepare plot
        logging.info("Plotting results")

        if self.plot_or_not:
            mpl_axes.plot(
                self.flux_list,
                (
                    self.dipole_moment_voltage_beta[:, 0] ** 2
                    + self.dipole_moment_voltage_beta[:, 1] ** 2
                )
                ** (1 / 2),
                label="1<->2",
                color="C6",
                dashes=[3, 3],
            )

            mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
            mpl_axes.set_ylabel(r"$\left|\left|\beta_{twin}\right|\right|$")
            plt.show()

    def experimental_data_load(self, mpl_axes, set_flux_list):
        """
        __ Parameters __
        mpl_axes: axes to plot the raw data on
        set_flux_list: True to set experimental data for the flux list

        __ Description __
        1 - load up experimental data
        2 - perform conversion from mA to flux
        3 - plot data
        4 - store the magnetic field points of the experiment

        Plots and finds differences between simulation and experiment
        """
        logging.info("'experimental_data_load' Importing data files")

        # 1 - data files to load
        base_data_name = "data/Qubit15_5JJ_Q2_"
        transition_12 = ["m3", "m2", "m1", "1", "2", "3"]
        transition_23 = ["m3b", "m2b", "m1b", "1b", "2b"]

        # 2 - for each data file
        self.flux_list_experimental_12 = []
        self.spectrum_experimental_12 = []
        self.flux_list_experimental_23 = []
        self.spectrum_experimental_23 = []

        for data_set in range(0, len(transition_12)):
            # a - generate file name
            temp_name = base_data_name + transition_12[data_set] + ".txt"

            # b - import, sort and convert to flux
            temp_data = self.import_and_sort_array(temp_name, 0)
            temp_data[0] = self.convert_mA_to_flux(temp_data[0], 0.125, 0.7)

            # c - plot data
            if self.plot_or_not:
                mpl_axes.plot(
                    temp_data[0],
                    temp_data[1],
                    marker="o",
                    color="#004BA8",
                    markeredgecolor="C2",
                    markersize=7,
                    alpha=0.95,
                    linestyle="",
                )

            # d - store imported data in 1 array
            self.flux_list_experimental_12.extend(temp_data[0])
            self.spectrum_experimental_12.extend(temp_data[1])

        for data_set in range(0, len(transition_23)):
            # a - generate file name
            temp_name = base_data_name + transition_23[data_set] + ".txt"

            # b - import, sort and convert to flux
            temp_data = self.import_and_sort_array(temp_name, 0)
            temp_data[0] = self.convert_mA_to_flux(temp_data[0], 0.125, 0.7)

            # c - plot data
            if self.plot_or_not:
                mpl_axes.plot(
                    temp_data[0],
                    temp_data[1],
                    marker="o",
                    color="C4",
                    markeredgecolor="#fb2c07",
                    markeredgewidth="0.4",
                    markersize=5,
                    alpha=0.95,
                    linestyle="",
                )

            # d - store imported fluxes
            self.flux_list_experimental_23.extend(temp_data[0])
            self.spectrum_experimental_23.extend(temp_data[1])

        # 3 - tidy up by sorting arrays
        temp_array = self.sort_array(
            np.array([self.flux_list_experimental_12, self.spectrum_experimental_12]), 0
        )
        self.flux_list_experimental_12 = temp_array[0]
        self.spectrum_experimental_12 = temp_array[1]

        temp_array = self.sort_array(
            np.array([self.flux_list_experimental_23, self.spectrum_experimental_23]), 0
        )
        self.flux_list_experimental_23 = temp_array[0]
        self.spectrum_experimental_23 = temp_array[1]

        logging.info(
            "Imported %i flux points",
            len(
                list(self.flux_list_experimental_12)
                + list(self.flux_list_experimental_23)
            ),
        )

        # 4 - set the flux array is required (then simulations are only done to compare with
        # experimental points)
        if set_flux_list:
            logging.info("Set experimental flux points for simulation")
            self.flux_list = np.array(
                list(self.flux_list_experimental_12)
                + list(self.flux_list_experimental_23)
            )
            self.flux_list.sort()

        if self.plot_or_not:
            plt.show()

        logging.info("'experimental_data_load' finished")

    def experimental_data_error(self):
        """
        __ Description __
        Error between the experimental data points and the simulation
        """
        try:
            logging.info(
                "==> 'experimental_data_error' comparing simulation to experiment"
            )

            if not hasattr(self, "spectrum_simulation_12"):
                raise RuntimeError(
                    "run the simulation first before calling 'experimental_data_error'"
                )
            if not hasattr(self, "spectrum_experimental_12"):
                raise RuntimeError(
                    "import experimental data before calling 'experimental_data_error'"
                )

            # 1 - prepare parameters
            error_cumulative = 0

            # 2 - transition12
            for i in range(0, len(self.flux_list_experimental_12)):
                # a - find the simulation entry corresponding to the experimental data point
                entry = np.where(self.flux_list == self.flux_list_experimental_12[i])[
                    0
                ][0]

                # b - compute the difference
                error_cumulative = (
                    error_cumulative
                    + (
                        self.spectrum_simulation_12[entry]
                        - self.spectrum_experimental_12[i]
                    )[0]
                    ** 2
                )

            # 3 - transition23
            for i in range(0, len(self.flux_list_experimental_23)):
                # a - find the simulation entry corresponding to the experimental data point
                entry = np.where(self.flux_list == self.flux_list_experimental_23[i])[
                    0
                ][0]

                # b - compute the difference
                error_cumulative = (
                    error_cumulative
                    + (
                        self.spectrum_simulation_23[entry]
                        - self.spectrum_experimental_23[i]
                    )[0]
                    ** 2
                )

            logging.info("'experimental_data_error' finished")

            return error_cumulative
        except TypeError as e:
            print(e)

    def import_and_sort_array(self, file_name, i):
        """
        __ Parameters __
        file_name: file to import, with path and extensions
        i: column to sort by

        __ Description __
        Import the 2D array and sort by the i-th column

        __ Return __
        return the sorted array
        """

        # 1 - import array
        array_to_return = np.loadtxt(file_name).transpose()

        # 2 - sort array by the i-th column
        sorted_index = np.argsort(array_to_return[i])
        array_to_return[0] = np.array(array_to_return[0])[sorted_index]
        array_to_return[1] = np.array(array_to_return[1])[sorted_index]

        return array_to_return

    def sort_array(self, array_to_sort, column_to_sort_by):
        """
        __ Parameters __
        array_to_sort: array to perform sorting for
        column_to_sort_by: which column to use for sorting

        __ Description __
        sorts the supplied array
        """
        # 1 - sort array by the i-th column
        sorted_index = np.argsort(array_to_sort[column_to_sort_by])
        for i in range(0, len(array_to_sort)):
            array_to_sort[i] = np.array(array_to_sort[i])[sorted_index]
        return array_to_sort

    def convert_index_to_island_state(self, index):
        """
        __ Parameters __
        index: number convert (0 <= index < states_total_number)

        __ Description __
        Takes an index and converts it to a unique cooper pair distribution
        across the 3 islands of the system - (n1, n2, n3).

        __ Changes __
        state_numerical_distribution:
            an ARRAY of the 3-island state
        state_cp_distribution:
            a MATRIX of the 3-island states, so that central index has a charge
            distribution (0,0,0)
        """

        # 0 - error checking
        if (index < 0) or (index >= self.states_total_number):
            print(
                "Index %data_set must be within 0 and %data_set -> converting to 0"
                % (index, self.states_total_number)
            )
            index = 0

        # 1 - evaluate the offset to apply, to get centering about (0,0,0)
        minimal_number_of_cp_on_island = (1 - self.states_total_number) / 2

        # 2 - represent 'index' by a trinary and cooper pair distribution
        self.state_numerical_distribution = []
        self.state_cp_distribution = []
        workingNumber = index

        # 3 - decomposition of a number to self.states_total_number base-system
        for island in [2, 1, 0]:

            # a - get the cp representation on this island
            island_no_cp = int(workingNumber / (self.states_total_number ** island))

            # b - store numerical and cp values
            self.state_numerical_distribution.append(island_no_cp)
            self.state_cp_distribution.append(
                island_no_cp + minimal_number_of_cp_on_island
            )

            # c - decrease the working number
            workingNumber = (
                workingNumber - island_no_cp * (self.states_total_number) ** island
            )

        # 4 - convert to numpy arrays
        self.state_numerical_distribution = np.array(self.state_numerical_distribution)
        self.state_cp_distribution = np.matrix(self.state_cp_distribution)

    def convert_numerical_state_to_index(self, state_numerical_distribution):
        """
        __ Parameters __
        state_numerical_distribution: array of the trinary state (NOT NORMALISED!)

        __ Description __
        Given a number in the 3 island representation
        convert it to an index for the matrix

        __ Return __
        index between 0 and total number of states that corresponds to supplied
        array
        """
        index = 0
        for data_set in range(0, 3):
            index = index + state_numerical_distribution[
                data_set
            ] * self.states_total_number ** (2 - data_set)

        return index

    def convert_mA_to_flux(
        self, mA_array: np.ndarray, offset: float, period: float
    ) -> np.ndarray:
        """
        __ Parameters __
        mA_array: array of mA values to convert
        offsset:  offset from 0
        period:   period in mA

        __ Description __
        Converts array measured in experimental mA to flux units
        """

        return (mA_array - offset) / period
