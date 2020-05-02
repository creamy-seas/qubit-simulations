import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from quantum_master import quantum_master
import honkler
import time
plt.style.use('ilya_plot')


class flux(quantum_master):
    """
    Quantum class to work with the twin, 5JJ qubit.
    """

    def __init__(self, alpha, states_per_island, flux_points, plot_or_not, message_or_not):
        """
        __ Parameters __
        EC:     charging energy
        EJ:     josephson energy
        alpha:  central island
        states_per_island: accuracy of simulation
        plot_or_not: show plot output or not
        message_or_not: display output messages or not

        __ Description __
        Perform the general setup
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """
        quantum_master.__init__(self, plot_or_not, message_or_not)

        # 1 - store user supplied parameters
        self.alpha = alpha
        self.delta = 0

        # 2 - global parameters
        self.flux_min = -3
        self.flux_max = 3
        self.flux_points = flux_points
        self.states_per_island = states_per_island
        self.states_total_number = self.states_per_island**3
        self.flux_list = np.linspace(
            self.flux_min, self.flux_max, self.flux_points)

        # 3 - simulation preparation
        self.prepare_correction()
        self.prepare_structure()
        self.prepare_hamiltonian_stage1()

        # 4 - plot preparation
        self.fig, self.ax = self.prepare_plot(1, 2)

    def override_parameters(self, EC, EJ, alpha):
        """
        __ Parameters __
        EC, EJ, alpha: parameters to set

        __ Description __
        For new simulations, set the new system parameters
        """
        if (self.message_or_not):
            print("==> 'override_parameters' with EC=%.4f\tEJ=%.4f\talpha=%.4f" %
                  (EC, EJ, alpha))

        self.EC = EC
        self.EJ = EJ
        self.alpha = alpha

    def prepare_correction(self):
        """
        Correct common errors before program begins
        """
        if ((self.states_per_island % 2) == 0):
            print("==> Changing states_per_island from %data_set -> %data_set" %
                  (self.states_per_island, self.states_per_island - 1))
            self.states_per_island = self.states_per_island - 1

    def prepare_structure(self):
        """
        Prepares the structure parameters
        """
        if (self.message_or_not):
            print("==> 'prepare_structure' creating energies and capacitances")

        # 1 - set jj dimensions
        self.param_jj_squares = 2
        self.param_jj_overlap_area = 200 * 200 * \
            self.param_jj_squares

        # 2 - set the energies EC and EJ
        self.energy_charging(self.param_jj_overlap_area)
        self.energy_josephson(self.param_jj_squares)

        # 3 - set capacitance matrix
        self.capacitance_normalised = np.matrix([
            [2, -1, 0],
            [-1, 2, -1],
            [0, -1, 1 + self.alpha]])

    def prepare_operators(self):
        """
        __ Description __
        Prepare operators for evaluating expectation values once the eigenstates are found

        __ Return __
        self.op_V
        self.op_Phi
        """

        # additional operator components
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
            voltage_element = self.EC * self.const_h * 10**9 / ((1 + 3 * self.alpha) * self.const_eCharge) * (
                self.state_cp_distribution[0, 0] +
                2 * self.state_cp_distribution[0, 1] +
                3 * self.state_cp_distribution[0, 2])
            self.op_V_row.append(x)
            self.op_V_col.append(x)
            self.op_V_elm.append(voltage_element)

            # 4 - phase operator phi_20
            if (self.state_numerical_distribution[1] < (self.states_per_island - 1)):
                # island 2 (element 1)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [0, 1, 0])

                self.op_Phi_row.append(x)
                self.op_Phi_col.append(y)
                self.op_Phi_elements.append(1)

        # 5 - finalise operators
        self.op_V = sp.coo_matrix((self.op_V_elm,
                                   (self.op_V_row, self.op_V_col))).tocsr()
        self.op_Phi = sp.coo_matrix((self.op_Phi_elements,
                                     (self.op_Phi_row, self.op_Phi_col))).tocsr()

    def prepare_hamiltonian_stage1(self):
        """
        Generate elements for a normalised Hamiltonian.

        During the simulation 'prepare_hamiltonian()' must be run, to change energies
        """
        if (self.message_or_not):
            print("==> 'prepare_normalised_hamiltonian' creating Hamiltonian entries")

        # 0 - Individual Hamiltonian components.
        # constant
        self.op_H_charging_elm = []
        self.op_H_charging_row = []
        self.op_H_charging_col = []
        self.op_H_diag_row = []
        self.op_H_diag_col = []
        # changing
        self.op_H_phi_row = []
        self.op_H_phi_col = []

        # 1 - generate matrix elements, by going across all the states
        for x in range(0, self.states_total_number):

            # 2 - convert the index to island occupation
            self.convert_index_to_island_state(x)

            # 3 - diagonal charging energy
            charging_energy = (self.state_cp_distribution *
                               (self.capacitance_normalised.I) *
                               (self.state_cp_distribution.transpose()))

            self.op_H_charging_row.append(x)
            self.op_H_charging_col.append(x)
            self.op_H_charging_elm.append(charging_energy[0, 0])

            # 4 - offdiagonal JJ energy - check that within matrix boundaries
            # offset y coordinate from diagonal, and fill out the symmetrical entries
            # ! 'diag_elm' array is not created, since all elements will be the same
            if (self.state_numerical_distribution[0] < (self.states_per_island - 1)):
                # cos(phi_10)
                y = self.convert_numerical_state_to_index(
                    self.state_numerical_distribution + [1, 0, 0])
                self.op_H_diag_row.extend([x, y])
                self.op_H_diag_col.extend([y, x])

                # cos(phi_21) cp exchange between 1 <-> 2
                if (self.state_numerical_distribution[1] > 0):
                    y = self.convert_numerical_state_to_index(
                        self.state_numerical_distribution + [1, -1, 0])
                    self.op_H_phi_row.extend([x, y])
                    self.op_H_phi_col.extend([y, x])

            if (self.state_numerical_distribution[1] < (self.states_per_island - 1)):
                if(self.state_numerical_distribution[0] < (self.states_per_island - 1)):
                    if(self.state_numerical_distribution[2] < (self.states_per_island - 1)):
                        # cos(phi_ext - phi_1 - phi_2 - phi_3)
                        y = self.convert_numerical_state_to_index(
                            self.state_numerical_distribution + [1, 1, 1])
                        self.op_H_diag_row.extend([x, y])
                        self.op_H_diag_col.extend([y, x])

            if (self.state_numerical_distribution[2] < (self.states_per_island - 1)):
                # cos(phi_32) cp exchange between 2 <-> 3
                if (self.state_numerical_distribution[1] > 0):
                    y = self.convert_numerical_state_to_index(
                        self.state_numerical_distribution + [0, -1, 1])
                    self.op_H_diag_row.extend([x, y])
                    self.op_H_diag_col.extend([y, x])

        if (self.message_or_not):
            print("  > Unchaning part of Hamiltonian has %i entries" %
                  (len(self.op_H_charging_row + self.op_H_diag_row)))
            print("  > Flux-dependent part of Hamiltonian has %i entries" %
                  (len(self.op_H_phi_row)))
            print("==> 'prepare_normalised_hamiltonian' finished")

    def prepare_hamiltonian_stage2(self):
        """
        __ Description __
        Using supplied EC, EJ, alpha, scale the normalised list created in 'prepare_normalised_hamiltonian' and
        combine them into 1.

        The lists are used in 'simulate' functions
        """

        if (self.message_or_not):
            print("==> 'prepare_hamiltonian' with EC=%.4f\tEJ=%.4f\talpha=%.4f" %
                  (self.EC, self.EJ, self.alpha))

        # 1 - main part of the Hamiltonian, which remains unchanged during simulation
        temp_charging_elm = self.EC * np.array(self.op_H_charging_elm)
        temp_diag_elm = - self.EJ / 2 * np.ones(len(self.op_H_diag_row))

        self.op_H_elm = list(temp_charging_elm) + list(temp_diag_elm)
        self.op_H_row = self.op_H_charging_row + self.op_H_diag_row
        self.op_H_col = self.op_H_charging_col + self.op_H_diag_col

        if((len(self.op_H_row) != len(self.op_H_elm)) or (len(self.op_H_col) != len(self.op_H_elm))):
            self.raise_error("Hamiltonin has %i rows, %i columns non zero coordinatres and only %i elements" % (
                len(self.op_H_row), len(self.op_H_col), len(self.op_H_elm)))

        # 2 - exchange elements - used by 'prepare_hamiltonian_stage3' to fill out different fluxes
        no_exchange = int(len(self.op_H_phi_row) / 2)
        self.op_H_SUPPORT_exchange_elm = np.ones(no_exchange)

        if (self.message_or_not):
            print("==> 'prepare_hamiltonian' finished")

    def prepare_hamiltonian_stage3(self, phi_external):
        """
        __ Parameters __
        phi_external: fluxes to the two loops

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

        temp_phi_p = -(self.EJ / 2 * self.alpha * np.exp(1j * phi_external)) * \
            self.op_H_SUPPORT_exchange_elm
        temp_phi_n = -(self.EJ / 2 * self.alpha * np.exp(-1j * phi_external)) * \
            self.op_H_SUPPORT_exchange_elm

        # 2 - interleave the lists, alternating between +ve and -ve fluxes, for the (x,y) and (y,x)
        # pairings in the Hamiltonian
        elm_list = list(np.ravel(np.column_stack((temp_phi_p, temp_phi_n))))

        return elm_list

    def simulate(self, simulate_expectations=[False, False]):
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

        if (self.message_or_not):
            print("==> 'simulate' running")

        self.spectrum_eigvals = []
        self.spectrum_eigvecs = []
        self.spectrum_simulation_12 = []
        self.spectrum_simulation_23 = []
        self.dipole_moment_voltage = []
        self.dipole_moment_voltage_beta = []
        self.dipole_moment_phi = []

        # copy main parts of Hamiltonian
        op_H_row = self.op_H_row.copy()
        op_H_col = self.op_H_col.copy()
        op_H_row.extend(self.op_H_phi_row)
        op_H_col.extend(self.op_H_phi_col)

        # 1 - simulation performed for each value of the external flux
        for ext_flux_number in range(0, len(self.flux_list)):

            # 2 - extraction of the onset phase, due to external flux
            phi_external = (self.flux_list[ext_flux_number]) * np.pi

            ####################
            # 3 - FINISH MAIN HAMILTONIAN for this particular bias
            ####################
            # a - copy constant terms
            op_H_elm = self.op_H_elm.copy()

            # b - add on the phase dependent elements
            op_H_elm.extend(self.prepare_hamiltonian_stage3(phi_external))

            if((len(op_H_row) != len(op_H_elm)) or (len(op_H_col) != len(op_H_elm))):
                self.raise_error("Hamiltonin lists have %i rows, %i columns  and only %i elements" % (
                    len(op_H_row), len(op_H_col), len(op_H_elm)))

            # 2 - construct sparse matrix
            self.op_H = sp.coo_matrix((op_H_elm,
                                       (op_H_row, op_H_col))).tocsr()

            # 3 - evaluate 3 lowest eigenenergies and eigenvectors, |1> |2> |3>
            eigvals, eigvecs = eigsh(self.op_H, 3, which='SA', tol=0)
            # sort in ascending order
            sort_index = np.argsort(eigvals)
            eigvals = np.array(eigvals)[sort_index]
            eigvecs = np.transpose(np.transpose(eigvecs)[sort_index])

            self.spectrum_eigvals.append(eigvals)
            self.spectrum_simulation_12.append([eigvals[1] - eigvals[0]])
            self.spectrum_simulation_23.append([eigvals[2] - eigvals[1]])

            # 4 - dipole transition
            if(simulate_voltage):
                # evalute voltage in a sandwhich between the ground and excited states
                state_0 = eigvecs[:, 0]
                state_1 = eigvecs[:, 1]
                temp_dipole_moment_voltage = state_0.dot(
                    self.op_V.dot(state_1))
                temp_dipole_moment_voltage_beta = temp_dipole_moment_voltage / \
                    self.const_Phi0 / \
                    self.spectrum_simulation_12[ext_flux_number]

                self.dipole_moment_voltage.append(
                    [temp_dipole_moment_voltage.real, temp_dipole_moment_voltage.imag])
                self.dipole_moment_voltage_beta.append(
                    [temp_dipole_moment_voltage_beta.real, temp_dipole_moment_voltage_beta.imag])

            self.track_progress(ext_flux_number, len(
                self.flux_list), 20, False)

            # time.sleep(0.005)

        # 4 - finalise arrays
        self.dipole_moment_voltage = np.array(self.dipole_moment_voltage)
        self.dipole_moment_voltage_beta = np.array(
            self.dipole_moment_voltage_beta)
        self.spectrum_eigvals = np.array(self.spectrum_eigvals)
        self.spectrum_simulation_12 = np.array(self.spectrum_simulation_12)
        self.spectrum_simulation_23 = np.array(self.spectrum_simulation_23)

        if (self.message_or_not):
            print("==> 'simulate' finished")

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
        if((int(completion * 1000) % int(increment * 1000)) == 0):
            if(heavy):
                current_stars = int(completion / 100 * no_stars)
                stars = ["*"] * current_stars
                space = ["-"] * (no_stars - current_stars)
                stars.extend(space)
                output = "".join(stars)
                print("[%s][%i/%i]" %
                      (output, current_number, total_number))
            else:
                print("  > [%i/%i]" %
                      (current_number, total_number))

    def plot_simulation(self, plotAxes):
        """
        Plot the eigenvalues and transition spectrum
        """
        # 1 - prepare plot
        if (self.message_or_not):
            print("==> Plotting results")

        if(self.plot_or_not):
            plotAxes.plot(self.flux_list,
                          self.spectrum_simulation_12, label="1<->2", color='#004BA8')
            plotAxes.plot(self.flux_list,
                          self.spectrum_simulation_23, label="2<->3", color='C4')
            plotAxes.set_ylim(0, 20)
            plotAxes.set_xlabel("Magnetic Flux ($\Phi$)")
            plotAxes.set_ylabel("$\omega/2\pi$ (GHz)")

            plt.show()

    def plot_dipole_moment_voltage(self, plotAxes):
        """
        __ Parameters __
        plotAxes: where to output result to

        __ Description __
        Plot the absolute value of the dipole moment
        """

        # 1 - prepare plot
        if (self.message_or_not):
            print("==> Plotting results")

        if(self.plot_or_not):
            plotAxes.plot(self.flux_list,
                          (self.dipole_moment_voltage[:, 0]**2 +
                           self.dipole_moment_voltage[:, 1]**2)**(1 / 2),
                          label="1<->2", color='C6')

            plotAxes.set_xlabel("Magnetic Flux ($\Phi$)")
            plotAxes.set_ylabel(
                r"$\left|\left|\langle 1|\hat{V}_{20}|2 \rangle\right|\right|$ ($\mu$V)")

            # real part
            # plotAxes.plot(self.flux_list,
            # self.dipole_moment_voltage[:, 0], label="1<->2re", color='C6')

            # imaginary part
            # axImag = plotFig.add_subplot(111, sharex=plotAxes, frameon=False)
            # axImag.yaxis.tick_right()
            # axImag.yaxis.set_label_position("right")
            # axImag.set_ylabel(
            #     r"$|\langle1|\hat{V}_{20}|2\rangle|^2$ ($\mu$V)")
            # axImag.plot(self.flux_list,
            #             (self.dipole_moment_voltage[:, 0]**2 +
            #              self.dipole_moment_voltage[:, 1]**2)**(1 / 2),
            #             label="1<->2im", color='C8')
            plt.show()

    def plot_dipole_moment_voltage_beta(self, plotAxes):
        """
        __ Parameters __
        plotAxes: where to output result to

        __ Functionality __
        Plot beta evaluated from the dipole moment
        """

        # 1 - prepare plot
        if (self.message_or_not):
            print("==> Plotting results")

        if(self.plot_or_not):
            plotAxes.plot(self.flux_list,
                          (self.dipole_moment_voltage_beta[:, 0]**2 +
                           self.dipole_moment_voltage_beta[:, 1]**2)**(1 / 2),
                          label="1<->2", color='C9')

            plotAxes.set_xlabel("Magnetic Flux ($\Phi$)")
            plotAxes.set_ylabel(
                r"$\left|\left|\beta_{4jj}\right|\right|$")
            plt.show()

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
        if((index < 0) or (index >= self.states_total_number)):
            print("Index %data_set must be within 0 and %data_set -> converting to 0" %
                  (index, self.states_total_number))
            index = 0

        # 1 - evaluate the offset to apply, to get centering about (0,0,0)
        minimal_number_of_cp_on_island = (1 - self.states_per_island) / 2

        # 2 - represent 'index' by a trinary and cooper pair distribution
        self.state_numerical_distribution = []
        self.state_cp_distribution = []
        workingNumber = index

        # 3 - decomposition of a number to self.states_per_island base-system
        for island in [2, 1, 0]:

            # a - get the cp representation on this island
            island_no_cp = int(
                workingNumber / (self.states_per_island**island))

            # b - store numerical and cp values
            self.state_numerical_distribution.append(island_no_cp)
            self.state_cp_distribution.append(
                island_no_cp + minimal_number_of_cp_on_island)

            # c - decrease the working number
            workingNumber = workingNumber - island_no_cp * \
                (self.states_per_island)**island

        # 4 - convert to numpy arrays
        self.state_numerical_distribution = np.array(
            self.state_numerical_distribution)
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
            index = index + \
                state_numerical_distribution[data_set] * \
                self.states_per_island**(2 - data_set)

        return index


if (__name__ == "__main__"):
    print("\nRunning 'twin.py'\n")
    start = time.time()

    # ############################################################
    # ################### Simulation  ############################
    ############################################################
    honkler.config_plot_size(0.05, 0.9, 0.15, 0.9)
    EC = 38
    EJ = 40
    alpha = 0.45
    test = flux(alpha, 7, 300, True, True)
    test.prepare_operators()
    test.override_parameters(EC, EJ, alpha)
    test.simulate([True, False])
    test.plot_simulation(test.ax[0])
    test.plot_dipole_moment_voltage_beta(test.ax[1])

    end = time.time()
    print("Time taken:\t%.3fs" % (end - start))
