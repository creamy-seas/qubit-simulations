from twin import twin
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
from skcuda import linalg
plt.style.use('ilya_plot')


class twin_cuda(twin):
    """
    Adds GPU functionality to the twin class
    """

    def __init__(self, alpha, assymetry, states_per_island, flux_points, plot_or_not, message_or_not):
        """
        __ Parameters __
        EC:     charging energy
        EJ:     josephson energy
        alpha:  central island
        assymetry: assymetry between two loops
        states_per_island: accuracy of simulation
        plot_or_not: show plot output or not

        __ Description __
        Perform the general setup
        - measured in nm
        - working with frequencies (normalised by hbar)
        - working in unit of Phi0
        """

        twin.__init__(self, alpha, assymetry, states_per_island,
                      flux_points, plot_or_not, message_or_not)
        linalg.init()

    def simulate(self):
        """
        __ Parameters __

        __ Description __
        Method performs the eigenvalue simulations:
        1 - finish building the Hamiltonian, with the flux-dependent elements
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        3 - plot out the spectrum

        """

        # 0 - prepare hamiltonian for this simulation
        self.prepare_hamiltonian_stage2()

        if (self.message_or_not):
            print("==> 'simulate' running")

        self.spectrum_eigvals = []
        self.spectrum_eigvecs = []
        self.spectrum_simulation_12 = []
        self.spectrum_simulation_23 = []

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
            # a - copy base elements
            op_H_elm = self.op_H_elm.copy()

            # b - add on the phase dependent elements
            op_H_elm.extend(self.prepare_hamiltonian_stage3(
                phi_external, phi_externalAss))

            if((len(op_H_row) != len(op_H_elm)) or (len(op_H_col) != len(op_H_elm))):
                self.raise_error("Hamiltonin lists have %i rows, %i columns  and only %i elements" % (
                    len(op_H_row), len(op_H_col), len(op_H_elm)))

            # 2 - construct numpy matrix (using intermediate sparse matrix)
            self.op_H = sp.coo_matrix((op_H_elm,
                                       (op_H_row, op_H_col))).tocsr()
            gpu_array = self.op_H.toarray(order="F")
            gpu_array = gpuarray.to_gpu(gpu_array)

            # 3 - evaluate 3 lowest eigenenergies and eigenvectors, |1> |2> |3>
            eigvecs, eigvals = linalg.eig(gpu_array, jobvl='N', jobvr='V')
            self.spectrum_eigvals.append(eigvals)
            self.spectrum_simulation_12.append([eigvals[1] - eigvals[0]])
            self.spectrum_simulation_23.append([eigvals[2] - eigvals[1]])

            # self.track_progress(ext_flux_number, len(
            #     self.flux_list), 20, False)

        # 4 - finalise arrays
        self.spectrum_eigvals = np.array(self.spectrum_eigvals)
        self.spectrum_simulation_12 = np.array(self.spectrum_simulation_12)
        self.spectrum_simulation_23 = np.array(self.spectrum_simulation_23)

        if (self.message_or_not):
            print("==> 'simulate' finished")


if __name__ == "__main__":
    print("\nRunning 'twin_cuda.py'\n")
    start = time.time()
    EC = 31
    EJ = 22
    alpha = 1.023
    assymetry = 1.011
    test = twin_cuda(alpha, assymetry, 5, 300, False)
    test.override_parameters(EC, EJ, test.alpha, test.assymetry)
    test.experimental_data_load(test.ax, True)
    test.simulate()
    print(test.experimental_data_error())
    end = time.time()
    print(end - start)
