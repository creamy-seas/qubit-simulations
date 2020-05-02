class SparseMatrixVisualiser(object):
    """Visualise the sparse matrix"""

    def __init__(self):
        pass

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
