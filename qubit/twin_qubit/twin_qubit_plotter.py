class TwinQubitPlotter(object):
    def __init__(self):
        pass

        def prepare_plot(self, nrows, ncols) -> Tuple[Figure, Axes]:

        plt.ioff()
        plt.close("all")
        fig = None
        mpl_axes = None

        if self.plot_or_not:
            logging.info("==> 'prepare_plot' is setting up figure and axes")

            # Interactive mode, to alow updating
            plt.ion()

            fig, mpl_axes = plt.subplots(nrows=nrows, ncols=ncols)
            try:
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(0, 30, 1280, 1600)
                fig.canvas.set_window_title("Run time data")

            except AttributeError:
                pass

        return fig, mpl_axes

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
