    def prepare_plot(self, nrows, ncols) -> Tuple[Figure, Axes]:

        plt.ioff()
        plt.close("all")
        fig = None
        mpl_axes = None

        if self.plot_or_not:
            logging.debug("==> 'prepare_plot' is setting up figure and axes")

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
