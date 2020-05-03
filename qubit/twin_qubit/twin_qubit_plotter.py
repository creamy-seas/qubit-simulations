from typing import Tuple, List
import logging

from matplotlib.pyplot import Axes, Figure
import matplotlib.pyplot as plt
import pinject


class TwinQubitPlotter:
    @pinject.copy_args_to_public_fields
    def __init__(self, flux_list, twin_qubit_simulator):
        pass

    def plot_transitions(self, mpl_axes: Axes):

        logging.info("📈 Plotting transitions")

        mpl_axes.plot(
            self.flux_list,
            self.twin_qubit_simulator.simulations["1-2"],
            label="1<->2",
            color="#004BA8",
            linewidth=1,
        )
        mpl_axes.plot(
            self.flux_list,
            self.twin_qubit_simulator.simulations["2-3"],
            label="2<->3",
            color="C4",
            linewidth=1,
        )

        mpl_axes.set_ylim(0, 20)
        mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
        mpl_axes.set_ylabel("$\omega/2\pi$ (GHz)")

    def plot_dipole_moment_voltage(self, mpl_axes: Axes):

        logging.info("📈 Plotting dipole moment")

        mpl_axes.plot(
            self.flux_list,
            (
                self.twin_qubit_simulator.simulations["dipole-voltage"][:, 0] ** 2
                + self.twin_qubit_simulator.simulations["dipole-voltage"][:, 1] ** 2
            )
            ** (1 / 2),
            label="1<->2",
            color="C6",
        )

        mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
        mpl_axes.set_ylabel(
            r"$\left|\left|\langle 1|\hat{V}_{20}|2 \rangle\right|\right|$ ($\mu$V)"
        )

    def plot_dipole_moment_voltage_beta(self, mpl_axes: Axes):

        logging.info("📈 Plotting dipole moment beta")

        mpl_axes.plot(
            self.flux_list,
            (
                self.twin_qubit_simulator.simulations["dipole-beta"][:, 0] ** 2
                + self.twin_qubit_simulator.simulations["dipole-beta"][:, 1] ** 2
            )
            ** (1 / 2),
            label="1<->2",
            color="C6",
        )

        mpl_axes.set_xlabel("Magnetic Flux ($\Phi$)")
        mpl_axes.set_ylabel(r"$\left|\left|\beta_{twin}\right|\right|$")
