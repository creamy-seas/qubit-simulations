import logging


class TwinQubit:
    """
    Quantum class to work with the twin, 5JJ qubit
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
        JJ_SQUARES = 2
        # --------------------
        self.delta = 0

        # 2 - very common parameters, which should only be set once
        self.verify_input()
        self.cp_offset_to_apply = (self.states_per_island - 1) // 2

        # 3 - simulation preparation
        self.set_derived_constants(JJ_SQUARES)
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

    def verify_input(self):
        """
        Correct common errors before program begins
        """
        if (self.states_total_number % 2) == 0:
            raise ValueError(
                f"Total number of states is {self.states_total_number} -> it needs to be odd"
            )

    def set_derived_constants(self, jj_squares: int):
        """Sets:
        EC
        EJ
        capacitance
        """
        logging.info("Setting energies and capacitance matrix")

        jj_overlap_area = 200 * 200 * jj_squares

        # Set the energies EC and EJ
        self.set_josephson_energy(jj_squares)
        self.set_charging_energy(jj_overlap_area)

        # Set capacitance matrix
        _capacitance = np.matrix(
            [[2, -1, 0], [-1, 2 + self.alpha + self.delta, -1], [0, -1, 2]]
        )
        self.capacitance = np.array(_capacitance)
        self.capacitance_inverse = np.array(_capacitance.I)
