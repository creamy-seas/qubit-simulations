class TwinQubitConverter(object):
    """Class that runs common conversions

    """
    def __init__(self):
        pass

        def convert_index_to_system_state(self, index) -> Tuple[np.ndarray, np.matrix]:
        """
        __ Parameters __
        [int] index:                    to convert into a one-to-one system state
                                        (0 <= index < states_total_number)

        [bool] set_class_variable:      set class variables to evaluate values?

        __ Description __
        Takes an index and converts it to a unique cooper pair distribution
        across the 3 islands of the system - (n1, n2, n3).

        __ Returns __
        [int, int, int] state_indx:            3-island index distribution
        [int, int, int] state_cp:              3-island cp distribution, centered about (0,0,0)
        """

        if (index < 0) or (index >= self.states_total_number):
            raise ValueError(
                "Index {index} must be within 0 and {self.states_total_number} -> converting to 0"
            )

        # 1 - convert index to  represetnation on the 3 islands
        state_indx = np.zeros(3, dtype=int)
        for idx_i, i in enumerate(
            np.base_repr(index, base=self.states_per_island)[::-1]
        ):
            state_indx[idx_i] = int(i)

        # 2 - generate cp distribution, centered around (0, 0, 0)
        state_cp = state_indx - self.cp_offset_to_apply

        return state_indx, state_cp

    def convert_system_state_to_index(self, state_indx: List) -> int:
        """
        __ Parameters __
        [int, int, int] -> inique index from 0 to states_per_island^3
        """

        return int("".join(map(str, state_indx)), base=self.states_per_island)

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

    def convert_energy_to_GHz(self, energy):
        return energy / (self.const_h * 10 ** (9))


    def set_charging_energy(self, area_nm2):
        """
        __ Parameters __
        area: JJ area in nm^2

        __ Description __
        Evaluating charging energy and capacitance
            (2e)^2/2C
            C = e*e_0*A/d

        __ Sets __
        self.EC:                the charging energy
        self.capacitance:       capacitance for normalisation
        """

        THICKNESS_ALOX = 2 * 10 ** (-9)
        PERMITIVITY_ALOX = 10
        PERMITIVITY_VACUUM = 8.85 * 10 ** (-12)
        CAPACITANCE_AREA_OFFSET = 1.15 * 10 ** (-18)  # convert nm2 to m2

        # 2 - evaluate and return
        self.param_capacitance = (
            PERMITIVITY_ALOX
            * PERMITIVITY_VACUUM
            * CAPACITANCE_AREA_OFFSET
            * area_nm2
            / THICKNESS_ALOX
        )
        self.EC = (2 * self.const_eCharge) ** 2 / (2 * self.param_capacitance)
        self.EC = self.convert_energy_to_GHz(self.EC)

    def set_josephson_energy(self, squares_JJ):
        """
        __ Parameters __
        squares: width of the JJ in square units

        __ Description __
        Evaluate the Josephson energy, critical current, and resistance.
        The more squares of JJ, the bigger the resistance, and the larger it's energy

        We assume that it is a standard JJ: 20nm - AlOx - 30nm
        Count in 100x100nm^2 squares

        __ Sets __
        self.EJ:                josephson energy in GHz
        self.critical_current:  critical current in A
        """

        TRANSITION_AL = 1.2
        BOLTZMAN = 1.38 * 10 ** (-23)
        DELTA_AL = 1.764 * TRANSITION_AL * BOLTZMAN
        RESISTANCE_CONSTANT = self.const_h / (4 * self.const_eCharge ** 2)
        RESISTANCE_OF_100x100_AL = 18.4 * 10 ** 3

        self.param_resistance = RESISTANCE_OF_100x100_AL / squares_JJ
        self.param_critical_current = (
            np.pi * DELTA_AL / (2 * self.const_eCharge * self.param_resistance)
        )
        self.EJ = RESISTANCE_CONSTANT * DELTA_AL / (2 * self.param_resistance)
        self.EJ = self.convert_energy_to_GHz(self.EJ)

    self.const_h = 6.64 * 10 ** (-34)
        self.const_eCharge = 1.6 * 10 ** (-19)
        self.const_Phi0 = self.const_h / (2 * self.const_eCharge)
