class TwinQubitOperatorBuilder(object):
    """Builds the voltage and phase operators as np.csr_matrix

    """

    def __init__(self):
        pass


        def generate_operators(self) -> Tuple[sp.csr_matrix, sp.csr_matrix]:
        """-> (op_V, op_Phi)"""

        op_V_dict = defaultdict(list)
        op_Phi_dict = defaultdict(list)
        op_V_constant = (
            self.const_h
            * 10 ** 9
            * self.EC
            / (2 * self.const_eCharge * (1 + self.alpha))
        )

        for x in range(0, self.states_total_number):

            # 2 - convert the index to island occupation
            (state_indx, state_cp) = self.convert_index_to_system_state(x)

            voltage_elm = op_V_constant * (np.dot(state_cp, [1, 2, 1]))
            op_V_dict["row-col"].append(x)
            op_V_dict["elm"].append(voltage_elm)

            # 4 - phase operator e^{i phi_20}
            if state_indx[1] < (self.states_per_island - 1):
                # island 2 (element 1)
                y = self.convert_system_state_to_index(state_indx + [0, 1, 0])
                op_Phi_dict["row"].append(x)
                op_Phi_dict["col"].append(y)
                op_Phi_dict["elm"].append(1)

        # 5 - finalise operators
        self.op_V = sp.coo_matrix(
            (op_V_dict["elm"], (op_V_dict["row-col"], op_V_dict["row-col"]))
        ).tocsr()
        self.op_Phi = sp.coo_matrix(
            (op_Phi_dict["elm"], (op_Phi_dict["row"], op_Phi_dict["col"]))
        ).tocsr()
