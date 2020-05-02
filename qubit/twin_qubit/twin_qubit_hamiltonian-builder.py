class TwinQubitHamiltonianBuilder(object):
    """Builds and stores a custom Hamiltonian for the twin qubit
    """

    def __init__(self):
        pass

    def prepare_hamiltonian_stage1(self):
        """
        __ Description __
        Generates building blocks of the Hamitltonian that need to be scaled later on
        |------------+------------+------------+---+-------------------------------|
        |            |            |            |   | Scaling                       |
        |------------+------------+------------+---+-------------------------------|
        | charge-elm | charge-row | charge-col |   | EC                            |
        | phi1-elm   | phi1-row   | phi1-col   |   | -EJ/2                         |
        | phi2-elm   | phi2-row   | phi2-col   |   | -alpha*EJ/2                   |
        | phi3-elm   | phi3-row   | phi3-col   |   | -EJ/2                         |
        | +phi21-elm | +phi21-row | +phi21-col |   | -EJ/2 * e^(+i * phi_ext)      |
        | -phi21-elm | -phi21-row | -phi21-col |   | -EJ/2 * e^(-i * phi_ext)      |
        | +phi32-elm | +phi32-row | +phi32-col |   | -EJ/2 * e^(-i * phi_ext * nu) |
        | -phi32-elm | -phi32-row | -phi32-col |   | -EJ/2 * e^(+i * phi_ext * nu) |
        |------------+------------+------------+---+-------------------------------|
        """

        logging.info("Creating hamiltonian for energies")

        self.op_H_build = defaultdict(lambda: defaultdict(list))
        self.op_H_build["fields"] = [
            "charge",
            "phi1",
            "phi2",
            "phi3",
            "+phi21",
            "-phi21",
            "+phi32",
            "-phi32",
        ]

        for x in range(0, self.states_total_number):

            (state_indx, state_cp) = self.convert_index_to_system_state(x)

            charging_energy = np.dot(
                state_cp, np.dot(self.capacitance_inverse, state_cp)
            )
            self.op_H_build["charge"]["row"].append(x)
            self.op_H_build["charge"]["col"].append(x)
            self.op_H_build["charge"]["elm"].append(charging_energy)

            # 4 - offdiagonal JJ energy - check within matrix boundaries
            # offset y coordinate from diagonal, and fill out the symmetrical entries
            if state_indx[0] < (self.states_per_island - 1):
                # cos(phi_1)
                y = self.convert_system_state_to_index(state_indx + [1, 0, 0])
                self.op_H_build["phi1"]["row"] += [x, y]
                self.op_H_build["phi1"]["col"] += [y, x]

                # cos(phi_2 - phi_1 - phi_ext), with cp exchange between 2 <-> 1
                if state_indx[1] > 0:
                    y = self.convert_system_state_to_index(state_indx + [1, -1, 0])
                    self.op_H_build["+phi21"]["row"] += [x]
                    self.op_H_build["+phi21"]["col"] += [y]
                    self.op_H_build["-phi21"]["row"] += [y]
                    self.op_H_build["-phi21"]["col"] += [x]

            if state_indx[1] < (self.states_per_island - 1):
                # alpha * cos(phi_2)
                y = self.convert_system_state_to_index(state_indx + [0, 1, 0])
                self.op_H_build["phi2"]["row"] += [x, y]
                self.op_H_build["phi2"]["col"] += [y, x]

            if state_indx[2] < (self.states_per_island - 1):
                # cos(phi_3)
                y = self.convert_system_state_to_index(state_indx + [0, 0, 1])
                self.op_H_build["phi3"]["row"] += [x, y]
                self.op_H_build["phi3"]["col"] += [y, x]

                # cos(phi_2 - phi_3 + nphi_ext), with cp exchange between 2 <-> 3
                if state_indx[1] > 0:
                    y = self.convert_system_state_to_index(state_indx + [0, -1, 1])
                    self.op_H_build["-phi32"]["row"] += [x]
                    self.op_H_build["-phi32"]["col"] += [y]
                    self.op_H_build["+phi32"]["row"] += [y]
                    self.op_H_build["+phi32"]["col"] += [x]

        # 5 - create arrays of elements
        self.op_H_build["charge"]["elm"] = np.array(self.op_H_build["charge"]["elm"])
        for i in self.op_H_build["fields"][1:]:
            self.op_H_build[i]["elm"] = np.ones(len(self.op_H_build[i]["row"]))

    def prepare_hamiltonian_stage2(self):
        """
        __ Description __
        - Scale the raw Hamiltonian elements by EC, EJ, alpha that are constant throughout
        simulations
        |------------+------------+------------+---+-------------------------------|
        |            |            |            |   | Scaling                       |
        |------------+------------+------------+---+-------------------------------|
        | charge-elm | charge-row | charge-col |   | EC                            |
        | phi1-elm   | phi1-row   | phi1-col   |   | -EJ/2                         |
        | phi2-elm   | phi2-row   | phi2-col   |   | -alpha * EJ/2                 |
        | phi3-elm   | phi3-row   | phi3-col   |   | -EJ/2                         |
        | ✘          | +phi21-row | +phi21-col |   |            ✘		   |
        | ✘          | -phi21-row | -phi21-col |   |            ✘		   |
        | ✘          | +phi32-row | +phi32-col |   |            ✘		   |
        | ✘          | -phi32-row | -phi32-col |   |            ✘		   |
        |------------+------------+------------+---+-------------------------------|
        | ['elm']    | ['row']    | ['col']    |   |                               |
        |------------+------------+------------+---+-------------------------------|
        """

        if self.message_or_not:
            print(
                f"⦿ 'prepare_hamiltonian_stage2' with EC={self.EC:.4f}\tEJ={self.EJ:.4f}\talpha={self.alpha:.4f}\tass={self.assymetry:.4f}"
            )

        self.op_H_simulate = defaultdict(partial(np.ndarray, 0))

        # 1 - elements are scaled by EC, -EJ/2, -alpha * EJ/2
        self.op_H_simulate["elm-const"] = np.hstack(
            (
                self.EC * self.op_H_build["charge"]["elm"],
                -self.EJ / 2 * self.op_H_build["phi1"]["elm"],
                -self.alpha * self.EJ / 2 * self.op_H_build["phi2"]["elm"],
                -self.EJ / 2 * self.op_H_build["phi3"]["elm"],
            )
        )

        # 2 - rows and cols are simply combined
        self.op_H_simulate["row"] = []
        self.op_H_simulate["col"] = []
        for i in self.op_H_build["fields"]:
            self.op_H_simulate["row"] += self.op_H_build[i]["row"]
            self.op_H_simulate["col"] += self.op_H_build[i]["col"]

        # 3 - error checks
        # 3a - number of elemets
        check = reduce(
            operator.sub,
            [
                len(self.op_H_simulate["elm-const"]),
                len(self.op_H_build["charge"]["elm"]),
                len(self.op_H_build["phi1"]["elm"]),
                len(self.op_H_build["phi2"]["elm"]),
                len(self.op_H_build["phi3"]["elm"]),
            ],
        )
        if check:
            self.raise_error(
                f"Inconsistency of {check} elements (op_H_simulate['elm] has length {len(self.op_H_simulate['elm-const'])})"
            )

        # 3b - row and colm numbers
        if len(self.op_H_simulate["row"]) != len(self.op_H_simulate["col"]):
            self.raise_error(
                f"Hamiltonian has {len(self.op_H_simulate['row'])} rows, {len(self.op_H_simulate['col'])} columns and only {len(self.op_H_simulate['elm-const'])} elements"
            )

        if self.message_or_not:
            print(
                f"  > Hamiltonian has:\n  {len(self.op_H_simulate['row'])} rows,\n  {len(self.op_H_simulate['col'])} columns\n  {len(self.op_H_simulate['elm-const'])} constant-elements"
            )
            print("⦿ 'prepare_hamiltonian_stage2' finished")

    def prepare_hamiltonian_stage3(self, phi_external, phi_externalAss):
        """
        __ Parameters __
        [float] phi_external, phi_externalAss:  fluxes to the two loops

        __ Description __
        Scale the flux-dependent part of the Hamiltonian changes as the field is being swept.
        |------------+------------+------------+---+-------------------------------|
        | +phi21-elm |            |            |   | -EJ/2 * e^(+i * phi_ext)      |
        | -phi21-elm |            |            |   | -EJ/2 * e^(-i * phi_ext)      |
        | +phi32-elm |            |            |   | -EJ/2 * e^(-i * phi_ext * nu) |
        | -phi32-elm |            |            |   | -EJ/2 * e^(+i * phi_ext * nu) |
        |------------+------------+------------+---+-------------------------------|
        """

        self.op_H_simulate["elm-var"] = np.hstack(
            (
                -self.EJ
                / 2
                * np.exp(1j * phi_external)
                * self.op_H_build["+phi21"]["elm"],
                -self.EJ
                / 2
                * np.exp(-1j * phi_external)
                * self.op_H_build["+phi21"]["elm"],
                -self.EJ
                / 2
                * np.exp(-1j * phi_externalAss)
                * self.op_H_build["+phi32"]["elm"],
                -self.EJ
                / 2
                * np.exp(1j * phi_externalAss)
                * self.op_H_build["-phi32"]["elm"],
            )
        )

    def simulate(self, simulate_dipole=False):
        """
        __ Parameters __
        [bool] simulate_dipole:   expectation values for the dipole

        __ Description __
        Method performs the eigenvalue simulations:
        1 - sweep fields
        2 - evaluate the lowest eigenstate and eigenenergies of the Hamiltonian
        3 - plot out the spectrum

        """

        if self.message_or_not:
            print("⦿ 'simulate' running")

        ########################################
        # ⦿ Build Hamiltonian and lists to fill out
        ########################################
        self.prepare_hamiltonian_stage2()
        self.simulation = defaultdict(list)
        self.simulation["eigvals"] = np.empty((0, 3))
        self.simulation["eigvecs"] = np.empty((0, self.states_total_number))

        ########################################
        # ⦿ 1 - simulation performed for each value of the external flux
        ########################################
        bar = pyprind.ProgBar(len(self.flux_list), bar_char="█")
        for idx, ext_flux_number in enumerate(self.flux_list):

            ########################################
            # ⦿ 2 - extraction of the onset phase, due to external flux
            ########################################
            phi_external = ext_flux_number * 2 * np.pi
            phi_externalAss = phi_external * self.assymetry

            ########################################
            # ⦿ 3 - build Hamiltonian for this particular flux
            ########################################
            self.prepare_hamiltonian_stage3(phi_external, phi_externalAss)
            self.op_H_simulate["elm"] = np.hstack(
                (self.op_H_simulate["elm-const"], self.op_H_simulate["elm-var"])
            )

            ########################################
            # ⦿ 4 - construct sparse matrix
            ########################################
            self.op_H = sp.coo_matrix(
                (
                    self.op_H_simulate["elm"],
                    (self.op_H_simulate["row"], self.op_H_simulate["col"]),
                )
            ).tocsr()

            ########################################
            # ⦿ 5 - evaluate 3 lowest eigenenergies and eigenvectors, |1> |2> |3>
            ##################################
            eigvals, eigvecs = eigsh(self.op_H, 3, which="SA", tol=0)
            # sort in ascending order
            sort_idx = np.argsort(eigvals)
            eigvals = np.array(eigvals)[sort_idx]
            eigvecs = np.transpose(eigvecs)[sort_idx]

            ########################################
            # ⦿ 6 - add onto the arrays
            ########################################
            self.simulation["eigvals"] = np.vstack(
                (self.simulation["eigvals"], eigvals)
            )
            self.simulation["eigvecs"] = np.vstack(
                (self.simulation["eigvecs"], eigvecs)
            )
            self.simulation["1-2"].append(eigvals[1] - eigvals[0])
            self.simulation["2-3"].append(eigvals[2] - eigvals[1])

            ########################################
            # ⦿ 7 - dipole transitions - evalute voltage in a sandwhich between the ground and excited states
            ########################################
            if simulate_dipole:
                state_0 = self.simulation["eigvecs"][0]
                state_1 = self.simulation["eigvecs"][1]
                dipole_moment_voltage = state_0.dot(self.op_V.dot(state_1))
                dipole_moment_beta = dipole_moment_voltage / (
                    self.const_Phi0 * self.simulation["1-2"][-1]
                )

                self.simulation["dipole-voltage"].append(np.abs(dipole_moment_voltage))
                self.simulation["dipole-beta"].append(np.abs(dipole_moment_beta))

            bar.update()

        if self.message_or_not:
            print("⦿ 'simulate' finished")
