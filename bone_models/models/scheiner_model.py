import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from ..parameters.scheiner_parameters import Scheiner_Parameters
from .pivonka_model import Pivonka_Model


class Scheiner_Model(Pivonka_Model):
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Scheiner_Parameters()
        self.initial_guess_root = np.array([6.196390627918603e-004, 5.583931899482344e-004, 8.069635262731931e-004, self.parameters.bone_volume.vascular_pore_fraction, self.parameters.bone_volume.bone_fraction])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def bone_cell_population_model(self, x, t=None):
        """ Calculates the system of ordinary differential equations for the bone cell population model.
        This function is inherited by the specific models. If a function is not relevant in the specific model
        (e.g. mechanical effects), it returns 0 (additive) or 1 (multiplicative) as neutral values.

        :param x: state variables of the model
        :type x: list
        :param t: time variable
        :type t: float
        :return: rate of change of state variables
        :rtype: list
        """
        OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction = x
        dOBpdt = ((self.parameters.differentiation_rate.OBu * self.calculate_TGFb_activation_OBu(OCa,t) -
                   self.parameters.differentiation_rate.OBp * OBp *
                   self.calculate_TGFb_repression_OBp(OCa,t) + self.calculate_external_injection_OBp(t)) +
                  self.apply_mechanical_effects(OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t) + self.apply_medication_effects_OBp())
        dOBadt = (self.parameters.differentiation_rate.OBp * OBp * self.calculate_TGFb_repression_OBp(OCa,t) -
                  self.parameters.apoptosis_rate.OBa * OBa * self.apply_medication_effects_OBa()
                  + self.calculate_external_injection_OBa(t))
        dOCadt = (self.parameters.differentiation_rate.OCp * self.calculate_RANKL_activation_OCp(OBp,OBa,t) -
                  self.parameters.apoptosis_rate.OCa * OCa * self.calculate_TGFb_activation_OCa(OCa,t) +
                  self.calculate_external_injection_OCa(t))
        dvascular_pore_fractiondt = self.parameters.bone_volume.resorption_rate * OCa - self.parameters.bone_volume.formation_rate * OBa
        dbone_volume_fractiondt = self.parameters.bone_volume.formation_rate * OBa - self.parameters.bone_volume.resorption_rate * OCa
        dxdt = [dOBpdt, dOBadt, dOCadt, dvascular_pore_fractiondt, dbone_volume_fractiondt]
        return dxdt

    def solve_bone_cell_population_model(self, tspan):
        """ Solve the bone cell population model using the ODE system over a given time interval.
        The initial conditions are set to the steady-state values.

        :param tspan: time span for the ODE solver
        :type tspan: numpy.ndarray with start and end time
        :return: solution of the ODE system
        :rtype: scipy.integrate._ivp.ivp.OdeResult
        """
        x0 = self.calculate_steady_state()
        print('Solving bone cell population model ...', end='')
        solution = solve_ivp(lambda t, x: self.bone_cell_population_model(x, t), tspan, [x0[0], x0[1], x0[2], self.parameters.bone_volume.vascular_pore_fraction, self.parameters.bone_volume.bone_fraction], rtol=1e-8, atol=1e-8)
        print('done')
        return solution

    def calculate_TGFb_concentration(self, OCa, t):
        """ Calculates the TGF-beta concentration based on the osteoclastic resorption, external injection and
        degradation rate.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: TGF-beta concentration
        :rtype: float"""
        TGFb = (self.parameters.bone_volume.stored_TGFb_content * OCa +
                self.calculate_external_injection_TGFb(t)) / self.parameters.degradation_rate.TGFb
        return TGFb

    def calculate_RANKL_concentration(self, OBp, OBa, t):
        """ Calculates the RANKL concentration based on the effective carrying capacity, RANKL-RANK-OPG binding,
        degradation rate, intrinsic RANKL production and external injection of RANKL.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: RANKL concentration
        :rtype: float"""
        RANKL_eff = self.calculate_effective_carrying_capacity_RANKL(OBp, OBa, t)
        RANKL_RANK_OPG = RANKL_eff / (1 + self.parameters.binding_constant.RANKL_OPG *
                                      self.calculate_OPG_concentration(OBp, OBa, t) +
                                      self.parameters.binding_constant.RANKL_RANK * self.parameters.concentration.RANK)
        RANKL = RANKL_RANK_OPG * ((self.parameters.production_rate.intrinsic_RANKL +
                                   self.calculate_external_injection_RANKL(t) + self.parameters.mechanics.RANKL_production) /
                                  (self.parameters.production_rate.intrinsic_RANKL +
                                  self.parameters.degradation_rate.RANKL * RANKL_eff))
        return RANKL

    def apply_mechanical_effects(self, OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t):
        if t is None:
            # no mechanical effects in steady-state
            return 0
        else:
            strain_effect_on_OBp = self.calculate_strain_effect_on_OBp(vascular_pore_fraction, bone_volume_fraction, t)
            if self.parameters.mechanics.update_OBp_proliferation_rate:
                self.parameters.proliferation_rate.OBp = ((self.parameters.differentiation_rate.OBu *
                                                           self.parameters.mechanics.fraction_of_OBu_differentiation_rate *
                                                           self.calculate_TGFb_activation_OBu(OCa, t)) /
                                                          self.steady_state.OBp / strain_effect_on_OBp)
            return (self.parameters.differentiation_rate.OBu*(-self.parameters.mechanics.fraction_of_OBu_differentiation_rate)*self.calculate_TGFb_activation_OBu(OCa, t) +
                    self.parameters.proliferation_rate.OBp * OBp * strain_effect_on_OBp)

    def calculate_strain_effect_on_OBp(self, vascular_pore_fraction, bone_volume_fraction, t):
        if t <= self.load_case.start_time:
            if self.parameters.mechanics.strain_energy_density_steady_state is None:
                self.parameters.mechanics.strain_energy_density_steady_state = self.calculate_strain_energy_density(vascular_pore_fraction, bone_volume_fraction,t)
            return self.parameters.mechanics.strain_effect_on_OBp_steady_state
        else:
            self.parameters.mechanics.update_OBp_proliferation_rate = False
            strain_energy_density = self.calculate_strain_energy_density(vascular_pore_fraction, bone_volume_fraction, t)
            if strain_energy_density > self.parameters.mechanics.strain_energy_density_steady_state:
                strain_effect_on_OBp = self.parameters.mechanics.strain_effect_on_OBp_steady_state * (
                            1 + 1.25 * (strain_energy_density
                                        - self.parameters.mechanics.strain_energy_density_steady_state)
                            / self.parameters.mechanics.strain_energy_density_steady_state)
                self.parameters.mechanics.RANKL_production = 0
                if strain_effect_on_OBp > 1:
                    return 1
                else:
                    return strain_effect_on_OBp
            elif strain_energy_density < self.parameters.mechanics.strain_energy_density_steady_state:
                self.parameters.mechanics.RANKL_production = (-10 ** 5 * (strain_energy_density -
                                                                          self.parameters.mechanics.strain_energy_density_steady_state)
                                                              / self.parameters.mechanics.strain_energy_density_steady_state)
                return self.parameters.mechanics.strain_effect_on_OBp_steady_state
            else:
                self.parameters.mechanics.RANKL_production = 0
                return self.parameters.mechanics.strain_effect_on_OBp_steady_state

    def calculate_strain_energy_density(self, vascular_pore_fraction, bone_volume_fraction, t):
        microscopic_strain_tensor = self.calculate_microscopic_strain_tensor(vascular_pore_fraction, bone_volume_fraction, t)
        return (1 / 2) * (microscopic_strain_tensor.T @ self.parameters.mechanics.stiffness_tensor_bone_matrix @ microscopic_strain_tensor)

    def calculate_microscopic_strain_tensor(self, vascular_pore_fraction, bone_volume_fraction, t):
        [strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores] = self.calculate_strain_concentration_tensors(bone_volume_fraction)
        return strain_concentration_tensor_bone_matrix @ self.calculate_macroscopic_strain_tensor(strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores, vascular_pore_fraction, bone_volume_fraction, t)

    def calculate_macroscopic_strain_tensor(self, strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores, vascular_pore_fraction, bone_volume_fraction, t):
        macroscopic_stiffness_tensor = self.calculate_macroscopic_stiffness_tensor(strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores, vascular_pore_fraction, bone_volume_fraction)
        return np.linalg.inv(macroscopic_stiffness_tensor) @ self.calculate_macroscopic_stress_vector(t).T

    def calculate_macroscopic_stiffness_tensor(self, strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores, vascular_pore_fraction, bone_volume_fraction):
        macroscopic_stiffness_tensor = ((vascular_pore_fraction / 100) * self.parameters.mechanics.stiffness_tensor_vascular_pores
                                        @ strain_concentration_tensor_vascular_pores +
                                        (bone_volume_fraction / 100) * self.parameters.mechanics.stiffness_tensor_bone_matrix
                                        @ strain_concentration_tensor_bone_matrix)
        return macroscopic_stiffness_tensor

    def calculate_strain_concentration_tensors(self, bone_volume_fraction):
        hill_tensor_cylindrical_inclusion = self.calculate_hill_tensor_cylindrical_inclusion()
        strain_concentration = np.linalg.inv((bone_volume_fraction / 100 *
                                np.linalg.inv(self.parameters.mechanics.unit_tensor_as_matrix +
                                              hill_tensor_cylindrical_inclusion @
                                              (self.parameters.mechanics.stiffness_tensor_vascular_pores -
                                               self.parameters.mechanics.stiffness_tensor_bone_matrix))) +
                                                bone_volume_fraction/100 * self.parameters.mechanics.unit_tensor_as_matrix)
        strain_concentration_tensor_bone_matrix = self.parameters.mechanics.unit_tensor_as_matrix @ strain_concentration
        strain_concentration_tensor_vascular_pores = (np.linalg.inv(self.parameters.mechanics.unit_tensor_as_matrix +
                                                                   hill_tensor_cylindrical_inclusion @
                                                                   (self.parameters.mechanics.stiffness_tensor_vascular_pores -
                                                                    self.parameters.mechanics.stiffness_tensor_bone_matrix))
                                                      @ strain_concentration)
        return strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores

    def calculate_macroscopic_stress_vector(self, t):
        if t <= self.load_case.start_time or t >= self.load_case.end_time:
            stress_matrix = self.parameters.mechanics.stress_tensor_normal_loading
        else:
            stress_matrix = self.load_case.stress_tensor
        stress_vector = np.array([stress_matrix[0, 0], stress_matrix[1, 1], stress_matrix[2, 2], 2 * stress_matrix[0, 1],
                           2 * stress_matrix[1, 2], 2 * stress_matrix[2, 0]])
        return stress_vector
    
    def calculate_hill_tensor_cylindrical_inclusion(self):
        if self.parameters.mechanics.hill_tensor_cylindrical_inclusion is not None:
            return self.parameters.mechanics.hill_tensor_cylindrical_inclusion
        else:
            stiffness_tensor_bone_matrix = self.stiffness_matrix_to_tensor()
            unit_vector_1 = np.array([1, 0, 0])
            unit_vector_2 = np.array([0, 1, 0])

            tensor_cylindrical_inclusion = np.zeros((3, 3, 3, 3))
            help_tensor_for_G = np.zeros((3, 3, 3))
            xi_G_xi = np.zeros((3, 3, 3, 3))
            Gamma_LOCAL_inc = np.zeros((3, 3, 3, 3))

            for phi_1 in np.arange(0, 2 * np.pi, self.parameters.mechanics.step_size_for_Hill_tensor_integration):
                K_help = np.zeros((3, 3, 3))
                K = np.zeros((3, 3))
                # unit length vector with theta = pi/2
                xi = np.cos(phi_1) * unit_vector_1 + np.sin(phi_1) * unit_vector_2
                # Compute K_help tensor
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            K_help[i, j, k] = xi[0] * stiffness_tensor_bone_matrix[0, i, j, k] + xi[1] * stiffness_tensor_bone_matrix[1, i, j, k] + xi[2] * stiffness_tensor_bone_matrix[2, i, j, k]
                # Compute K matrix (acoustic tensor)
                for i in range(3):
                    for j in range(3):
                        K[i, j] = K_help[i, j, 0] * xi[0] + K_help[i, j, 1] * xi[1] + K_help[i, j, 2] * xi[2]
                # Compute K_inv matrix (= G)
                K_inv = np.linalg.inv(K)
                # Compute help tensor - first nested loop for K_inv * xi
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            help_tensor_for_G[i, j, k] = K_inv[i, j] * xi[k]
                # Second nested loop for xi * K_inv
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                xi_G_xi[i, j, k, l] = xi[i] * help_tensor_for_G[j, k, l]
                # Third nested loop for Gamma in each increment
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                Gamma_LOCAL_inc[i, j, k, l] = (xi_G_xi[i, j, k, l] + xi_G_xi[i, j, l, k] + xi_G_xi[
                                    j, i, k, l] + xi_G_xi[j, i, l, k]) / 4.0
                # Update Pcyl_LOCAL with the increment for numerical integration
                tensor_cylindrical_inclusion += Gamma_LOCAL_inc * self.parameters.mechanics.step_size_for_Hill_tensor_integration / (2 * np.pi)
                # calculate matrix equivalent of the Hill tensor
            Pcyl = np.array([
                [
                    1 * tensor_cylindrical_inclusion[0, 0, 0, 0], 1 * tensor_cylindrical_inclusion[0, 0, 1, 1], 1 * tensor_cylindrical_inclusion[0, 0, 2, 2],
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 0, 1, 2], np.sqrt(2) * tensor_cylindrical_inclusion[0, 0, 2, 0],
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 0, 0, 1]
                ],
                [
                    1 * tensor_cylindrical_inclusion[1, 1, 0, 0], 1 * tensor_cylindrical_inclusion[1, 1, 1, 1], 1 * tensor_cylindrical_inclusion[1, 1, 2, 2],
                    np.sqrt(2) * tensor_cylindrical_inclusion[1, 1, 1, 2], np.sqrt(2) * tensor_cylindrical_inclusion[1, 1, 2, 0],
                    np.sqrt(2) * tensor_cylindrical_inclusion[1, 1, 0, 1]
                ],
                [
                    1 * tensor_cylindrical_inclusion[2, 2, 0, 0], 1 * tensor_cylindrical_inclusion[2, 2, 1, 1], 1 * tensor_cylindrical_inclusion[2, 2, 2, 2],
                    np.sqrt(2) * tensor_cylindrical_inclusion[2, 2, 1, 2], np.sqrt(2) * tensor_cylindrical_inclusion[2, 2, 2, 0],
                    np.sqrt(2) * tensor_cylindrical_inclusion[2, 2, 0, 1]
                ],
                [
                    np.sqrt(2) * tensor_cylindrical_inclusion[1, 2, 0, 0], np.sqrt(2) * tensor_cylindrical_inclusion[1, 2, 1, 1],
                    np.sqrt(2) * tensor_cylindrical_inclusion[1, 2, 2, 2],
                    2 * tensor_cylindrical_inclusion[1, 2, 1, 2], 2 * tensor_cylindrical_inclusion[1, 2, 0, 2], 2 * tensor_cylindrical_inclusion[1, 2, 0, 1]
                ],
                [
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 2, 0, 0], np.sqrt(2) * tensor_cylindrical_inclusion[0, 2, 1, 1],
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 2, 2, 2],
                    2 * tensor_cylindrical_inclusion[0, 2, 1, 2], 2 * tensor_cylindrical_inclusion[0, 2, 0, 2], 2 * tensor_cylindrical_inclusion[0, 2, 0, 1]
                ],
                [
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 1, 0, 0], np.sqrt(2) * tensor_cylindrical_inclusion[0, 1, 1, 1],
                    np.sqrt(2) * tensor_cylindrical_inclusion[0, 1, 2, 2],
                    2 * tensor_cylindrical_inclusion[0, 1, 1, 2], 2 * tensor_cylindrical_inclusion[0, 1, 0, 2], 2 * tensor_cylindrical_inclusion[0, 1, 0, 1]
                ]
            ])
            self.parameters.mechanics.hill_tensor_cylindrical_inclusion = Pcyl
            return Pcyl
        
    def stiffness_matrix_to_tensor(self):
        # Reshape cbm into a 3x3x3x3 tensor
        stiffness_matrix = self.parameters.mechanics.stiffness_tensor_bone_matrix
        stiffness_tensor = np.zeros((3, 3, 3, 3))

        # Group 1: (1,1,*,*)
        stiffness_tensor[0,0,0,0] = stiffness_matrix[0,0]
        stiffness_tensor[0,0,1,1] = stiffness_matrix[0,1]
        stiffness_tensor[0,0,2,2] = stiffness_matrix[0,2]
        stiffness_tensor[0,0,1,2] = stiffness_matrix[0,3]/np.sqrt(2)
        stiffness_tensor[0,0,2,1] = stiffness_matrix[0,3]/np.sqrt(2)
        stiffness_tensor[0,0,0,2] = stiffness_matrix[0,4]/np.sqrt(2)
        stiffness_tensor[0,0,2,0] = stiffness_matrix[0,4]/np.sqrt(2)
        stiffness_tensor[0,0,0,1] = stiffness_matrix[0,5]/np.sqrt(2)
        stiffness_tensor[0,0,1,0] = stiffness_matrix[0,5]/np.sqrt(2)
        
        # Group 2: (2,2,*,*)
        stiffness_tensor[1,1,0,0] = stiffness_matrix[1,0]
        stiffness_tensor[1,1,1,1] = stiffness_matrix[1,1]
        stiffness_tensor[1,1,2,2] = stiffness_matrix[1,2]
        stiffness_tensor[1,1,1,2] = stiffness_matrix[1,3]/np.sqrt(2)
        stiffness_tensor[1,1,2,1] = stiffness_matrix[1,3]/np.sqrt(2)
        stiffness_tensor[1,1,0,2] = stiffness_matrix[1,4]/np.sqrt(2)
        stiffness_tensor[1,1,2,0] = stiffness_matrix[1,4]/np.sqrt(2)
        stiffness_tensor[1,1,0,1] = stiffness_matrix[1,5]/np.sqrt(2)
        stiffness_tensor[1,1,1,0] = stiffness_matrix[1,5]/np.sqrt(2)
        
        # Group 3: (3,3,*,*)
        stiffness_tensor[2,2,0,0] = stiffness_matrix[2,0]
        stiffness_tensor[2,2,1,1] = stiffness_matrix[2,1]
        stiffness_tensor[2,2,2,2] = stiffness_matrix[2,2]
        stiffness_tensor[2,2,1,2] = stiffness_matrix[2,3]/np.sqrt(2)
        stiffness_tensor[2,2,2,1] = stiffness_matrix[2,3]/np.sqrt(2)
        stiffness_tensor[2,2,0,2] = stiffness_matrix[2,4]/np.sqrt(2)
        stiffness_tensor[2,2,2,0] = stiffness_matrix[2,4]/np.sqrt(2)
        stiffness_tensor[2,2,0,1] = stiffness_matrix[2,5]/np.sqrt(2)
        stiffness_tensor[2,2,1,0] = stiffness_matrix[2,5]/np.sqrt(2)
        
        # Group 4: (2,3,*,*) and (3,2,*,*)
        stiffness_tensor[1,2,0,0] = stiffness_matrix[3,0]/np.sqrt(2)
        stiffness_tensor[2,1,0,0] = stiffness_matrix[3,0]/np.sqrt(2)
        stiffness_tensor[1,2,1,1] = stiffness_matrix[3,1]/np.sqrt(2)
        stiffness_tensor[2,1,1,1] = stiffness_matrix[3,1]/np.sqrt(2)
        stiffness_tensor[1,2,2,2] = stiffness_matrix[3,2]/np.sqrt(2)
        stiffness_tensor[2,1,2,2] = stiffness_matrix[3,2]/np.sqrt(2)
        stiffness_tensor[1,2,1,2] = stiffness_matrix[3,3]/2
        stiffness_tensor[1,2,2,1] = stiffness_matrix[3,3]/2
        stiffness_tensor[2,1,2,1] = stiffness_matrix[3,3]/2
        stiffness_tensor[2,1,1,2] = stiffness_matrix[3,3]/2
        stiffness_tensor[1,2,0,2] = stiffness_matrix[3,4]/2
        stiffness_tensor[1,2,2,0] = stiffness_matrix[3,4]/2
        stiffness_tensor[2,1,2,0] = stiffness_matrix[3,4]/2
        stiffness_tensor[2,1,0,2] = stiffness_matrix[3,4]/2
        stiffness_tensor[1,2,0,1] = stiffness_matrix[3,5]/2
        stiffness_tensor[1,2,1,0] = stiffness_matrix[3,5]/2
        stiffness_tensor[2,1,1,0] = stiffness_matrix[3,5]/2
        stiffness_tensor[2,1,0,1] = stiffness_matrix[3,5]/2
        
        # Group 5: (1,3,*,*) and (3,1,*,*)
        stiffness_tensor[0,2,0,0] = stiffness_matrix[4,0]/np.sqrt(2)
        stiffness_tensor[2,0,0,0] = stiffness_matrix[4,0]/np.sqrt(2)
        stiffness_tensor[0,2,1,1] = stiffness_matrix[4,1]/np.sqrt(2)
        stiffness_tensor[2,0,1,1] = stiffness_matrix[4,1]/np.sqrt(2)
        stiffness_tensor[0,2,2,2] = stiffness_matrix[4,2]/np.sqrt(2)
        stiffness_tensor[2,0,2,2] = stiffness_matrix[4,2]/np.sqrt(2)
        stiffness_tensor[0,2,1,2] = stiffness_matrix[4,3]/2
        stiffness_tensor[0,2,2,1] = stiffness_matrix[4,3]/2
        stiffness_tensor[2,0,2,1] = stiffness_matrix[4,3]/2
        stiffness_tensor[2,0,1,2] = stiffness_matrix[4,3]/2
        stiffness_tensor[0,2,0,2] = stiffness_matrix[4,4]/2
        stiffness_tensor[0,2,2,0] = stiffness_matrix[4,4]/2
        stiffness_tensor[2,0,2,0] = stiffness_matrix[4,4]/2
        stiffness_tensor[2,0,0,2] = stiffness_matrix[4,4]/2
        stiffness_tensor[0,2,0,1] = stiffness_matrix[4,5]/2
        stiffness_tensor[0,2,1,0] = stiffness_matrix[4,5]/2
        stiffness_tensor[2,0,1,0] = stiffness_matrix[4,5]/2
        stiffness_tensor[2,0,0,1] = stiffness_matrix[4,5]/2
        
        # Group 6: (1,2,*,*) and (2,1,*,*)
        stiffness_tensor[0,1,0,0] = stiffness_matrix[5,0]/np.sqrt(2)
        stiffness_tensor[1,0,0,0] = stiffness_matrix[5,0]/np.sqrt(2)
        stiffness_tensor[0,1,1,1] = stiffness_matrix[5,1]/np.sqrt(2)
        stiffness_tensor[1,0,1,1] = stiffness_matrix[5,1]/np.sqrt(2)
        stiffness_tensor[0,1,2,2] = stiffness_matrix[5,2]/np.sqrt(2)
        stiffness_tensor[1,0,2,2] = stiffness_matrix[5,2]/np.sqrt(2)
        stiffness_tensor[0,1,1,2] = stiffness_matrix[5,3]/2
        stiffness_tensor[0,1,2,1] = stiffness_matrix[5,3]/2
        stiffness_tensor[1,0,2,1] = stiffness_matrix[5,3]/2
        stiffness_tensor[1,0,1,2] = stiffness_matrix[5,3]/2
        stiffness_tensor[0,1,0,2] = stiffness_matrix[5,4]/2
        stiffness_tensor[0,1,2,0] = stiffness_matrix[5,4]/2
        stiffness_tensor[1,0,2,0] = stiffness_matrix[5,4]/2
        stiffness_tensor[1,0,0,2] = stiffness_matrix[5,4]/2
        stiffness_tensor[0,1,0,1] = stiffness_matrix[5,5]/2
        stiffness_tensor[0,1,1,0] = stiffness_matrix[5,5]/2
        stiffness_tensor[1,0,1,0] = stiffness_matrix[5,5]/2
        stiffness_tensor[1,0,0,1] = stiffness_matrix[5,5]/2
        return stiffness_tensor
        