import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ..parameters.martinez_reina_parameters import Martinez_Reina_Parameters
from .scheiner_model import Scheiner_Model


class Martinez_Reina_Model(Scheiner_Model):
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Martinez_Reina_Parameters()
        self.initial_guess_root = np.array([6.196390627918603e-004, 5.583931899482344e-004, 8.069635262731931e-004, self.parameters.bone_volume.vascular_pore_fraction, self.parameters.bone_volume.bone_fraction])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None
        self.ageing_queue = None
        self.denosumab_concentration_over_time = None
        self.time_for_denosumab = None
        self.bone_apparent_density = np.array([])

    def calculate_strain_energy_density(self, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t):
        stress_vector = self.calculate_macroscopic_stress_vector(t)

        compliance_matrix = self.calculate_compliance_matrix(OBa, OCa, bone_volume_fraction, t)
        strain_matrix = compliance_matrix @ stress_vector.T

        stiffness_tensor = np.linalg.inv(compliance_matrix)
        strain_energy_density = ((1 / 2) * strain_matrix.T @ stiffness_tensor @ strain_matrix)
        return strain_energy_density

    def calculate_compliance_matrix(self, OBa, OCa, bone_volume_fraction, t):
        youngs_modulus = self.calculate_youngs_modulus(OBa, OCa, bone_volume_fraction)
        # Initialize compliance tensors as 6x6 matrices of zeros
        compliance_matrix = np.zeros((6, 6))
        compliance_matrix[0, 0] = 1 / youngs_modulus
        compliance_matrix[0, 1] = -self.parameters.mechanics.poissons_ratio / youngs_modulus
        compliance_matrix[0, 2] = -self.parameters.mechanics.poissons_ratio / youngs_modulus
        compliance_matrix[1, 0] = compliance_matrix[0, 1]
        compliance_matrix[1, 1] = compliance_matrix[0, 0]
        compliance_matrix[1, 2] = -self.parameters.mechanics.poissons_ratio / youngs_modulus
        compliance_matrix[2, 0] = -self.parameters.mechanics.poissons_ratio / youngs_modulus
        compliance_matrix[2, 1] = -self.parameters.mechanics.poissons_ratio / youngs_modulus
        compliance_matrix[2, 2] = compliance_matrix[0, 0]
        compliance_matrix[3, 3] = (2 + 2 * self.parameters.mechanics.poissons_ratio) / youngs_modulus
        compliance_matrix[4, 4] = (2 + 2 * self.parameters.mechanics.poissons_ratio) / youngs_modulus
        compliance_matrix[5, 5] = (2 + 2 * self.parameters.mechanics.poissons_ratio) / youngs_modulus
        return compliance_matrix

    def calculate_youngs_modulus(self, OBa, OCa, bone_volume_fraction):
        ash_fraction = self.calculate_ash_fraction(OBa, OCa, bone_volume_fraction)
        youngs_modulus = 84.37 * ((bone_volume_fraction / 100) ** 2.58) * (ash_fraction ** 2.74)
        return youngs_modulus

    def calculate_ash_fraction(self, OBa, OCa, bone_volume_fraction):
        volume_fraction_mineral = self.calculate_average_mineral_content(OBa, OCa, bone_volume_fraction)
        ash_fraction = (self.parameters.mineralisation.density_mineral * volume_fraction_mineral /
                        (self.parameters.mineralisation.density_mineral * volume_fraction_mineral +
                         self.parameters.mineralisation.density_organic *
                         self.parameters.mineralisation.volume_fraction_organic))
        bone_material_density = (1 + (self.parameters.mineralisation.density_mineral - 1) * volume_fraction_mineral +
                                 (self.parameters.mineralisation.density_organic - 1) * self.parameters.mineralisation.volume_fraction_organic)
        bone_apparent_density = bone_material_density * self.parameters.bone_volume.bone_fraction / 100
        self.bone_apparent_density = np.append(self.bone_apparent_density, bone_apparent_density)
        return ash_fraction

    def calculate_average_mineral_content(self, OBa, OCa, bone_volume_fraction):
        # update mineralization queue
        if self.ageing_queue is None:
            self.initialise_ageing_queue(OBa, OCa, bone_volume_fraction/100)
        else:
            self.update_ageing_queue(OBa, OCa, bone_volume_fraction/100)
        sum_mineral_content = 0
        for j in range(len(self.ageing_queue) - 1):
            sum_mineral_content += self.ageing_queue[j] * self.calculate_mineralisation_law(j)
        sum_mineral_content += self.ageing_queue[-1] * self.parameters.mineralisation.maximum_mineral_content
        average_mineral_content = sum_mineral_content / (bone_volume_fraction / 100)
        return average_mineral_content

    def update_ageing_queue(self, OBa, OCa, bone_volume_fraction):
        # update mineralization queue
        resorbed_bone_fraction = OCa * self.parameters.bone_volume.resorption_rate
        formed_bone_fraction = OBa * self.parameters.bone_volume.formation_rate

        summed_queue_bone_volume_content = 0
        for i in range(len(self.ageing_queue) - 2, int(self.parameters.mineralisation.lag_time), -1):
            self.ageing_queue[i] = self.ageing_queue[i - 1] * (1 - resorbed_bone_fraction / bone_volume_fraction)
            if self.ageing_queue[i] < 1e-13:
                self.ageing_queue[i] = 0
            summed_queue_bone_volume_content += self.ageing_queue[i]
        # We assume that the tissue in the mineralisation lag time is not resorbed
        for j in range(int(self.parameters.mineralisation.lag_time), 0, -1):
            self.ageing_queue[j] = self.ageing_queue[j - 1]
            if self.ageing_queue[j] < 1e-13:
                self.ageing_queue[j] = 0
            summed_queue_bone_volume_content += self.ageing_queue[j]
        self.ageing_queue[0] = formed_bone_fraction
        summed_queue_bone_volume_content += self.ageing_queue[0]
        # The last element of the queue stores the volume needed for all the elements of VFPREV to sum (1-p)=vb
        self.ageing_queue[-1] = bone_volume_fraction + (formed_bone_fraction - resorbed_bone_fraction) - summed_queue_bone_volume_content
        if self.ageing_queue[-1] < 1e-13:
            self.ageing_queue[-1] = 0
            #print('Take care, volume of last element in the queue <0')
        #print('Updated ageing queue.')
        pass

    def initialise_ageing_queue(self, OBa, OCa, bone_volume_fraction):
        self.ageing_queue = np.zeros(self.parameters.mineralisation.length_of_queue)
        j = 0
        while j < self.parameters.mineralisation.length_of_queue:
            self.update_ageing_queue(OBa, OCa, bone_volume_fraction)
            j += 1
        print('Initialised ageing queue.')
        pass

    def calculate_mineralisation_law(self, t):
        if t <= self.parameters.mineralisation.lag_time:
            mineral_content = 0
        elif t <= self.parameters.mineralisation.primary_phase_duration + self.parameters.mineralisation.lag_time:
            mineral_content = self.parameters.mineralisation.primary_mineral_content * (
                    t - self.parameters.mineralisation.lag_time) / self.parameters.mineralisation.primary_phase_duration
        else:
            mineral_content = (self.parameters.mineralisation.maximum_mineral_content + (self.parameters.mineralisation.primary_mineral_content
                                                                                 - self.parameters.mineralisation.maximum_mineral_content) *
                       np.exp(-self.parameters.mineralisation.rate * (
                               t - self.parameters.mineralisation.primary_phase_duration - self.parameters.mineralisation.lag_time)))
        return mineral_content

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
        denosumab_effect = self.parameters.denosumab.accessibility_factor * self.parameters.binding_constant.RANKL_denosumab * self.calculate_denosumab_concentration(t)
        RANKL_eff = self.calculate_effective_carrying_capacity_RANKL(OBp, OBa, t)
        RANKL_RANK_OPG = RANKL_eff / (1 + self.parameters.binding_constant.RANKL_OPG *
                                      self.calculate_OPG_concentration(OBp, OBa, t) +
                                      self.parameters.binding_constant.RANKL_RANK * self.parameters.concentration.RANK + denosumab_effect)
        RANKL = RANKL_RANK_OPG * ((self.parameters.production_rate.intrinsic_RANKL +
                                   self.calculate_external_injection_RANKL(t)) /
                                  (self.parameters.production_rate.intrinsic_RANKL +
                                  self.parameters.degradation_rate.RANKL * RANKL_eff))
        return RANKL

    def calculate_denosumab_concentration(self, t):
        if self.denosumab_concentration_over_time is None:
                self.solve_for_denosumab_concentration()
        if t is None or t <= self.load_case.start_denosumab_treatment or t >= self.load_case.end_denosumab_treatment:
            return 0
        else:
            # calculate how many injections were given already
            number_of_injections_given = np.floor((t - self.load_case.start_denosumab_treatment) / self.load_case.treatment_period)
            # translate the current time to the interval [0, treatment_period]
            time_translation = (t - self.load_case.start_denosumab_treatment - number_of_injections_given * self.load_case.treatment_period)
            # find the closest index in the solution of the denosumab ODE
            closest_index = np.argmin(np.abs(self.time_for_denosumab - time_translation))
            # get C_den in ng/ml and calculate it to pmol/L using the molar mass (in ng/mol) of denosumab
            denosumb_concentration = self.denosumab_concentration_over_time[closest_index] / (10 ** -3) * (1 / (self.parameters.denosumab.molar_mass * 10 ** 12))
            return denosumb_concentration

    def solve_for_denosumab_concentration(self):
        initial_denosumb_concentration = 0
        t_span = [0, self.load_case.treatment_period]
        sol = solve_ivp(self.pharmacokinetics_pharmocodynamics_denosumab, t_span, [initial_denosumb_concentration], rtol=1e-10, atol=1e-10, max_step=1)
        self.denosumab_concentration_over_time = sol.y[0]  # ng/ml
        self.time_for_denosumab = sol.t
        plt.figure()
        plt.plot(sol.t, sol.y[0])
        plt.xlabel('Time [days]')
        plt.ylabel('Denosumab concentration [ng/ml]')
        # plt.show()
        print('Denosumab concentration initialized')

    def pharmacokinetics_pharmocodynamics_denosumab(self, t, concentration):
        dCdt = ((self.parameters.denosumab.absorption_rate * (self.load_case.denosumab_dose / (self.parameters.denosumab.volume_central_compartment
                                                                                              / self.parameters.denosumab.bioavailability) *
                                                             np.exp(-self.parameters.denosumab.absorption_rate * t)) -
                (concentration / (self.parameters.denosumab.michaelis_menten_constant + concentration)) *
                (self.parameters.denosumab.maximum_volume / (self.parameters.denosumab.volume_central_compartment /
                                                             self.parameters.denosumab.bioavailability))) -
                self.parameters.denosumab.elimination_rate * concentration)
        return dCdt

    def calculate_external_injection_RANKL(self, t):
        if t is None or self.load_case.start_postmenopausal_osteoporosis >= t or t >= self.load_case.end_postmenopausal_osteoporosis:
            return 0
        elif self.load_case.start_postmenopausal_osteoporosis <= t <= self.load_case.end_postmenopausal_osteoporosis:
            external_injection_RANKL = self.parameters.PMO.increase_in_RANKL * (
                    self.parameters.PMO.reduction_factor ** 2 / (self.parameters.PMO.reduction_factor ** 2 +
                                                                 ((t - self.load_case.start_postmenopausal_osteoporosis)
                                                                  / self.parameters.PMO.characteristic_time) ** 2))
            return external_injection_RANKL
        else:
            return 0


































