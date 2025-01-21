import numpy as np
import matplotlib.pyplot as plt
from parameters.martonova_parameters import Parameters
from scipy.integrate import solve_ivp


class Martonova_Model:
    """ This class defines the model by Martonova et al., 2023.
    Parameters: parameters (object) - an object of the Parameters class.
    Inputs: load_case (object) - an object of the Load_Case class, defines healthy/ disease state and treatment.
    Outputs: cellular_activity (array) - the cellular activity of the model, integrated_activity (float),
            cellular_responsiveness (float)."""
    def __init__(self, load_case):
        self.parameters = Parameters()
        # calculate active complex activity constant
        self.parameters.activity.active_complex = ((self.parameters.activity.active_receptor *
                                                    self.parameters.kinematics.receptor + self.parameters.activity.inactive_receptor)
                                                   / (self.parameters.kinematics.receptor + 1) * (
                                                           self.parameters.kinematics.complex + 1) - self.parameters.activity.inactive_complex) / self.parameters.kinematics.complex
        # define basal PTH pulse parameters based on load case
        self.initialise_basal_PTH_pulse(load_case)
        # calculate drug PTH pulse based on load case
        self.calculate_drug_PTH_pulse(load_case)
        self.initial_condition = np.array([0.9, 0, 0])
        self.number_of_periods = 30
        self.period_for_activity_constants = 20

    def initialise_basal_PTH_pulse(self, load_case):
        """ This function initialises the basal PTH pulse parameters based on the load case.
        The parameters include min/ max (= base and pulse component) of the pulse, on/off duration, and period."""
        self.parameters.basal_PTH_pulse.min = load_case.basal_PTH_pulse_min
        self.parameters.basal_PTH_pulse.max = load_case.basal_PTH_pulse_max
        self.parameters.basal_PTH_pulse.off_duration = load_case.basal_PTH_pulse_off_duration
        self.parameters.basal_PTH_pulse.on_duration = load_case.basal_PTH_pulse_on_duration
        self.parameters.basal_PTH_pulse.period = (self.parameters.basal_PTH_pulse.off_duration +
                                                  self.parameters.basal_PTH_pulse.on_duration)

    def calculate_drug_PTH_pulse(self, load_case):
        """ This function calculates the drug PTH pulse based on the load case. """
        if load_case.drug_dose is None:
            # no drug PTH injected, only basal PTH pulses are present
            pass
        else:
            # TODO: Implement drug PTH pulse calculation
            pass

    def calculate_cellular_activity(self):
        """This function calculates the cellular activity of the model based on receptor and complex concentrations.
        It corresponds to the term alpha in equation (4) in the Li & Goldbeter paper."""
        receptor_complex_concentrations = self.calculate_receptor_complex_concentrations()
        [active_receptor, active_complex, inactive_complex, inactive_receptor, time] = receptor_complex_concentrations
        # calculate cellular activity according to equation (4) in the Li & Goldbeter paper
        cellular_activity = (self.parameters.activity.active_receptor * active_receptor +
                             self.parameters.activity.active_complex * active_complex +
                             self.parameters.activity.inactive_complex * inactive_complex +
                             self.parameters.activity.inactive_receptor * inactive_receptor)
        return cellular_activity, time

    def calculate_receptor_complex_concentrations(self):
        """ This function calculates the concentrations of receptors and complexes over time."""
        sol = solve_ivp(self.receptor_ligand_model,
                        [0, self.number_of_periods * self.parameters.basal_PTH_pulse.period],
                        self.initial_condition, method="Radau", max_step=0.1)
        time = sol.t
        active_receptor = sol.y[0, :]
        active_complex = sol.y[1, :]
        inactive_complex = sol.y[2, :]
        # normalized concentrations sum up to one
        inactive_receptor = 1 - active_receptor - active_complex - inactive_complex
        return [active_receptor, active_complex, inactive_complex, inactive_receptor, time]

    def receptor_ligand_model(self, t, x):
        """ This function defines the receptor-ligand model based on an ODE system.
        It corresponds to equation (2) in the Li & Goldbeter paper."""
        active_receptor = x[0]
        active_complex = x[1]
        inactive_complex = x[2]
        # calculate PTH concentration depending on time (on or off phase)
        PTH_concentration = self.calculate_PTH_concentration(t)
        # calculate the ODE system based on equation (2) in the Li & Goldbeter paper
        d_active_receptor_dt = (
                -self.parameters.kinematics.active_complex_binding * active_receptor * PTH_concentration -
                self.parameters.kinematics.receptor_desensitized * active_receptor +
                self.parameters.kinematics.active_complex_unbinding * active_complex +
                self.parameters.kinematics.receptor_resensitized *
                (1 - active_receptor - inactive_complex - active_complex))
        d_active_complex_dt = (self.parameters.kinematics.active_complex_binding * active_receptor * PTH_concentration -
                               self.parameters.kinematics.active_complex_unbinding * active_complex -
                               self.parameters.kinematics.complex_desensitized * active_complex +
                               self.parameters.kinematics.complex_resensitized * inactive_complex)
        d_inactive_complex_dt = (self.parameters.kinematics.complex_desensitized * active_complex -
                                 self.parameters.kinematics.complex_resensitized * inactive_complex -
                                 self.parameters.kinematics.inactive_complex_unbinding * inactive_complex +
                                 self.parameters.kinematics.inactive_complex_binding * PTH_concentration *
                                 (1 - active_receptor - inactive_complex - active_complex))
        # return ODE system
        dxdt = [d_active_receptor_dt, d_active_complex_dt, d_inactive_complex_dt]
        return dxdt

    def calculate_PTH_concentration(self, t):
        """ This function calculates the PTH concentration based on the time t."""
        if self.parameters.injected_PTH_pulse.min is not None:
            # TODO: Implement injected PTH pulse calculation
            pass
        else:
            # no injected PTH, only basal PTH pulse is present
            number_of_pulses = np.floor(t / self.parameters.basal_PTH_pulse.period)
            if (number_of_pulses * self.parameters.basal_PTH_pulse.period <= t <=
                    (number_of_pulses * self.parameters.basal_PTH_pulse.period + self.parameters.basal_PTH_pulse.on_duration)):
                # gland on phase
                PTH = ((self.parameters.basal_PTH_pulse.max + self.parameters.basal_PTH_pulse.min) *
                     self.parameters.kinematics.active_binding_unbinding)
            else:  # gland off phase
                PTH = self.parameters.basal_PTH_pulse.min * self.parameters.kinematics.active_binding_unbinding
        return PTH

    def calculate_activity_constants(self, cellular_activity, sol_t):
        # calculate basal activity alpha_B
        basal_activity = self.calculate_basal_activity()
        # calculate integrated activity for step increase (alphaTstep, equation (27))
        integrated_activity_for_step_increase = self.calculate_integrated_activity_for_step_increase()
        # calculate integrated activity (alphaT, equation (26))
        chosen_activity_pulse = []
        chosen_activity_pulse_time = []
        if self.parameters.injected_PTH_pulse.min is not None:
            # TODO: Implement injected PTH pulse calculation
            pass
        else:
            for i in range(len(sol_t)):
                if (self.period_for_activity_constants * self.parameters.basal_PTH_pulse.period <= sol_t[i] <=
                        self.period_for_activity_constants * self.parameters.basal_PTH_pulse.period + self.parameters.basal_PTH_pulse.on_duration):
                    chosen_activity_pulse.append(cellular_activity[i])
                    chosen_activity_pulse_time.append(sol_t[i])
        basal_integrated_activity = np.linalg.norm(np.trapz(np.array(chosen_activity_pulse) - basal_activity, chosen_activity_pulse_time, axis=0))
        basal_cellular_responsiveness = (basal_integrated_activity / integrated_activity_for_step_increase) * (basal_integrated_activity / self.parameters.basal_PTH_pulse.period)
        return basal_activity, basal_integrated_activity, basal_cellular_responsiveness

    def calculate_basal_activity(self):
        """ This function calculates the basal activity of the cellular activity.
        It corresponds to the term alpha_0 in equation (9) in the Li & Goldbeter paper. """
        basal_activity = 1 / (1 + self.parameters.kinematics.receptor) * (
                self.parameters.activity.active_receptor * self.parameters.kinematics.receptor + self.parameters.activity.inactive_receptor)
        return basal_activity

    def calculate_integrated_activity_for_step_increase(self):
        """ This function calculates the integrated activity for a step increase in the cellular activity.
        It corresponds to the term alpha_Tstep in equation (28) in the Li & Goldbeter paper. """
        difference_in_receptor_fraction = self.calculate_difference_in_receptor_fraction(self.parameters.basal_PTH_pulse.min + self.parameters.basal_PTH_pulse.max)
        difference_in_weights = self.calculate_difference_in_weights(self.parameters.basal_PTH_pulse.min + self.parameters.basal_PTH_pulse.max)
        adaptation_time = self.calculate_adaptation_time(self.parameters.basal_PTH_pulse.min + self.parameters.basal_PTH_pulse.max)
        # eq 27
        amplitude_of_cellular_activity = difference_in_receptor_fraction * difference_in_weights
        # eq 28
        integrated_activity_for_step_increase = amplitude_of_cellular_activity * adaptation_time
        return integrated_activity_for_step_increase

    def calculate_difference_in_receptor_fraction(self, stimulus_concentration):
        """ This function calculates the difference in receptor fraction for a stimulus concentration.
        It corresponds to the term Q in equation (20) in the Li & Goldbeter paper. """
        # calculate receptor desensitisation/ resensitisation constants for basal PTH pulse minimum
        desensitised_receptors_after_adaptation_basal_min = self.calculate_desensitised_receptors_after_adaptation(
            self.parameters.basal_PTH_pulse.min)
        # calculate receptor desensitisation/ resensitisation constants for basal PTH pulse maximum
        desensitised_receptors_after_adaptation_basal_max = self.calculate_desensitised_receptors_after_adaptation(
            stimulus_concentration)
        # calculate difference between fraction of receptors for min and max PTH pulse (eq 20)
        difference_in_receptor_fraction = desensitised_receptors_after_adaptation_basal_max - desensitised_receptors_after_adaptation_basal_min
        return difference_in_receptor_fraction

    def calculate_adaptation_time(self, stimulus_concentration):
        """ This function calculates the adaptation time of the cellular activity to a stimulus.
         It corresponds to the term tau_a in equation (11) in the Li & Goldbeter paper. """
        contribution_of_receptor_desensitation_basal_max = self.calculate_contribution_of_receptor_desensitation(stimulus_concentration)
        contribution_of_receptor_resensitisation_basal_max = self.calculate_contribution_of_receptor_resensitisation(stimulus_concentration)
        adaptation_time = (1 / (contribution_of_receptor_desensitation_basal_max + contribution_of_receptor_resensitisation_basal_max))
        return adaptation_time

    def calculate_difference_in_weights(self, stimulus_concentration):
        # calculate apparent weights of active/ desensitised receptors
        weight_active_receptor = self.calculate_weight_active_receptor(stimulus_concentration)
        weight_desensitised_receptor = self.calculate_weight_desensitised_receptor(stimulus_concentration)
        difference_in_weights = weight_active_receptor - weight_desensitised_receptor
        return difference_in_weights

    def calculate_contribution_of_receptor_desensitation(self, stimulus_concentration):
        """ This function calculates the contribution of receptor desensitisation to the cellular activity.
         It corresponds to the term u in equation (12a) in the Li & Goldbeter paper. """
        contribution_of_receptor_desensitation_basal_min = (self.parameters.kinematics.receptor_desensitized +
                                                            self.parameters.kinematics.complex_desensitized *
                                                            stimulus_concentration) / (1 + stimulus_concentration)
        return contribution_of_receptor_desensitation_basal_min

    def calculate_contribution_of_receptor_resensitisation(self, stimulus_concentration):
        """ This function calculates the contribution of receptor resensitisation to the cellular activity.
         It corresponds to the term v in equation (12b) in the Li & Goldbeter paper. """
        kinetic_constant = self.parameters.kinematics.receptor / self.parameters.kinematics.complex
        contribution_of_receptor_resensitisation_basal_min = ((self.parameters.kinematics.receptor_resensitized +
                                                               self.parameters.kinematics.complex_resensitized *
                                                               stimulus_concentration * kinetic_constant) /
                                                              (1 + stimulus_concentration * kinetic_constant))
        return contribution_of_receptor_resensitisation_basal_min

    def calculate_desensitised_receptors_after_adaptation(self, stimulus_concentration):
        """ This function calculates the total fraction of desensitised receptors after adaptation to a stimulus.
         It corresponds to the term Ds in equation (8) in the Li & Goldbeter paper. """
        contribution_of_receptor_desensitation = self.calculate_contribution_of_receptor_desensitation(stimulus_concentration)
        contribution_of_receptor_resensitisation = self.calculate_contribution_of_receptor_resensitisation(stimulus_concentration)
        desensitised_receptors_after_adaptation = (contribution_of_receptor_desensitation /
                                                   (contribution_of_receptor_desensitation + contribution_of_receptor_resensitisation))
        return desensitised_receptors_after_adaptation

    def calculate_weight_active_receptor(self, stimulus_concentration):
        """ This function calculates the weight of active receptor in the cellular activity.
         It corresponds to the term a in equation (13a) in the Li & Goldbeter paper. """
        weight_active_receptor = (self.parameters.activity.active_receptor + self.parameters.activity.active_complex * stimulus_concentration
                                  / (1 + stimulus_concentration))
        return weight_active_receptor

    def calculate_weight_desensitised_receptor(self, stimulus_concentration):
        """ This function calculates the weight of desensitised receptor in the cellular activity.
         It corresponds to the term b in equation (13b) in the Li & Goldbeter paper. """
        kinetic_constant = self.parameters.kinematics.receptor / self.parameters.kinematics.complex
        weight_desensitised_receptor = (self.parameters.activity.inactive_receptor + self.parameters.activity.inactive_complex * stimulus_concentration * kinetic_constant) / (
                                               1 + stimulus_concentration * kinetic_constant)
        return weight_desensitised_receptor
