import numpy as np
import matplotlib.pyplot as plt
from bone_models.bone_cell_population_models.parameters.martonova_parameters import Martonova_Parameters
from scipy.integrate import solve_ivp


class Martonova_Model:
    """ This class defines the two-state receptor model by Martonova et al. (2023) for pulsatile endogenous PTH in
    healthy and disease state including drug administration.

    .. note::
       **Source Publication**:
       Martonova D., Lavaill M., Forwood M.R., Robling A., Cooper D.M.L., Leyendecker S., Pivonka P. (2023).
       *Effects of PTH glandular and external dosing patterns on bone cell activity using a two-state receptor model
       —Implications for bone disease progression and treatment*
       PLOS ONE, 18(3), e0283544.
       :doi:`10.1371/journal.pone.0283544`

       It is based on the model by Li and Goldbeter (1986) for receptor-ligand binding, which will be referred to in
       some of the implemented model functions.

       Li Y., Goldbeter A. (1989).
       *Frequency specificity in intercellular communication. Influence of patterns of periodic signaling on target cell responsiveness.*
       Biophysical journal, 55(1), 125–145.
       :doi:`10.1016/S0006-3495(89)82785-7`

    :param load_case: Load case for the model
    :type load_case: Martonova_Load_Case
    :param parameters: Model parameters
    :type parameters: Martonova_Parameters
    :param initial_condition: Initial condition for the ODE system for receptor-ligand binding.
    :type initial_condition: numpy.ndarray
    :param number_of_periods: Number of periods for basal PTH pulse.
    :type number_of_periods: int
    :param period_for_activity_constants: Period for calculating activity constants after cellular adaptation to basal PTH.
    :type period_for_activity_constants: int """
    def __init__(self, load_case):
        """ Constructor method for the Martonova_Model class.
        It calculates the active complex activity constant and the drug PTH pulse based on the load case.

        :param load_case: Load case for the model, see :class:`Martonova_Load_Case` for details.
        :type load_case: object """
        self.parameters = Martonova_Parameters()
        # calculate active complex activity constant
        self.parameters.activity.active_complex = ((self.parameters.activity.active_receptor *
                                                    self.parameters.kinematics.receptor + self.parameters.activity.inactive_receptor)
                                                   / (self.parameters.kinematics.receptor + 1) * (
                                                           self.parameters.kinematics.complex + 1) - self.parameters.activity.inactive_complex) / self.parameters.kinematics.complex
        # define basal PTH pulse parameters based on load case
        self.load_case = load_case
        # calculate drug PTH pulse based on load case
        self.calculate_drug_PTH_pulse(load_case)
        self.initial_condition = np.array([0.9, 0, 0])
        self.number_of_periods = 200
        self.period_for_activity_constants = 100

    def calculate_drug_PTH_pulse(self, load_case):
        """ This function calculates the drug PTH pulse based on the load case. If there is a drug PTH injection, the
        ODE system (pharmacokinetics-pharmacodynamics model) is solved to determine the drug PTH pulse parameters for a
        square-wave pulse approximation (on/off-phase, pulse height).

        :param load_case: Load case for the model
        :type load_case: Martonova_Load_Case
        :return: None
        :rtype: None """
        if load_case.drug_dose is None:
            # no drug PTH injected, only basal PTH pulses are present
            pass
        else:
            # [pmol] = mg * 10^6 pg/mg / 4117.8 pmol/pg
            init_dose = load_case.drug_dose * (10 ** 6) / 4117.8  # convert mg to pM
            init_PTH = self.load_case.basal_PTH_pulse.min * 1000  # convert nM to pM
            init_area = 0
            x0 = np.array([init_dose, init_PTH, init_area])

            # solve ODE for one compartment model, stop when PTH is at baseline again (event function)
            def event(t, x):
                return x[1] - init_PTH

            event.terminal = True
            event.direction = -1
            sol = solve_ivp(self.PK_PD_model, [0, load_case.injection_frequency * 60], x0, rtol=1e-7, events=event)
            dose = sol.y[0, :]
            drug_concentration = sol.y[1, :]
            area = sol.y[2, :]
            maximum_concentration = np.max(drug_concentration)
            area_without_basal_PTH = area[-1] - init_PTH * sol.t[-1]
            # calculate drug PTH pulse parameters for square wave pulse approximation
            self.load_case.injected_PTH_pulse.max = (maximum_concentration - init_PTH) / 1000  # convert pM to nM
            self.load_case.injected_PTH_pulse.on_duration = (area_without_basal_PTH / (
                        maximum_concentration - init_PTH)) * 60  # convert h to min
            self.load_case.injected_PTH_pulse.off_duration = load_case.injection_frequency * 60 - self.load_case.injected_PTH_pulse.on_duration
            self.load_case.injected_PTH_pulse.period = self.load_case.injected_PTH_pulse.on_duration + self.load_case.injected_PTH_pulse.off_duration
            pass

    def PK_PD_model(self, t, x):
        """ This function defines the pharmacokinetics-pharmacodynamics model based on an ODE system for drug PTH pulse.
        The ODE system describes the change of dose, drug concentration, and area under the curve over time depending on
        absorption, elimination, distribution rates and bioavailability.

        :param t: time variable
        :type t: float
        :param x: concentrations of dose, drug concentration, and area under the curve
        :type x: list
        :return: dxdt (change of concentrations of dose, drug concentration, and area under the curve with t)
        :rtype: list """
        dose = x[0]
        drug_concentration = x[1]
        d_dose_dt = -self.parameters.pharmacokinetics.absorption_rate * dose * self.parameters.pharmacokinetics.bioavailability
        d_drug_concentration_dt = ((self.parameters.pharmacokinetics.bioavailability / self.parameters.pharmacokinetics.volume_of_distribution) *
                                   self.parameters.pharmacokinetics.absorption_rate * dose - self.parameters.pharmacokinetics.elimination_rate * drug_concentration)
        d_area_dt = drug_concentration
        d_x_dt = [d_dose_dt, d_drug_concentration_dt, d_area_dt]
        return d_x_dt

    def solve_for_activity(self):
        """ This function calculates the time dependent cellular activity and the activity constants basal activity,
        cellular responsiveness, and integrated activity.

        :return: cellular activity, time, basal activity, integrated activity, cellular responsiveness
        :rtype: list """
        cellular_activity, time = self.calculate_cellular_activity()
        basal_activity, integrated_activity, cellular_responsiveness = self.calculate_activity_constants(cellular_activity, time)
        return cellular_activity, time, basal_activity, integrated_activity, cellular_responsiveness

    def calculate_cellular_activity(self):
        """This function calculates the cellular activity of the model based on receptor and complex concentrations.
        It corresponds to the term alpha(t) in equation (4) in the Li & Goldbeter paper.

        :return: cellular activity and time
        :rtype: list"""
        receptor_complex_concentrations = self.calculate_receptor_complex_concentrations()
        [active_receptor, active_complex, inactive_complex, inactive_receptor, time] = receptor_complex_concentrations
        # calculate cellular activity according to equation (4) in the Li & Goldbeter paper
        cellular_activity = (self.parameters.activity.active_receptor * active_receptor +
                             self.parameters.activity.active_complex * active_complex +
                             self.parameters.activity.inactive_complex * inactive_complex +
                             self.parameters.activity.inactive_receptor * inactive_receptor)
        return cellular_activity, time

    def calculate_receptor_complex_concentrations(self):
        """ This function calculates the concentrations of receptors and complexes over time by solving the ODE system
        for the receptor-ligand model.

        :return: concentrations of active receptor, active complex, inactive complex, inactive receptor, and time
        :rtype: list"""
        sol = solve_ivp(self.receptor_ligand_model,
                        [0.0, self.number_of_periods * self.load_case.basal_PTH_pulse.period],
                        self.initial_condition, method="Radau", max_step=0.1)
        time = sol.t
        active_receptor = sol.y[0, :]
        active_complex = sol.y[1, :]
        inactive_complex = sol.y[2, :]
        # normalized concentrations sum up to one
        inactive_receptor = 1 - active_receptor - active_complex - inactive_complex
        return [active_receptor, active_complex, inactive_complex, inactive_receptor, time]

    def receptor_ligand_model(self, t, x):
        """ This function defines the receptor-ligand model based on an ODE system depending on the current PTH concentration.
        The system describes the change of concentrations of active receptor, active complex, and inactive complex over
        time. The concentration of inactive receptor is calculated based on the sum of the other concentrations.
        The receptor-ligand model corresponds to equation (2) in the Li & Goldbeter paper.

        :param t: time variable
        :type t: float
        :param x: concentrations of active receptor, active complex, and inactive complex
        :type x: list
        :return: ODE system (change of concentrations of active receptor, active complex, and inactive complex)
        :rtype: list"""
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
        """ This function calculates the PTH concentration depending on time, basal PTH pulse and injected PTH pulse
        (if present).

        :param t: time variable
        :type t: float

        :return: PTH concentration
        :rtype: float"""
        # determine the number of pulses for basal PTH
        glandular_pulse = np.floor(t / self.load_case.basal_PTH_pulse.period)
        if self.load_case.injected_PTH_pulse.max is not None:
            # injected PTH is present
            PTH = self.calculate_PTH_concentration_from_basal_and_injected(t, glandular_pulse)
        else:
            # no injected PTH, only basal PTH pulse is present
            PTH = self.calculate_PTH_concentration_from_basal(t, glandular_pulse)
        return PTH

    def calculate_PTH_concentration_from_basal(self, t, glandular_pulse):
        """ This function calculates the PTH concentration based on the basal PTH pulse if no injected PTH is present.
        The PTH concentration depends on non-pulsatile (off-phase) or pulsatile (on-phase) share of the basal PTH pulse.

        :param t: time variable
        :type t: float
        :param glandular_pulse: number of pulses for basal PTH
        :type glandular_pulse: int

        :return: PTH concentration
        :rtype: float"""
        if (glandular_pulse * self.load_case.basal_PTH_pulse.period <= t <=
                (glandular_pulse * self.load_case.basal_PTH_pulse.period + self.load_case.basal_PTH_pulse.on_duration)):
            PTH = ((self.load_case.basal_PTH_pulse.max + self.load_case.basal_PTH_pulse.min) *
                   self.parameters.kinematics.active_binding_unbinding)
        else:
            PTH = self.load_case.basal_PTH_pulse.min * self.parameters.kinematics.active_binding_unbinding
        return PTH

    def calculate_PTH_concentration_from_basal_and_injected(self, t, glandular_pulse):
        """ This function calculates the PTH concentration based on the basal and injected PTH pulse if injected PTH is
        present. The PTH concentration depends on non-pulsatile (off-phase) or pulsatile (on-phase) share of the basal
        PTH pulse and injected PTH pulse.

        :param t: time variable
        :type t: float
        :param glandular_pulse: number of pulses for basal PTH
        :type glandular_pulse: int

        :return: PTH concentration
        :rtype: float"""
        injected_pulse = np.floor(t / self.load_case.injected_PTH_pulse.period)
        # determine if injected PTH and/or basal PTH pulse is active (on-phase or off-phase)
        if (injected_pulse * self.load_case.injected_PTH_pulse.period <= t <= injected_pulse *
                self.load_case.injected_PTH_pulse.period + self.load_case.injected_PTH_pulse.on_duration):
            if (glandular_pulse * self.load_case.basal_PTH_pulse.period <= t <= glandular_pulse *
                    self.load_case.basal_PTH_pulse.period + self.load_case.basal_PTH_pulse.on_duration):
                PTH = (self.load_case.basal_PTH_pulse.max + self.load_case.basal_PTH_pulse.min +
                       self.load_case.injected_PTH_pulse.max) * self.parameters.kinematics.active_binding_unbinding
            else:
                PTH = ((self.load_case.basal_PTH_pulse.min + self.load_case.injected_PTH_pulse.max) *
                       self.parameters.kinematics.active_binding_unbinding)
        else:
            if (glandular_pulse * self.load_case.basal_PTH_pulse.period <= t <= glandular_pulse *
                    self.load_case.basal_PTH_pulse.period + self.load_case.basal_PTH_pulse.on_duration):
                PTH = ((self.load_case.basal_PTH_pulse.min + self.load_case.basal_PTH_pulse.max) *
                       self.parameters.kinematics.active_binding_unbinding)
            else:
                PTH = self.load_case.basal_PTH_pulse.min * self.parameters.kinematics.active_binding_unbinding
        return PTH

    def calculate_activity_constants(self, cellular_activity, time):
        """ This function calculates the activity constants of the cellular activity: basal activity (alpha_0), integrated activity (alpha_T), and cellular responsiveness (alpha_R).
        In the case of injected PTH, the activity constants are calculated for basal and injected PTH pulses. Otherwise, only the basal PTH pulse is considered.

        :param cellular_activity: cellular activity (alpha(t) in original publication)
        :type cellular_activity: list
        :param time: time variable
        :type time: list

        :return: basal activity, integrated activity, cellular responsiveness
        :rtype: list"""
        basal_activity = self.calculate_basal_activity()
        integrated_activity_for_step_increase = self.calculate_integrated_activity_for_step_increase(self.load_case.basal_PTH_pulse.min + self.load_case.basal_PTH_pulse.max)
        if self.load_case.injected_PTH_pulse.max is not None:
            integrated_activity, cellular_responsiveness = self.calculate_activity_constants_for_basal_and_injection(basal_activity, integrated_activity_for_step_increase, cellular_activity, time)
        else:
            integrated_activity, cellular_responsiveness = self.calculate_activity_constants_for_basal(basal_activity, integrated_activity_for_step_increase, cellular_activity, time)
        return basal_activity, integrated_activity, cellular_responsiveness

    def calculate_activity_constants_for_basal(self, basal_activity, integrated_activity_for_step_increase, cellular_activity, time):
        """ This function calculates the activity constants integrated activity and cellular responsiveness for the
        basal PTH pulse if no injected PTH pulse is present. It first chooses one pulse of the cellular activity after
        adaptation for the calculation of the constants.

        :param basal_activity: basal activity
        :type basal_activity: float
        :param integrated_activity_for_step_increase: integrated activity for a step increase of the same magnitude
        :type integrated_activity_for_step_increase: float
        :param cellular_activity: cellular activity depending on time
        :type cellular_activity: list
        :param time: time variable
        :type time: list

        :return: integrated activity and cellular responsiveness
        :rtype: list"""
        chosen_basal_activity_pulse = []
        chosen_basal_activity_pulse_time = []
        for i in range(len(time)):
            if (self.period_for_activity_constants * self.load_case.basal_PTH_pulse.period <= time[i] <=
                    self.period_for_activity_constants * self.load_case.basal_PTH_pulse.period + self.load_case.basal_PTH_pulse.on_duration):
                chosen_basal_activity_pulse.append(cellular_activity[i])
                chosen_basal_activity_pulse_time.append(time[i])
        integrated_activity = self.calculate_integrated_activity(chosen_basal_activity_pulse, chosen_basal_activity_pulse_time, basal_activity)
        cellular_responsiveness = self.calculate_cellular_responsiveness(integrated_activity, integrated_activity_for_step_increase, self.load_case.basal_PTH_pulse.period)
        return integrated_activity, cellular_responsiveness

    def calculate_activity_constants_for_basal_and_injection(self, basal_activity, integrated_activity_for_step_increase, cellular_activity, time):
        """ This function calculates the activity constants integrated activity and cellular responsiveness for the basal
        and injected PTH pulse. It first chooses one pulse of the cellular activity after adaptation for both basal and
        injected for the calculation of the constants. The constants are calculated separately for basal and injected PTH pulses and then summed up.

        :param basal_activity: basal activity
        :type basal_activity: float
        :param integrated_activity_for_step_increase: integrated activity for a step increase of the same magnitude
        :type integrated_activity_for_step_increase: float
        :param cellular_activity: cellular activity depending on time
        :type cellular_activity: list
        :param time: time variable
        :type time: list

        :return: integrated activity and cellular responsiveness
        :rtype: list
        """
        chosen_basal_activity_pulse = []
        chosen_basal_activity_pulse_time = []
        chosen_injected_activity_pulse = []
        chosen_injected_activity_pulse_time = []
        # find basal activity in one pulse
        for i in range(len(time)):
            if (time[-1] - 0.5 * self.load_case.basal_PTH_pulse.off_duration - (time[-1] - 0.5 * self.load_case.basal_PTH_pulse.off_duration) % self.load_case.basal_PTH_pulse.period
                    <= time[i] <=
                    time[-1] - 0.5 * self.load_case.basal_PTH_pulse.off_duration - (time[-1] - 0.5 * self.load_case.basal_PTH_pulse.off_duration) % self.load_case.basal_PTH_pulse.period + self.load_case.basal_PTH_pulse.on_duration):
                chosen_basal_activity_pulse.append(cellular_activity[i])
                chosen_basal_activity_pulse_time.append(time[i])
            if self.period_for_activity_constants * self.load_case.injected_PTH_pulse.period * 60 <= time[i] <= self.period_for_activity_constants * self.load_case.injected_PTH_pulse.period * 60 + self.load_case.injected_PTH_pulse.on_duration:
                chosen_injected_activity_pulse.append(cellular_activity[i])
                chosen_injected_activity_pulse_time.append(time[i])

        integrated_activity_for_injection = self.calculate_integrated_activity(chosen_injected_activity_pulse, chosen_injected_activity_pulse_time, basal_activity)
        integrated_activity_for_step_increase_for_injection = self.calculate_integrated_activity_for_step_increase(self.load_case.basal_PTH_pulse.min + self.load_case.injected_PTH_pulse.max)
        cellular_responsiveness_for_injection = self.calculate_cellular_responsiveness(integrated_activity_for_injection, integrated_activity_for_step_increase_for_injection, self.load_case.injected_PTH_pulse.period*60)

        basal_integrated_activity = self.calculate_integrated_activity(chosen_basal_activity_pulse, chosen_basal_activity_pulse_time, basal_activity)
        basal_cellular_responsiveness = self.calculate_cellular_responsiveness(basal_integrated_activity, integrated_activity_for_step_increase, self.load_case.basal_PTH_pulse.period)

        integrated_activity = basal_integrated_activity + integrated_activity_for_injection
        cellular_responsiveness = basal_cellular_responsiveness + cellular_responsiveness_for_injection
        return integrated_activity, cellular_responsiveness

    def calculate_integrated_activity(self, chosen_activity_pulse, chosen_activity_pulse_time, basal_activity):
        """ This function numerically calculates the integrated activity of one pulse of the cellular activity above baseline.
        It corresponds to the term alpha_T in equation (26) in the Li & Goldbeter paper.

        :param chosen_activity_pulse: cellular activity of one pulse after adaptation
        :type chosen_activity_pulse: list
        :param chosen_activity_pulse_time: time variable for the cellular activity
        :type chosen_activity_pulse_time: list
        :param basal_activity: basal activity
        :type basal_activity: float

        :return: integrated activity
        :rtype: float"""
        integrated_activity = np.linalg.norm(
            np.trapz(np.array(chosen_activity_pulse) - basal_activity, chosen_activity_pulse_time, axis=0))
        return integrated_activity

    def calculate_cellular_responsiveness(self, integrated_activity, integrated_activity_for_step_increase, period):
        """ This function calculates the cellular responsiveness of the cellular activity depending on the integrated activity.
        It corresponds to the term alpha_R in equation (31) in the Li & Goldbeter paper.

        :param integrated_activity: integrated activity of the cellular activity
        :type integrated_activity: float
        :param integrated_activity_for_step_increase: integrated activity for a step increase of the same magnitude
        :type integrated_activity_for_step_increase: float
        :param period: period of one PTH pulse
        :type period: float

        :return: cellular responsiveness
        :rtype: float"""
        cellular_responsiveness = (integrated_activity / integrated_activity_for_step_increase) * (integrated_activity / period)
        return cellular_responsiveness

    def calculate_basal_activity(self):
        """ This function calculates the basal activity of the cellular activity.
        It corresponds to the term alpha_0 in equation (9) in the Li & Goldbeter paper.

        :return: basal activity
        :rtype: float"""
        basal_activity = 1 / (1 + self.parameters.kinematics.receptor) * (
                self.parameters.activity.active_receptor * self.parameters.kinematics.receptor + self.parameters.activity.inactive_receptor)
        return basal_activity

    def calculate_integrated_activity_for_step_increase(self, stimulus_concentration):
        """ This function calculates the integrated activity for a step increase of same magnitude in the cellular activity.
        It corresponds to the term alpha_Tstep in equation (28) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: integrated activity for step increase
        :rtype: float"""
        difference_in_receptor_fraction = self.calculate_difference_in_receptor_fraction(stimulus_concentration)
        difference_in_weights = self.calculate_difference_in_weights(stimulus_concentration)
        adaptation_time = self.calculate_adaptation_time(stimulus_concentration)

        amplitude_of_cellular_activity = difference_in_receptor_fraction * difference_in_weights
        integrated_activity_for_step_increase = amplitude_of_cellular_activity * adaptation_time
        return integrated_activity_for_step_increase

    def calculate_difference_in_receptor_fraction(self, stimulus_concentration):
        """ This function calculates the difference in receptor fraction for a stimulus concentration.
        It corresponds to the term Q in equation (20) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: difference in receptor fraction for minimum and maximum of PTH pulse
        :rtype: float"""
        # calculate receptor desensitisation/ resensitisation constants for basal PTH pulse minimum
        desensitised_receptors_after_adaptation_basal_min = self.calculate_desensitised_receptors_after_adaptation(
            self.load_case.basal_PTH_pulse.min)
        # calculate receptor desensitisation/ resensitisation constants for basal PTH pulse maximum
        desensitised_receptors_after_adaptation_basal_max = self.calculate_desensitised_receptors_after_adaptation(
            stimulus_concentration)
        # calculate difference between fraction of receptors for min and max PTH pulse (eq 20)
        difference_in_receptor_fraction = desensitised_receptors_after_adaptation_basal_max - desensitised_receptors_after_adaptation_basal_min
        return difference_in_receptor_fraction

    def calculate_adaptation_time(self, stimulus_concentration):
        """ This function calculates the adaptation time of the cellular activity to a stimulus depending on desensitised
        and resensitised receptors. It corresponds to the term tau_a in equation (11) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: adaptation time
        :rtype: float"""
        contribution_of_receptor_desensitation_basal_max = self.calculate_contribution_of_receptor_desensitation(
            stimulus_concentration)
        contribution_of_receptor_resensitisation_basal_max = self.calculate_contribution_of_receptor_resensitisation(
            stimulus_concentration)
        adaptation_time = (1 / (
                contribution_of_receptor_desensitation_basal_max + contribution_of_receptor_resensitisation_basal_max))
        return adaptation_time

    def calculate_difference_in_weights(self, stimulus_concentration):
        """ This function calculates the difference in weights of active and desensitised receptors for a stimulus concentration.
        It corresponds to the term P in equation (10) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: difference in weights of active and desensitised receptors
        :rtype: float"""
        weight_active_receptor = self.calculate_weight_active_receptor(stimulus_concentration)
        weight_desensitised_receptor = self.calculate_weight_desensitised_receptor(stimulus_concentration)
        difference_in_weights = weight_active_receptor - weight_desensitised_receptor
        return difference_in_weights

    def calculate_contribution_of_receptor_desensitation(self, stimulus_concentration):
        """ This function calculates the contribution of receptor desensitisation to the cellular activity.
        It corresponds to the term u in equation (12a) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: contribution of receptor desensitisation
        :rtype: float"""
        contribution_of_receptor_desensitation_basal_min = (self.parameters.kinematics.receptor_desensitized +
                                                            self.parameters.kinematics.complex_desensitized *
                                                            stimulus_concentration) / (1 + stimulus_concentration)
        return contribution_of_receptor_desensitation_basal_min

    def calculate_contribution_of_receptor_resensitisation(self, stimulus_concentration):
        """ This function calculates the contribution of receptor resensitisation to the cellular activity.
        It corresponds to the term v in equation (12b) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: contribution of receptor resensitisation
        :rtype: float"""
        kinetic_constant = self.parameters.kinematics.receptor / self.parameters.kinematics.complex
        contribution_of_receptor_resensitisation_basal_min = ((self.parameters.kinematics.receptor_resensitized +
                                                               self.parameters.kinematics.complex_resensitized *
                                                               stimulus_concentration * kinetic_constant) /
                                                              (1 + stimulus_concentration * kinetic_constant))
        return contribution_of_receptor_resensitisation_basal_min

    def calculate_desensitised_receptors_after_adaptation(self, stimulus_concentration):
        """ This function calculates the total fraction of desensitised receptors after adaptation to a stimulus.
        It depends on the contribution of receptor desensitisation and resensitisation. It corresponds to the term Ds
        in equation (8) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: desensitised receptors after adaptation
        :rtype: float"""
        contribution_of_receptor_desensitation = self.calculate_contribution_of_receptor_desensitation(stimulus_concentration)
        contribution_of_receptor_resensitisation = self.calculate_contribution_of_receptor_resensitisation(stimulus_concentration)
        desensitised_receptors_after_adaptation = (contribution_of_receptor_desensitation /
                                                   (contribution_of_receptor_desensitation + contribution_of_receptor_resensitisation))
        return desensitised_receptors_after_adaptation

    def calculate_weight_active_receptor(self, stimulus_concentration):
        """ This function calculates the weight of active receptor in the cellular activity.
        It corresponds to the term a in equation (13a) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: weight of active receptor
        :rtype: float"""
        weight_active_receptor = ((self.parameters.activity.active_receptor + self.parameters.activity.active_complex * stimulus_concentration)
                                  / (1 + stimulus_concentration))
        return weight_active_receptor

    def calculate_weight_desensitised_receptor(self, stimulus_concentration):
        """ This function calculates the weight of desensitised receptor in the cellular activity.
        It corresponds to the term b in equation (13b) in the Li & Goldbeter paper.

        :param stimulus_concentration: stimulus/ ligand concentration
        :type stimulus_concentration: float

        :return: weight of desensitised receptor
        :rtype: float"""
        kinetic_constant = self.parameters.kinematics.receptor / self.parameters.kinematics.complex
        weight_desensitised_receptor = (self.parameters.activity.inactive_receptor + self.parameters.activity.inactive_complex
                                        * stimulus_concentration * kinetic_constant) / (1 + stimulus_concentration * kinetic_constant)
        return weight_desensitised_receptor
