import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ..parameters.martinez_reina_parameters import Martinez_Reina_Parameters
from .scheiner_model import Scheiner_Model


class Martinez_Reina_Model(Scheiner_Model):
    """ This class implements the Martinez-Reina et al., 2019 model. It is a bone cell population model that includes
    precursor, active osteoblasts and active osteoclasts, the mechanical effects of bone remodelling in line with
    Scheiner et al., 2013, the mineralisation of bone with a queuing algorithm, and the effects of denosumab treatment
    and postmenopausal osteoporosis. The mechanical model is adjusted to account for trabecular bone and a Youngs modulus
    depending on the ash fraction. The bone volume resorbed and formed then agan depends on the bone cells stimulated by
    PMO, denosumab and mechnical stimulation.

    .. note::
        **Source Publication**:
        Martinez-Reina, J., & Pivonka, P. (2019).
        * Effects of long-term treatment of denosumab on bone mineral density: insights from an in-silico model of bone
        mineralization.*
        Bone, 125, 87-95.
        :doi:`10.1016/j.bone.2019.04.022`

    :param load_case: load case of the model
    :type load_case: Martinez_Reina_Load_Case
    :param parameters: model parameters
    :type parameters: Martinez_Reina_Parameters
    :param initial_guess_root: initial guess for the root-finding algorithm for steady-state for OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction
    :type initial_guess_root: numpy.ndarray
    :param steady_state: steady state values of the model
    :type steady_state: object
    :param ageing_queue: queue for the mineralisation of bone
    :type ageing_queue: numpy.ndarray
    :param denosumab_concentration_over_time: denosumab concentration over time
    :type denosumab_concentration_over_time: numpy.ndarray
    :param time_for_denosumab: time for denosumab concentration
    :type time_for_denosumab: numpy.ndarray
    :param bone_apparent_density: apparent density of bone over time
    :type bone_apparent_density: numpy.ndarray
    :param bone_material_density: material density of bone over time
    :type bone_material_density: numpy.ndarray
    """
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Martinez_Reina_Parameters()
        self.initial_guess_root = np.array([6.196390627918603e-004, 5.583931899482344e-004, 8.069635262731931e-004,
                                            self.parameters.bone_volume.vascular_pore_fraction,
                                            self.parameters.bone_volume.bone_fraction])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None
        self.ageing_queue = None
        self.denosumab_concentration_over_time = None
        self.time_for_denosumab = None
        self.bone_apparent_density = np.array([])
        self.bone_material_density = np.array([])

    def solve_bone_cell_population_model(self, tspan):
        """ Solve the bone cell population model and volume fractions using the ODE system over a given time interval.
        The initial conditions are set to the steady-state values.
        The queuing algorithm should only be updated every day, which prohibits the evaluation in each time step of the
        solver. Thus, the bone cell population model is solved for each interval [day, day+1] and the queuing algorithm
        is updated at the end of each interval with the current formed and resorbed bone volume fraction. Additionally,
        the bone material/ apparent density are calculated.
        The return array includes the time vector, the concentrations of OBp, OBa, OCa, the vascular pore fraction and
        the bone volume fraction.

        :param tspan: time span for the ODE solver
        :type tspan: numpy.ndarray with start and end time
        :return: solution of the ODE system [time, OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction]
        :rtype: list of arrays
        """
        x0 = self.calculate_steady_state()
        time_vector = np.array([0])
        OBp_vector = np.array([x0[0]])
        OBa_vector = np.array([x0[1]])
        OCa_vector = np.array([x0[2]])
        vascular_pore_fraction_vector = np.array([x0[3]])
        bone_volume_fraction_vector = np.array([x0[4]])
        self.initialise_ageing_queue(OBa_vector[-1], OCa_vector[-1], bone_volume_fraction_vector[-1]/100)
        for time in np.arange(tspan[0], tspan[1], 1):
            solution = solve_ivp(lambda t, x: self.bone_cell_population_model(x, t), [time_vector[-1], time], [OBp_vector[-1], OBa_vector[-1], OCa_vector[-1],
                                                                                                               vascular_pore_fraction_vector[-1], bone_volume_fraction_vector[-1]])
            time_vector = np.append(time_vector, solution.t)
            OBp_vector = np.append(OBp_vector, solution.y[0])
            OBa_vector = np.append(OBa_vector, solution.y[1])
            OCa_vector = np.append(OCa_vector, solution.y[2])
            vascular_pore_fraction_vector = np.append(vascular_pore_fraction_vector, solution.y[3])
            bone_volume_fraction_vector = np.append(bone_volume_fraction_vector, solution.y[4])

            self.update_ageing_queue(OBa_vector[-1], OCa_vector[-1], bone_volume_fraction_vector[-1]/100)
            bone_material_density = (1 + (self.parameters.mineralisation.density_mineral - 1) * self.volume_fraction_mineral +
                        (self.parameters.mineralisation.density_organic - 1) * self.parameters.mineralisation.volume_fraction_organic)
            bone_apparent_density = bone_material_density * self.parameters.bone_volume.bone_fraction / 100
            self.bone_material_density = np.append(self.bone_material_density, bone_material_density)
            self.bone_apparent_density = np.append(self.bone_apparent_density, bone_apparent_density)
        return [time_vector, OBp_vector, OBa_vector, OCa_vector, vascular_pore_fraction_vector, bone_volume_fraction_vector]

    def calculate_strain_energy_density(self, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t):
        """ Calculate the strain energy density based on macroscopic stress vector, compliance matrix, strain matrix and
        stiffness tensor.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param vascular_pore_fraction: vascular pore fraction
        :type vascular_pore_fraction: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float
        :param t: time variable
        :type t: float

        :return: strain energy density
        :rtype: float"""
        stress_vector = self.calculate_macroscopic_stress_vector(t)

        compliance_matrix = self.calculate_compliance_matrix(OBa, OCa, bone_volume_fraction, t)
        strain_matrix = compliance_matrix @ stress_vector.T

        stiffness_tensor = np.linalg.inv(compliance_matrix)
        strain_energy_density = ((1 / 2) * strain_matrix.T @ stiffness_tensor @ strain_matrix)
        return strain_energy_density

    def calculate_compliance_matrix(self, OBa, OCa, bone_volume_fraction, t):
        """ Calculate the compliance matrix based on the current Young's modulus and Poisson's ratio. The Youngs modulus
        is updated before the calculation.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float
        :param t: time variable
        :type t: float

        :return: compliance matrix
        :rtype: numpy.ndarray"""
        youngs_modulus = self.calculate_youngs_modulus(OBa, OCa, bone_volume_fraction, t)
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

    def calculate_youngs_modulus(self, OBa, OCa, bone_volume_fraction, t):
        """ Calculate the Young's modulus based on the bone volume fraction and ash fraction.
        The ash fraction is updated before the calculation.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float
        :param t: time variable
        :type t: float

        :return: Young's modulus
        :rtype: float
        """
        ash_fraction = self.calculate_ash_fraction(OBa, OCa, bone_volume_fraction, t)
        youngs_modulus = 84.37 * ((bone_volume_fraction/100) ** 2.58) * (ash_fraction ** 2.74)
        return youngs_modulus

    def calculate_ash_fraction(self, OBa, OCa, bone_volume_fraction, t):
        """ Calculate the ash fraction based on the average mineral content, the volume fraction of mineral and organic
        material and the respective densities. The volume fraction of mineral is updated before the calculation.

        : param OBa: active osteoblast cell concentration
        : type OBa: float
        : param OCa: active osteoclast cell concentration
        : type OCa: float
        : param bone_volume_fraction: bone volume fraction
        : type bone_volume_fraction: float
        : param t: time variable
        : type t: float

        : return: ash fraction
        : rtype: float"""
        volume_fraction_mineral = self.calculate_average_mineral_content(OBa, OCa, bone_volume_fraction, t)
        ash_fraction = (self.parameters.mineralisation.density_mineral * volume_fraction_mineral /
                        (self.parameters.mineralisation.density_mineral * volume_fraction_mineral +
                         self.parameters.mineralisation.density_organic *
                         self.parameters.mineralisation.volume_fraction_organic))
        self.volume_fraction_mineral = volume_fraction_mineral
        return ash_fraction

    def calculate_average_mineral_content(self, OBa, OCa, bone_volume_fraction, t):
        """ Calculate the average mineral content based on the mineralisation law and the mineralisation queue.
        Each element of the queue is multiplied by the mineral content determined by mineralisation law and summed up.
        The last element is assigned the maximum mineral content.
        The average mineral content is then calculated by dividing the sum by the bone volume fraction.

        : param OBa: active osteoblast cell concentration
        : type OBa: float
        : param OCa: active osteoclast cell concentration
        : type OCa: float
        : param bone_volume_fraction: bone volume fraction
        : type bone_volume_fraction: float
        : param t: time variable
        : type t: float

        : return: average mineral content
        : rtype: float"""
        sum_mineral_content = 0
        for j in range(len(self.ageing_queue) - 1):
            sum_mineral_content += self.ageing_queue[j] * self.calculate_mineralisation_law(j)
        sum_mineral_content += self.ageing_queue[-1] * self.parameters.mineralisation.maximum_mineral_content
        average_mineral_content = sum_mineral_content / (bone_volume_fraction/100)
        return average_mineral_content

    def update_ageing_queue(self, OBa, OCa, bone_volume_fraction):
        """ Update the ageing queue based on the current resorbed and formed bone volume fraction. The queue is updated
        by shifting all elements by one position, resorb volume based on the elements size and adding the new formed bone
        volume fraction at the beginning. The last element of the queue stores the volume needed for all the elements to
        sum up to the bone volume fraction.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float

        :return: None"""
        # update mineralization queue
        resorbed_bone_fraction = OCa * self.parameters.bone_volume.resorption_rate/100
        formed_bone_fraction = OBa * self.parameters.bone_volume.formation_rate/100

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
        self.ageing_queue[-1] = bone_volume_fraction + (
                    formed_bone_fraction - resorbed_bone_fraction) - summed_queue_bone_volume_content
        if self.ageing_queue[-1] < 1e-13:
            self.ageing_queue[-1] = 0
            print('Take care, volume of last element in the queue <0')
        #print('Updated ageing queue.')
        pass

    def initialise_ageing_queue(self, OBa, OCa, bone_volume_fraction):
        """ Initialise the ageing queue with the initial bone volume fraction and resorbed/ formed fractions. It is
        initialised with zeros and the update_ageing_queue function is called until the queue is filled.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float

        :return: None"""
        self.ageing_queue = np.zeros(self.parameters.mineralisation.length_of_queue)
        j = 0
        while j < self.parameters.mineralisation.length_of_queue:
            self.update_ageing_queue(OBa, OCa, bone_volume_fraction)
            j += 1
        print('Initialised ageing queue.')
        pass

    def calculate_mineralisation_law(self, t):
        """ Calculate the mineralisation law based on the time variable. In the lag time, the mineral content is zero.
        In the primary phase, the mineral content increases linearly from zero to the primary mineral content.
        In the secondary phase, the mineral content increases exponentially to the maximum mineral content.
        The mineral content is based on the age of the patch in the queue.

        :param t: time variable
        :type t: float

        :return: mineral content
        :rtype: float"""
        if t <= self.parameters.mineralisation.lag_time:
            mineral_content = 0
        elif t <= self.parameters.mineralisation.primary_phase_duration + self.parameters.mineralisation.lag_time:
            mineral_content = self.parameters.mineralisation.primary_mineral_content * (
                    t - self.parameters.mineralisation.lag_time) / self.parameters.mineralisation.primary_phase_duration
        else:
            mineral_content = (self.parameters.mineralisation.maximum_mineral_content + (
                        self.parameters.mineralisation.primary_mineral_content
                        - self.parameters.mineralisation.maximum_mineral_content) *
                               np.exp(-self.parameters.mineralisation.rate * (
                                       t - self.parameters.mineralisation.primary_phase_duration - self.parameters.mineralisation.lag_time)))
        return mineral_content

    def calculate_RANKL_concentration(self, OBp, OBa, t):
        """ Calculates the RANKL concentration based on the effective carrying capacity, RANKL-RANK-OPG binding,
        degradation rate, intrinsic RANKL production and external injection of RANKL (PMO) and denosumab effect.
        Postmenopausal osteoporosis is simulated via an external injection of RANKL (see Eq. (18) in the source paper).
        Denosumab effect is included based on the current concentration in the compartment, binding constant of denosumab
        and RANKL and an accessibility factor.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float

        :return: RANKL concentration
        :rtype: float"""
        denosumab_effect = self.parameters.denosumab.accessibility_factor * self.parameters.binding_constant.RANKL_denosumab * self.calculate_denosumab_concentration(
            t)
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
        """ Calculates the denosumab concentration based on the current time and the pharmacokinetics-pharmacodynamics
        model. If the time is outside of the treatment period, the concentration is zero. Otherwise, the concentration
        is calculated using the solution of the PK-PD model after translating the time back to the solution time interval.
        The PK-PD model returns the concentration in ng/ml, which is then converted to pmol/L using the molar mass of denosumab.

        :param t: time variable
        :type t: float

        :return: denosumab concentration
        :rtype: float"""
        if self.denosumab_concentration_over_time is None:
            self.solve_for_denosumab_concentration()
        if t is None or t <= self.load_case.start_denosumab_treatment or t >= self.load_case.end_denosumab_treatment:
            return 0
        else:
            # calculate how many injections were given already
            number_of_injections_given = np.floor(
                (t - self.load_case.start_denosumab_treatment) / self.load_case.treatment_period)
            # translate the current time to the interval [0, treatment_period]
            time_translation = (
                        t - self.load_case.start_denosumab_treatment - number_of_injections_given * self.load_case.treatment_period)
            # find the closest index in the solution of the denosumab ODE
            closest_index = np.argmin(np.abs(self.time_for_denosumab - time_translation))
            # get C_den in ng/ml and calculate it to pmol/L using the molar mass (in ng/mol) of denosumab
            denosumab_concentration = self.denosumab_concentration_over_time[closest_index]/self.parameters.denosumab.molar_mass * (10 ** 6)
            return denosumab_concentration

    def solve_for_denosumab_concentration(self):
        """ Solve the PK-PD model for denosumab concentration over time. The initial concentration is set to zero and
        the ODe is solved over the treatment period.

        :return: None"""
        initial_denosumb_concentration = 0
        t_span = [0, self.load_case.treatment_period]
        sol = solve_ivp(self.pharmacokinetics_pharmocodynamics_denosumab, t_span, [initial_denosumb_concentration],
                        rtol=1e-10, atol=1e-10, max_step=1)
        self.denosumab_concentration_over_time = sol.y[0]  # ng/ml
        self.time_for_denosumab = sol.t
        plt.figure()
        plt.plot(sol.t, sol.y[0])
        plt.xlabel('Time [days]')
        plt.ylabel('Denosumab concentration [ng/ml]')
        # plt.show()
        print('Denosumab concentration initialized')
        pass

    def pharmacokinetics_pharmocodynamics_denosumab(self, t, concentration):
        """ The PK-PD model for denosumab concentration over time. The ODE is based on the absorption rate,
        elimination rate, volume of the central compartment, bioavailability, Michaelis-Menten constant, maximum volume
        and injected dose.

        :param t: time variable
        :type t: float
        :param concentration: injected denosumab concentration
        :type concentration: float

        :return: derivative of the denosumab concentration
        :rtype: float"""
        dCdt = ((self.parameters.denosumab.absorption_rate * (
                    self.load_case.denosumab_dose / (self.parameters.denosumab.volume_central_compartment
                                                     / self.parameters.denosumab.bioavailability) *
                    np.exp(-self.parameters.denosumab.absorption_rate * t)) -
                 (concentration / (self.parameters.denosumab.michaelis_menten_constant + concentration)) *
                 (self.parameters.denosumab.maximum_volume / (self.parameters.denosumab.volume_central_compartment /
                                                              self.parameters.denosumab.bioavailability))) -
                self.parameters.denosumab.elimination_rate * concentration)
        return dCdt

    def calculate_external_injection_RANKL(self, t):
        """ Calculate the external injection of RANKL based on the time variable for simulation of PMO.
        The injection is zero if the time is outside of the PMO simulation period. Otherwise, the injection is calculated
        based on the increase in RANKL, the reduction factor, the characteristic time and the time variable.
        The RANKL increase due to PMO is assumed to start at a fixed value and decrease over time.

        :param t: time variable
        :type t: float

        :return: external injection of RANKL (increase due to PMO)
        :rtype: float"""
        if t is None or self.load_case.start_postmenopausal_osteoporosis > t or t > self.load_case.end_postmenopausal_osteoporosis:
            external_injection_RANKL = 0
        elif self.load_case.start_postmenopausal_osteoporosis <= t <= self.load_case.end_postmenopausal_osteoporosis:
            external_injection_RANKL = self.parameters.PMO.increase_in_RANKL * (
                    self.parameters.PMO.reduction_factor ** 2 / (self.parameters.PMO.reduction_factor ** 2 +
                                                                 ((t - self.load_case.start_postmenopausal_osteoporosis)
                                                                  / self.parameters.PMO.characteristic_time) ** 2))
        else:
            external_injection_RANKL = 0
        return external_injection_RANKL
