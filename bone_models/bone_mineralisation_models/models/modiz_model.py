import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from bone_models.bone_cell_population_models.models.legacy.lerebours_model_v1 import Lerebours_Model_V1
from bone_models.bone_cell_population_models.models.martinez_reina_model import Martinez_Reina_Model
from bone_models.bone_mineralisation_models.parameters.modiz_parameters import Modiz_Parameters


class Modiz_Model(Martinez_Reina_Model):
    """
    Implementation of the Modiz model for explicit bone mineralisation.

    Based on the work of Modiz et al. (2026), this model simulates the
    mineralisation process in bone tissue as a function of porosity. It
    integrates bone ageing, volume fraction dynamics, turnover rates, and
    specific surface area to track mineral density over time.

    .. note::

        **Source Publication**:
        Corinna Modiz, Natalia M. Castoldi, Jose L. Calvo-Gallego, Stefan Scheiner, Vittorio Sansalone, Saulo Martelli, Javier Martínez-Reina, Peter Pivonka,
        *Factors of mineralization dynamics — A mathematical model for microstructural bone remodeling*
        **Bone**, 209, 117905.
        :doi:`10.1016/j.bone.2026.117905`

    .. note::
        This model is a specialized subclass of :ref:`Martinez_Reina_Model`
        and requires a :ref:`Lerebours_Model_V1` instance to handle
        cell population steady-state calculations.

    Key Components:
        - Ageing Queue: Tracks the temporal evolution of bone patches.
        - Porosity Dependence: Mineralisation kinetics scale with bone volume fraction.
        - Hypothesis Framework: Supports multiple mechanistic assumptions (1, 2.1, 2.2, 3)
          regarding resorption bias and mineral apposition.
    """

    def __init__(self, load_case, porosity, specific_surface_multiplier=1, hypothesis=3):
        """
        Initialise the Modiz model with biological and geometric parameters.

        :param load_case: Load case defining steady-state cell concentrations.
        :type load_case: Lerebours_Load_Case
        :param porosity: Initial porosity of the bone (percentage, e.g., 20.0).
        :type porosity: float
        :param specific_surface_multiplier: Scaling factor for the specific surface area,
                                            used for dispersion/boundary curves.
        :type specific_surface_multiplier: float
        :param hypothesis: The model hypothesis to apply (1 (turnover), 2.1 (age bias), 2.2 (reduced resorption), or 3 (mineralization kinetics)).
        :type hypothesis: float
        """
        super().__init__(load_case)
        self.lerebours_model = Lerebours_Model_V1(load_case, porosity / 100, specific_surface_multiplier)
        self.parameters = Modiz_Parameters()
        self.hypothesis = hypothesis

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
        pass

    def update_ageing_queue(self, OBa, OCa, bone_volume_fraction):
        """
        Update the bone ageing queue for a single time step.
        This method shifts the queue applies resorption and formation rates. The update follows these steps:

        1. **Shift & Formation**: All bone patches age by one day; new bone
           (formed fraction) is added to the start of the queue.
        2. **Resorption**: Existing bone volume is reduced across the queue
           based on the active osteoclast (OCa) concentration.
        3. **Hypothesis Handling**:
            - If ``hypothesis 2.1``: Resorption rate is biased by age.
            - If ``hypothesis 2.2``: Resorption rate is scaled with the bone volume fraction.
            - Otherwise: Resorption follows standard lag-time rules.
        4. **Mass Balance**: The final queue element acts as a residual to ensure
           the total sum matches the current ``bone_volume_fraction``.

        :param OBa: Active osteoblast cell concentration.
        :type OBa: float
        :param OCa: Active osteoclast cell concentration.
        :type OCa: float
        :param bone_volume_fraction: Current fraction of bone volume (0.0 to 1.0).
        :type bone_volume_fraction: float

        :notes: This method modifies ``self.ageing_queue`` in-place.
        """
        resorbed_bone_fraction = OCa * self.parameters.bone_volume.resorption_rate / 100
        formed_bone_fraction = OBa * self.parameters.bone_volume.formation_rate / 100
        if self.hypothesis == 2.2:
            resorption_rate = self.parameters.bone_volume.resorption_rate / 100
            if bone_volume_fraction <= 0.3:
                resorption_rate = resorption_rate / 2
            elif bone_volume_fraction < 0.8:
                resorption_rate = resorption_rate * (bone_volume_fraction + 0.2)
            else:
                resorption_rate = resorption_rate
            resorbed_bone_fraction = OCa * resorption_rate
        threshold = 1e-13
        temp_queue = self.ageing_queue.copy()
        scaling_factor = 1 - resorbed_bone_fraction / bone_volume_fraction
        if self.hypothesis == 2.1:
            resorption_bias_age = int(self.parameters.mineralisation.resorption_bias_age)
            self.ageing_queue[1:resorption_bias_age + 1] = temp_queue[0:resorption_bias_age]
            self.ageing_queue[resorption_bias_age + 1:-1] = temp_queue[resorption_bias_age:-2] * scaling_factor
            self.ageing_queue[self.ageing_queue < threshold] = 0
        else:
            lag_time = int(self.parameters.mineralisation.lag_time)
            self.ageing_queue[1:lag_time + 1] = temp_queue[0:lag_time]
            self.ageing_queue[lag_time + 1:-1] = temp_queue[lag_time:-2] * scaling_factor
            self.ageing_queue[self.ageing_queue < threshold] = 0
        self.ageing_queue[0] = formed_bone_fraction
        summed_queue_bone_volume_content = np.sum(self.ageing_queue[:-1])
        self.ageing_queue[-1] = bone_volume_fraction + (
                    formed_bone_fraction - resorbed_bone_fraction) - summed_queue_bone_volume_content
        # Check if the last element is negative
        if self.ageing_queue[-1] < threshold:
            self.ageing_queue[-1] = 0
        pass

    def calculate_average_mineral_content(self, OBa, OCa, bone_volume_fraction):
        """ Calculate the average mineral content based on the ageing queue after initialization, current bone volume
        fraction and mineralization law.

        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param bone_volume_fraction: bone volume fraction
        :type bone_volume_fraction: float

        :return: average mineral content
        :rtype: float"""
        sum_mineral_content = 0
        for j in range(len(self.ageing_queue) - 1):
            sum_mineral_content += self.ageing_queue[j] * self.calculate_mineralisation_law(j, bone_volume_fraction)
        sum_mineral_content += self.ageing_queue[-1] * self.parameters.mineralisation.maximum_mineral_content
        average_mineral_content = sum_mineral_content / bone_volume_fraction
        return average_mineral_content

    def calculate_mineralisation_law(self, t, bone_volume_fraction):
        """ Calculate the mineralisation law based on the age of the bone patch (queue element) and initial bone volume
        fraction. The law is used to calculate the mineral content in the bone tissue and includes lag time and
        exponential increase to the maximum mineral content. The mineral apposition rate depends on the initial bone
        volume fraction and prescribed hypothesis.

        :param t: age of the bone patch (queue element)
        :type t: int
        :param bone_volume_fraction: initial bone volume fraction
        :type bone_volume_fraction: float
        :return: mineral content
        :rtype: float"""
        if t <= self.parameters.mineralisation.lag_time:
            mineral = 0
        else:
            mineral_apposition_rate = self.calculate_mineral_apposition_rate(bone_volume_fraction)
            mineral = self.parameters.mineralisation.maximum_mineral_content * (1 - np.exp(-mineral_apposition_rate *
                                                                                           (t - self.parameters.mineralisation.lag_time)))
        return mineral

    def calculate_mineral_apposition_rate(self, bone_volume_fraction):
        """ Calculate the mineral apposition rate based on the initial bone volume fraction and prescribed hypothesis.

        :param bone_volume_fraction: initial bone volume fraction
        :type bone_volume_fraction: float
        :return: mineral apposition rate
        :rtype: float"""
        mineral_apposition_rate = self.parameters.mineralisation.reference_apposition_rate / bone_volume_fraction
        if self.hypothesis == 1 or self.hypothesis == 2.1 or self.hypothesis == 2.2:
            mineral_apposition_rate = self.parameters.mineralisation.reference_apposition_rate
        return mineral_apposition_rate

    def calculate_mineral_densities(self, bone_volume_fraction):
        """ Calculate the mineral densities based on the initial bone volume fraction and prescribed hypothesis.
        This function calculates the steady-state concentrations based on the Lerebours model, initialises the ageing
        queue and calculates the material and apparent density. Note: we use the turnover curve the other way around
        compared to the Lerebours model, so we need to put in bone volume fraction instead of porosity.

        :param bone_volume_fraction: initial bone volume fraction
        :type bone_volume_fraction: float
        :return: bone material density, bone apparent density
        :rtype: tuple of floats"""
        self.lerebours_model.calculate_steady_state(bone_volume_fraction)
        OBa_steady_state = self.lerebours_model.steady_state.OBa
        OCa_steady_state = self.lerebours_model.steady_state.OCa

        self.initialise_ageing_queue(OBa_steady_state, OCa_steady_state, bone_volume_fraction)

        volume_fraction_mineral = self.calculate_average_mineral_content(OBa_steady_state, OCa_steady_state,
                                                                         bone_volume_fraction)
        self.volume_fraction_mineral = volume_fraction_mineral
        self.ash_fraction = ((self.parameters.mineralisation.density_mineral * volume_fraction_mineral) /
                             (self.parameters.mineralisation.density_mineral * volume_fraction_mineral +
                              self.parameters.mineralisation.density_organic * self.parameters.mineralisation.volume_fraction_organic))
        bone_material_density = (1 + (self.parameters.mineralisation.density_mineral - 1) * volume_fraction_mineral +
                                 (self.parameters.mineralisation.density_organic - 1) * self.parameters.mineralisation.volume_fraction_organic)
        bone_apparent_density = bone_material_density * bone_volume_fraction
        return bone_material_density, bone_apparent_density

    def calculate_mineral_queue(self, ageing_queue, bone_volume_fraction):
        """ Calculate the mineral queue based on the ageing queue after initialization and initial bone volume fraction.
        The queue is calculated by multiplying the ageing queue with mineral content from the mineralisation law.

        :param ageing_queue: ageing queue
        :type ageing_queue: numpy array
        :param bone_volume_fraction: initial bone volume fraction
        :type bone_volume_fraction: float

        :return: mineral queue
        :rtype: numpy array"""
        mineral_queue = np.zeros(len(ageing_queue))
        for j in range(len(self.ageing_queue - 1)):
            mineral_queue[j] = self.ageing_queue[j] * self.calculate_mineralisation_law(j, bone_volume_fraction)
        mineral_queue[-1] = self.ageing_queue[-1] * self.parameters.mineralisation.maximum_mineral_content
        return mineral_queue

    def fit_specific_surface_to_turnover_data(self, show_plot=False):
        """ Fit the specific surface to the turnover data from Parfitt et al.
        The specific surface is calculated using the Lerebours model and the turnover data is fitted using a
        least-squares optimization and error function.

        :param show_plot: whether to show and save the plot of the fitted curve, data and specific surface
        :type show_plot: bool
        :return: scaling parameter
        :rtype: float """
        parfitt_cortical_bone_volume_fraction = 0.95
        parfitt_trabecular_bone_volume_fraction = 0.2
        parfitt_cortical_turnover = 0.77e-4  # 1/day
        parfitt_trabecular_turnover = 1.43e-4  # 1/day
        parfitt_bone_volume_fraction = np.array(
            [parfitt_cortical_bone_volume_fraction, parfitt_trabecular_bone_volume_fraction])
        parfitt_turnover = np.array([parfitt_cortical_turnover, parfitt_trabecular_turnover])

        initial_guess = 1e-3
        result = sc.optimize.minimize(self.calculate_error_for_turnover, initial_guess,
                                      args=(parfitt_bone_volume_fraction, parfitt_turnover))
        fitted_scaling_parameter = result.x[0]
        print("Fitted scaling parameter for specific surface to turnover data: ", fitted_scaling_parameter)
        if show_plot:
            self.plot_turnover_fit(fitted_scaling_parameter, parfitt_bone_volume_fraction, parfitt_turnover)
        return fitted_scaling_parameter

    def calculate_error_for_turnover(self, fitting_parameter, bone_volume_fraction_data, experimental_turnover_data):
        """ Calculate the predicted turnover and then the error to the experimental turnover data for the specific bone volume fraction.
        The error is calculated as the sum of squared differences between the predicted and experimental turnover data.

        :param fitting_parameter: scaling parameter for the specific surface
        :type fitting_parameter: float
        :param bone_volume_fraction_data: bone volume fraction data
        :type bone_volume_fraction_data: numpy array
        :param experimental_turnover_data: experimental turnover data
        :type experimental_turnover_data: numpy array
        :return: error
        :rtype: float"""
        predicted_turnover_data = fitting_parameter * self.lerebours_model.specific_surface(bone_volume_fraction_data)
        error = np.sum((predicted_turnover_data - experimental_turnover_data) ** 2)  # Sum of squared errors
        return error

    def plot_turnover_fit(self, fitted_scaling_parameter, parfitt_bone_volume_fraction, parfitt_turnover):
        """ Plot the fitted turnover curve along with the experimental data and specific surface.

        :param fitted_scaling_parameter: fitted scaling parameter for the specific surface
        :type fitted_scaling_parameter: float
        :param parfitt_bone_volume_fraction: bone volume fraction data from Parfitt et al.
        :type parfitt_bone_volume_fraction: numpy array
        :param parfitt_turnover: turnover data from Parfitt et al.
        :type parfitt_turnover: numpy array
        :return: None"""

        bone_volume_fraction = np.linspace(0, 1, 200)
        specific_surface_values = self.lerebours_model.specific_surface(bone_volume_fraction)
        turnover_values = fitted_scaling_parameter * specific_surface_values

        fig, ax1 = plt.subplots()
        line1, = ax1.plot(bone_volume_fraction, specific_surface_values, label=r'$S_\lambda$', linestyle='dashed',
                          # custom dash
                          color='#4a4a4aff', zorder=3, linewidth=1.7, alpha=0.6)
        ax1.set_xlabel('Bone Volume Fraction [-]')
        ax1.set_ylabel(r'Specific Surface [mm$^{-1}$]', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.set_xlim(0, 1)
        ax1.grid(True, linewidth=0.5, alpha=0.5)
        ax2 = ax1.twinx()
        line2, = ax2.plot(bone_volume_fraction, turnover_values, color='#d448a3ff', label='TR - Model',
                          linestyle='solid', alpha=0.85, zorder=2)
        scatter = ax2.scatter(parfitt_bone_volume_fraction, parfitt_turnover, color='#d448a3ff',
                              label='TR - Experimental', zorder=1)
        ax2.set_ylabel(r'Bone Turnover Rate [day$^{-1}$]', color='black')
        ax2.tick_params(axis='y', labelcolor='black')
        # Combine legends from both axes
        lines = [line1, line2, scatter]
        labels = [line.get_label() for line in lines]
        ax1.legend(lines, labels, )
        ax1.set_zorder(2)
        ax2.set_zorder(1)
        ax1.patch.set_visible(False)
        plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
        plt.tight_layout()
        plt.savefig('specific_surface_and_turnover_shared_axes.pdf', bbox_inches='tight')
        plt.show()
