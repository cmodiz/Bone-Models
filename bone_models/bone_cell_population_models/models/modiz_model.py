import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from bone_models.bone_cell_population_models.models.lemaire_model import Lemaire_Model
from bone_models.bone_cell_population_models.models.martonova_model import Martonova_Model
from bone_models.bone_cell_population_models.load_cases.martonova_load_cases import *
from bone_models.bone_cell_population_models.parameters.modiz_parameters import Modiz_Parameters
import matplotlib as mpl


class Modiz_Model(Lemaire_Model):
    """ This class defines the bone cell population model semi-coupled with the 2-state receptor PTH model by
    Modiz et al., 2025. The model extends the bone cell population model by Lemaire et al., 2004, with activation of
    osteoblasts by PTH calculated by a 2-state receptor model with pulsatile PTH (Martonova et al., 2023).
    It is a subclass of the Lemaire model (see :class:`Lemaire_Model`), inherits most methods,
    but modifies the calculate_PTH_activation_OB method. The Martonova model is used to calculate the activation by
    PTH separately for healthy and disease states.

    .. note::
        **Source Publication:**
        Modiz C., Castoldi N.M., Scheiner S., Martinez Reina J., Gallego J. L. C., Sansalone V., Martelli S., Pivonka P. (2025).
        *Computational Simulations of Endocrine Bone Diseases Related to Pathological Glandular PTH Secretion Using a
        Multi-Scale Bone Cell Population Model*
        Journal of Applied Mathematical Modelling (submitted)

        The model is based on the following publications:
        Lemaire, V., Tobin, F. L., Greller, L. D., Cho, C. R., & Suva, L. J. (2004).
        *Modeling the interactions between osteoblast and osteoclast activities in bone remodeling.*
        Journal of Theoretical Biology, 229(3), 293-309.
        :doi:`10.1016/j.jtbi.2004.03.023`

        Martonova D., Lavaill M., Forwood M.R., Robling A., Cooper D.M.L., Leyendecker S., Pivonka P. (2023).
        *Effects of PTH glandular and external dosing patterns on bone cell activity using a two-state receptor model
        —Implications for bone disease progression and treatment*
        PLOS ONE, 18(3), e0283544.
        :doi:`10.1371/journal.pone.0283544`

    :param load_case: load case for the model, include both load cases for Lemaire and Martonova model
    :type load_case: Modiz_Load_Cases
    :param model_type: type of the model (which activity constant is used) either 'cellular responsiveness' or 'integrated activity'
    :type model_type: str
    :param calibration_type: type of calibration (alignment of activity constants and old activation using either all states or only healthy state), either 'all' or 'only for healthy state'
    :type calibration_type: str
    :param parameters: instance of the Modiz_Parameters class
    :type parameters: Modiz_Parameters

    :raises ValueError: If model_type is not 'cellular responsiveness' or 'integrated activity'.
    :raises ValueError: If calibration_type is not 'all' or 'only for healthy state'. """
    def __init__(self, load_case, model_type='cellular responsiveness', calibration_type='all'):
        """ Constructor method. Initialises parent class with respective load case and sets model type and calibration.
        Asserts that model type and calibration type are valid. Calculates the activity constants for healthy and
        disease states (in load case) using the Martonova model. """

        super().__init__(load_case.lemaire)
        assert model_type in ['cellular responsiveness', 'integrated activity'], \
            "Invalid model_type. Must be 'cellular responsiveness' or 'integrated activity'."
        assert calibration_type in ['all', 'only for healthy state'], \
            "Invalid calibration_type. Must be 'all' or 'only for healthy state'."

        self.parameters = Modiz_Parameters()
        self.model_type = model_type
        self.calibration_type = calibration_type

        martonova_healthy_model = Martonova_Model(Martonova_Healthy())
        _, _, _, self.parameters.healthy_integrated_activity, self.parameters.healthy_cellular_responsiveness = (
            martonova_healthy_model.solve_for_activity())
        martonova_disease_model = Martonova_Model(load_case.martonova)
        _, _, _, self.parameters.disease_integrated_activity, self.parameters.disease_cellular_responsiveness = (
            martonova_disease_model.solve_for_activity())

    def calculate_PTH_activation_OB(self, t):
        """ Calculates the activation of osteoblasts by PTH. The calibration parameter is selected depending on the
        calibration type. The activation is calculated using either cellular responsiveness or integrated activity
        depending on the model type and the calibration parameter. Either healthy or disease activity constants are
        returned depending on the time and load case.

        :param t: time at which the activation is calculated, if None, the activation is calculated for the steady state
        :type t: float

        :return: activation of osteoblasts by PTH
        :rtype: float """
        if self.calibration_type == 'all':
            calibration_parameter_cellular_responsiveness = self.parameters.calibration.cellular_responsiveness
            calibration_parameter_integrated_activity = self.parameters.calibration.integrated_activity
        elif self.calibration_type == 'only for healthy state':
            calibration_parameter_cellular_responsiveness = self.parameters.calibration_only_for_healthy_state.cellular_responsiveness
            calibration_parameter_integrated_activity = self.parameters.calibration_only_for_healthy_state.integrated_activity
        else:
            sys.exit("Invalid calibration type")

        if self.model_type == 'cellular responsiveness':
            if (t is not None) and self.load_case.start_time <= t <= self.load_case.end_time:
                PTH_activation = calibration_parameter_cellular_responsiveness * self.parameters.disease_cellular_responsiveness
            else:
                PTH_activation = calibration_parameter_cellular_responsiveness * self.parameters.healthy_cellular_responsiveness
        elif self.model_type == 'integrated activity':
            if (t is not None) and self.load_case.start_time <= t <= self.load_case.end_time:
                PTH_activation = calibration_parameter_integrated_activity * self.parameters.disease_integrated_activity
            else:
                PTH_activation = calibration_parameter_integrated_activity * self.parameters.healthy_integrated_activity
        return PTH_activation


class Reference_Lemaire_Model(Lemaire_Model):
    """ This class is used to modify the Lemaire model to include a disease load case with multiplicative elevation.
    It is used a reference model for the calibration of the Modiz model.

    :param load_case: load case for the model
    :type load_case: Modiz_Load_Case"""
    def __init__(self, load_case):
        """ Constructor method. Initialises parent class with respective load case.  """
        super().__init__(load_case)
        self.load_case = load_case

    def calculate_PTH_concentration(self, t):
        """ Calculates the PTH concentration. In a disease state (load case) the PTH concentration is elevated/
        decreased by a respective multiplicative factor saved as load case parameter.

        :param t: time at which the PTH concentration is calculated
        :type t: float

        :return: PTH concentration
        :rtype: float """
        if (t is not None) and self.load_case.start_time <= t <= self.load_case.end_time:
            PTH = ((self.parameters.production_rate.intrinsic_PTH * self.load_case.PTH_elevation) /
                   self.parameters.degradation_rate.PTH)
        else:
            PTH = self.parameters.production_rate.intrinsic_PTH / self.parameters.degradation_rate.PTH
        return PTH


def identify_calibration_parameters():
    """ This function identifies the calibration (elevation/decrease) parameters for the Modiz model.
    It calculates integrated activity, cellular responsiveness and old activation for healthy and disease states using
    the Martonova model and the Lemaire model. The old activation of the Lemaire model is made comparable using the
    elevation/decrease parameter. It then performs an optimization to find the calibration parameters for
    cellular responsiveness and integrated activity using all states. The calibration parameters are printed.
    The calibration parameters are then saved as parameters for the Modiz model.

    :print: calibration parameters for cellular responsiveness and integrated activity """
    diseases = [Martonova_Healthy(), Martonova_Hyperparathyroidism(), Martonova_Osteoporosis(), Martonova_Postmenopausal_Osteoporosis(),
                Martonova_Hypercalcemia(), Martonova_Hypocalcemia(), Martonova_Glucocorticoid_Induced_Osteoporosis()]

    lemaire_activation_PTH_list = []
    martonova_integrated_activity = []
    martonova_cellular_responsiveness = []
    lemaire_model = Lemaire_Model(load_case=None)
    lemaire_PTH_activation = lemaire_model.calculate_PTH_activation_OB(t=None)
    for disease in diseases:
        elevation_parameter = calculate_elevation_parameter(disease)
        print(elevation_parameter)

        lemaire_activation_PTH = lemaire_PTH_activation * elevation_parameter
        lemaire_activation_PTH_list.append(lemaire_activation_PTH)

        martonova_model = Martonova_Model(load_case=disease)
        cellular_activity, time = martonova_model.calculate_cellular_activity()
        basal_activity, integrated_activity, cellular_responsiveness = martonova_model.calculate_activity_constants(
            cellular_activity, time)
        martonova_integrated_activity.append(integrated_activity)
        martonova_cellular_responsiveness.append(cellular_responsiveness)

    initial_guess = 0.1
    lemaire_activation_PTH_array = np.array(lemaire_activation_PTH_list)
    martonova_integrated_activity_array = np.array(martonova_integrated_activity)
    martonova_cellular_responsiveness_array = np.array(martonova_cellular_responsiveness)
    # Perform optimization for cellular responsiveness
    result = minimize(objective_function, initial_guess,
                      args=(lemaire_activation_PTH_list, martonova_cellular_responsiveness,))
    calibration_parameter_cellular_responsiveness = result.x[0]
    # Calculate relRMSE
    rmse = (1 / len(lemaire_activation_PTH_array) * np.sum((lemaire_activation_PTH_array - calibration_parameter_cellular_responsiveness * martonova_cellular_responsiveness_array)**2))**(1/2)
    rel_rmse = rmse/ (np.max(lemaire_activation_PTH_list)-np.min(lemaire_activation_PTH_list))
    print(result.message)
    print(result)
    print("Calibration parameter for cellular responsiveness:", calibration_parameter_cellular_responsiveness)
    print("Absolute squared error for cellular responsiveness:", np.abs(result.fun))
    print("RMSE for cellular responsiveness:", rmse)
    print("Relative RMSE for cellular responsiveness:", rel_rmse)

    result = minimize(objective_function, initial_guess,
                      args=(lemaire_activation_PTH_list, martonova_integrated_activity,))
    calibration_parameter_integrated_activity = result.x[0]
    # Calculate relRMSE
    rmse = (1 / len(lemaire_activation_PTH_array) * np.sum((lemaire_activation_PTH_array - calibration_parameter_integrated_activity * martonova_integrated_activity_array)**2))**(1/2)
    rel_rmse = rmse / (np.max(lemaire_activation_PTH_list)-np.min(lemaire_activation_PTH_list))
    print(result.message)
    print(result)
    print("Calibration parameter for integrated activity:", calibration_parameter_integrated_activity)
    print("Absolute squared error for integrated activity:", np.abs(result.fun))
    print("RMSE for integrated activity:", rmse)
    print("Relative RMSE for integrated activity:", rel_rmse)
    pass


def objective_function(parameter, lemaire_activation_PTH, martonova_activation_PTH):
    """ Objective function for the optimization to find the calibration parameters for cellular responsiveness and
    integrated activity.

    :param parameter: calibration parameter
    :type parameter: float
    :param lemaire_activation_PTH: activation of osteoblasts by PTH calculated by the Lemaire model
    :type lemaire_activation_PTH: list
    :param martonova_activation_PTH: activation of osteoblasts by PTH calculated by the Martonova model
    :type martonova_activation_PTH: list

    :return: error
    :rtype: float """
    error = np.sum((lemaire_activation_PTH - parameter * martonova_activation_PTH)**2)
    return error


def calculate_elevation_parameter(disease):
    """ Calculates the elevation/ reduction parameter for a disease state. The elevation parameter is the ratio of the
    minimum and maximum basal PTH pulse of the disease state to the minimum and maximum basal PTH pulse of the healthy state.

    :param disease: disease state
    :type disease: Load_Case

    :return: elevation parameter
    :rtype: float """
    healthy = Martonova_Healthy()
    min_healthy = healthy.basal_PTH_pulse.min
    max_healthy = healthy.basal_PTH_pulse.max
    min_disease = disease.basal_PTH_pulse.min
    max_disease = disease.basal_PTH_pulse.max
    return (min_disease + max_disease) / (min_healthy + max_healthy)


def identify_calibration_parameters_only_for_healthy_state():
    """ This function identifies the calibration/alignment (elevation/decrease) parameters for the Modiz model. It calculates
    integrated activity, cellular responsiveness and old activation for only healthy state using the Martonova
    model and the Lemaire model. The parameter is then calculated for the healthy state only.

    :print: calibration parameters for cellular responsiveness and integrated activity """
    lemaire_model = Lemaire_Model(load_case=None)
    lemaire_PTH_activation = lemaire_model.calculate_PTH_activation_OB(t=None)
    elevation_parameter = calculate_elevation_parameter(Martonova_Healthy())
    print(elevation_parameter)

    lemaire_activation_PTH = lemaire_PTH_activation * elevation_parameter

    martonova_model = Martonova_Model(load_case=Martonova_Healthy())
    cellular_activity, time = martonova_model.calculate_cellular_activity()
    _, integrated_activity, cellular_responsiveness = martonova_model.calculate_activity_constants(
        cellular_activity, time)

    calibration_parameter_cellular_responsiveness = lemaire_activation_PTH / cellular_responsiveness
    print("Calibration parameter for cellular responsiveness:", calibration_parameter_cellular_responsiveness)

    calibration_parameter_integrated_activity = lemaire_activation_PTH / integrated_activity
    print("Calibration parameter for integrated activity:", calibration_parameter_integrated_activity)


def analyse_effect_of_different_pulse_characteristics(plot=True):
    """ This function analyses the effect of different pulse characteristics on the integrated activity and cellular responsiveness.
    It calculates the integrated activity and cellular responsiveness for healthy and disease states using the Martonova
    model for different pulse characteristics. The pulse characteristics are varied by changing the on duration of the
    pulse on adn off phases. The integrated activity and cellular responsiveness are then plotted against the log of the
    ratio of the on and off duration of the pulse. The integrated activity and cellular responsiveness for the disease
    states are also plotted as a comparison.

    :param plot: if True, the results are plotted
    :type plot: bool
    :return: -
    :rtype: - """
    diseases = [Martonova_Healthy(), Martonova_Hyperparathyroidism(), Martonova_Osteoporosis(), Martonova_Postmenopausal_Osteoporosis(),
                Martonova_Hypercalcemia(), Martonova_Hypocalcemia(), Martonova_Glucocorticoid_Induced_Osteoporosis()]
    integrated_activity_list = []
    cellular_responsiveness_list = []
    integrated_activity_list_for_disease = []
    cellular_responsiveness_list_for_disease = []
    log_on_off_phase_disease_list = []

    load_case = Martonova_Healthy()
    model = Martonova_Model(load_case)
    period = model.load_case.basal_PTH_pulse.period
    on_duration_of_pulse_list = np.linspace(0, period, num=30)

    for on_duration_of_pulse in on_duration_of_pulse_list:
        model.load_case.basal_PTH_pulse.on_duration = on_duration_of_pulse
        model.load_case.basal_PTH_pulse.off_duration = period - on_duration_of_pulse
        cellular_activity, time = model.calculate_cellular_activity()
        _, integrated_activity, cellular_responsiveness = model.calculate_activity_constants(cellular_activity, time)
        integrated_activity_list.append(integrated_activity)
        cellular_responsiveness_list.append(cellular_responsiveness)

    for disease in diseases:
        model = Martonova_Model(load_case=disease)
        cellular_activity, time = model.calculate_cellular_activity()
        _, integrated_activity, cellular_responsiveness = model.calculate_activity_constants(
            cellular_activity, time)
        integrated_activity_list_for_disease.append(integrated_activity)
        cellular_responsiveness_list_for_disease.append(cellular_responsiveness)
        log_on_off_phase_disease_list.append(np.log(model.load_case.basal_PTH_pulse.on_duration /
                                                    model.load_case.basal_PTH_pulse.off_duration))

    if plot:
        colors=['k', '#59a89c', '#0b81a2', '#7E4794', '#e25759', '#595959', '#9d2c00']
        markers = ['o', 's', '^', 'D', 'P', '*', 'X']
        labels = ['Healthy', 'HPT', 'OP', 'PMO', 'HyperC', 'HypoC', 'GIO']

        plt.figure(figsize=(7, 6))
        mpl.rcParams['font.size'] = 16
        plt.plot(np.log(on_duration_of_pulse_list/(period - on_duration_of_pulse_list)), integrated_activity_list, label='Varied healthy', color='k')
        for i in range(len(log_on_off_phase_disease_list)):
            plt.scatter(log_on_off_phase_disease_list[i],
                        integrated_activity_list_for_disease[i],
                        color=colors[i],
                        marker=markers[i],
                        label=labels[i])
        plt.xlabel(r'$\log(\frac{\tau_{on}}{\tau_{off}}) [-]$')
        plt.ylabel(r'Integrated Activity $\alpha_T$ [-]')
        plt.legend()
        plt.grid(True)
        # plt.title('Integrated Activity Over Different Pulse Characteristics')
        plt.show()

        plt.figure(figsize=(7, 6))
        mpl.rcParams['font.size'] = 16
        plt.plot(np.log(on_duration_of_pulse_list/(period - on_duration_of_pulse_list)), cellular_responsiveness_list, label='Varied healthy', color='k')
        for i in range(len(log_on_off_phase_disease_list)):
            plt.scatter(log_on_off_phase_disease_list[i],
                        cellular_responsiveness_list_for_disease[i],
                        color=colors[i],
                        marker=markers[i],
                        label=labels[i])
        plt.xlabel(r'$\log(\frac{\tau_{on}}{\tau_{off}}) [-]$')
        plt.ylabel(r'Cellular Responsiveness $\alpha_R$ [-]')
        plt.legend()
        plt.grid(True)
        # plt.title('Integrated Activity Over Different Pulse Characteristics')
        plt.show()

    pass


