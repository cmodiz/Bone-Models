import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from models.lemaire_model import Lemaire_Model
from models.martonova_model import Martonova_Model
from load_cases.martonova_load_cases import Healthy, Hyperparathyroidism, Osteoporosis, Postmenopausal_Osteoporosis, \
    Hypercalcemia, Hypocalcemia, Glucocorticoid_Induced_Osteoporosis
from parameters.modiz_parameters import Parameters


class Modiz_Model(Lemaire_Model):
    """ This class defines the bone cell population model semi-coupled with the 2-state receptor PTH model by
    Modiz et al., 2025. The model extends the bone cell population model by Lemaire et al., 2004, with activation of
    osteoblasts by PTH calculated by a 2-state receptor model with pulsatile PTH (Martonova et al., 2023).
    It is a subclass of the Lemaire model (see :class:`Lemaire_Model`), inherits most methods,
    but modifies the calculate_PTH_activation_OB method. The Martonova model is used to calculate the activation by
    PTH separately for healthy and disease states. """
    def __init__(self, load_case, model_type='cellular responsiveness', calibration_type='all'):
        """ Constructor method. Initialises parent class with respective load case and sets model type and calibration.
        Asserts that model type and calibration type are valid. Calculates the activity constants for healthy and
        disease states (in load case) using the Martonova model.

        :param load_case: load case for the model, include both loadcases for Lemaire and Martonova model, see :class:`Load_Case` for details
        :type load_case: Load_Case
        :param model_type: type of the model (which activity constant is used) either 'cellular responsiveness' or 'integrated activity'
        :type model_type: str
        :param calibration_type: type of calibration (alignment of activity constants and old activation using either all
        states or only healthy state), either 'all' or 'only for healthy state'
        :type calibration_type: str
        :param parameters: instance of the Parameters class, see :class:`Parameters` for details
        :type parameters: Parameters

        :raises ValueError: If model_type is not 'cellular responsiveness' or 'integrated activity'.
        :raises ValueError: If calibration_type is not 'all' or 'only for healthy state'. """
        super().__init__(load_case.lemaire)
        assert model_type in ['cellular responsiveness', 'integrated activity'], \
            "Invalid model_type. Must be 'cellular responsiveness' or 'integrated activity'."
        assert calibration_type in ['all', 'only for healthy state'], \
            "Invalid calibration_type. Must be 'all' or 'only for healthy state'."

        self.parameters = Parameters()
        self.model_type = model_type
        self.calibration_type = calibration_type

        martonova_healthy_model = Martonova_Model(Healthy())
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
    It is used a reference model for the calibration of the Modiz model. """

    def __init__(self, load_case):
        """ Constructor method. Initialises parent class with respective load case.

        :param load_case: load case for the model, see :class:`Load_Case` for details
        :type load_case: Load_Case """
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
    cellular responsiveness and integrated activity using all states.

    :print: calibration parameters for cellular responsiveness and integrated activity """
    diseases = [Healthy(), Hyperparathyroidism(), Osteoporosis(), Postmenopausal_Osteoporosis(),
                Hypercalcemia(), Hypocalcemia(), Glucocorticoid_Induced_Osteoporosis()]

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
    # Perform optimization for cellular responsiveness
    result = minimize(objective_function, initial_guess,
                      args=(lemaire_activation_PTH_list, martonova_cellular_responsiveness,))
    calibration_parameter_cellular_responsiveness = result.x[0]
    print(result.message)
    print("Calibration parameter for cellular responsiveness:", calibration_parameter_cellular_responsiveness)

    result = minimize(objective_function, initial_guess,
                      args=(lemaire_activation_PTH_list, martonova_integrated_activity,))
    calibration_parameter_integrated_activity = result.x[0]
    print(result.message)
    print("Calibration parameter for integrated activity:", calibration_parameter_integrated_activity)
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
    error = np.sum((lemaire_activation_PTH - parameter * martonova_activation_PTH) ** 2)
    return error


def calculate_elevation_parameter(disease):
    """ Calculates the elevation/ reduction parameter for a disease state. The elevation parameter is the ratio of the
    minimum and maximum basal PTH pulse of the disease state to the minimum and maximum basal PTH pulse of the healthy state.

    :param disease: disease state
    :type disease: Load_Case

    :return: elevation parameter
    :rtype: float """
    healthy = Healthy()
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
    elevation_parameter = calculate_elevation_parameter(Healthy())
    print(elevation_parameter)

    lemaire_activation_PTH = lemaire_PTH_activation * elevation_parameter

    martonova_model = Martonova_Model(load_case=Healthy())
    cellular_activity, time = martonova_model.calculate_cellular_activity()
    _, integrated_activity, cellular_responsiveness = martonova_model.calculate_activity_constants(
        cellular_activity, time)

    calibration_parameter_cellular_responsiveness = lemaire_activation_PTH / cellular_responsiveness
    print("Calibration parameter for cellular responsiveness:", calibration_parameter_cellular_responsiveness)

    calibration_parameter_integrated_activity = lemaire_activation_PTH / integrated_activity
    print("Calibration parameter for integrated activity:", calibration_parameter_integrated_activity)
