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
    def __init__(self, load_case, model_type='cellular responsiveness'):
        super().__init__(load_case.lemaire)
        self.parameters = Parameters()
        self.model_type = model_type

        martonova_healthy_model = Martonova_Model(Healthy())
        _, _, _, self.parameters.healthy_integrated_activity, self.parameters.healthy_cellular_responsiveness = (
            martonova_healthy_model.solve_for_activity())
        martonova_disease_model = Martonova_Model(load_case.martonova)
        _, _, _, self.parameters.disease_integrated_activity, self.parameters.disease_cellular_responsiveness = (
            martonova_disease_model.solve_for_activity())

    def calculate_PTH_activation_OB(self, t):
        if self.model_type == 'cellular responsiveness':
            if (t is not None) and self.load_case.start_time <= t <= self.load_case.end_time:
                PTH_activation = self.parameters.calibration_parameters.cellular_responsiveness * self.parameters.disease_cellular_responsiveness
            else:
                PTH_activation = self.parameters.calibration_parameters.cellular_responsiveness * self.parameters.healthy_cellular_responsiveness
        elif self.model_type == 'integrated activity':
            if (t is not None) and self.load_case.start_time <= t <= self.load_case.end_time:
                PTH_activation = self.parameters.calibration_parameters.integrated_activity * self.parameters.disease_integrated_activity
            else:
                PTH_activation = self.parameters.calibration_parameters.integrated_activity * self.parameters.healthy_integrated_activity
        return PTH_activation


def identify_calibration_parameters():
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
    error = np.sum((lemaire_activation_PTH - parameter * martonova_activation_PTH) ** 2)
    return error


def calculate_elevation_parameter(disease):
    healthy = Healthy()
    min_healthy = healthy.basal_PTH_pulse.min
    max_healthy = healthy.basal_PTH_pulse.max
    min_disease = disease.basal_PTH_pulse.min
    max_disease = disease.basal_PTH_pulse.max
    return (min_disease + max_disease) / (min_healthy + max_healthy)
