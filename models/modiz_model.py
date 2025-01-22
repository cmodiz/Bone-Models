import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from models.lemaire_model import Lemaire_Model
from models.martonova_model import Martonova_Model
import load_cases.martonova_load_cases as load_cases
from parameters.modiz_parameters import Parameters


class Modiz_Model(Lemaire_Model, Martonova_Model):
    def __init__(self, load_case):
        Lemaire_Model.__init__(self, load_case)
        Martonova_Model.__init__(self, load_case)
        self.parameters = Parameters()

        pass

    def calculate_PTH_activation_OB(self, t):
        cellular_activity, time = self.calculate_cellular_activity()
        basal_activity, integrated_activity, cellular_responsiveness = self.calculate_activity_constants(cellular_activity, time)
        PTH_activation_OB = self.parameters.calibration_parameters.cellular_responsiveness * cellular_responsiveness
        return PTH_activation_OB


def identify_calibration_parameters():
    disease_states = ['Healthy', 'Hyperparathyroidism', 'Osteoporosis', 'Postmenopausal_Osteoporosis', 'Hypercalcemia',
                      'Hypocalcemia', 'Glucocorticoid_Induced_Osteoporosis']

    lemaire_activation_PTH_list = []
    martonova_integrated_activity = []
    martonova_cellular_responsiveness = []
    for disease in disease_states:
        elevation_parameter = calculate_elevation_parameter(disease)

        lemaire_model = Lemaire_Model(load_case=None)
        lemaire_model.parameters.production_rate.intrinsic_PTH = lemaire_model.parameters.production_rate.intrinsic_PTH * elevation_parameter
        lemaire_activation_PTH = lemaire_model.calculate_PTH_activation_OB(t=None)
        lemaire_activation_PTH_list.append(lemaire_activation_PTH)

        load_case_class = getattr(load_cases, disease)
        martonova_model = Martonova_Model(load_case=load_case_class())
        cellular_activity, time = martonova_model.calculate_cellular_activity()
        basal_activity, integrated_activity, cellular_responsiveness = martonova_model.calculate_activity_constants(cellular_activity, time)
        martonova_integrated_activity.append(integrated_activity)
        martonova_cellular_responsiveness.append(cellular_responsiveness)

    initial_guess = 0.1
    # Perform optimization for cellular responsiveness
    result = minimize(error_function, initial_guess, args=(lemaire_activation_PTH_list, martonova_cellular_responsiveness,))
    calibration_parameter_cellular_responsiveness = result.x[0]
    print(result.message)
    print("Calibration parameter for cellular responsiveness:", calibration_parameter_cellular_responsiveness)

    result = minimize(error_function, initial_guess, args=(lemaire_activation_PTH_list, martonova_integrated_activity,))
    calibration_parameter_integrated_activity = result.x[0]
    print(result.message)
    print("Calibration parameter for integrated activity:", calibration_parameter_integrated_activity)
    pass


def error_function(parameter, list_1, list_2):
    total_error = 0
    for i in range(len(list_1)):
        total_error += (list_1[i] - parameter * list_2[i]) ** 2
    return total_error


def calculate_elevation_parameter(disease):
    min_healthy = 0.00332  # [nM]
    max_healthy = 0.00276  # [nM]
    if disease == 'Healthy':
        return 1
    elif disease == 'Hyperparathyroidism':
        min_hpt = 0.01381
        max_hpt = 0.00977
        return (min_hpt + max_hpt) / (min_healthy + max_healthy)
    elif disease == 'Osteoporosis':
        min_op = 0.003321
        max_op = 0.0016967
        return (min_op + max_op) / (min_healthy + max_healthy)
    elif disease == 'Postmenopausal_Osteoporosis':
        min_pmo = 0.00332 * 0.9
        max_pmo = 0.00276 * 0.9
        return (min_pmo + max_pmo) / (min_healthy + max_healthy)
    elif disease == 'Hypercalcemia':
        min_hyperc = 0.00332 * 0.25
        max_hyperc = 0.00276 * 0.12
        return (min_hyperc + max_hyperc) / (min_healthy + max_healthy)
    elif disease == 'Hypocalcemia':
        min_hypoc = 0.00332 * 2.57
        max_hypoc = 0.00276 * 13.1
        return (min_hypoc + max_hypoc) / (min_healthy + max_healthy)
    elif disease == 'Glucocorticoid_Induced_Osteoporosis':
        min_gio = 0.00332 * 0.48
        max_gio = 0.00276 * 1.75
        return (min_gio + max_gio) / (min_healthy + max_healthy)
