from bone_models.models import Modiz_Model
from bone_models.models import Lemaire_Model
from bone_models.models.modiz_model import identify_calibration_parameters, identify_calibration_parameters_only_for_healthy_state, analyse_effect_of_different_pulse_characteristics
from bone_models.models.modiz_model import Reference_Lemaire_Model
from bone_models.load_cases.modiz_load_cases import (Modiz_Healthy_to_Hyperparathyroidism, Modiz_Healthy_to_Osteoporosis, Modiz_Healthy_to_Postmenopausal_Osteoporosis, Modiz_Healthy_to_Hypercalcemia, Modiz_Healthy_to_Hypocalcemia, Modiz_Healthy_to_Glucocorticoid_Induced_Osteoporosis)
from bone_models.load_cases.modiz_load_cases import Modiz_Reference_Healthy_to_Hyperparathyroidism, Modiz_Reference_Healthy_to_Osteoporosis, Modiz_Reference_Healthy_to_Postmenopausal_Osteoporosis, Modiz_Reference_Healthy_to_Hypercalcemia, Modiz_Reference_Healthy_to_Hypocalcemia, Modiz_Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis
import matplotlib.pyplot as plt
from bone_models.utils.plots import *
import pandas as pd
import numpy as np

analyse_effect_of_different_pulse_characteristics(plot=True)

run_calibration = False
if run_calibration:
    identify_calibration_parameters()
    identify_calibration_parameters_only_for_healthy_state()

tspan = [0, 100]
Disease_Load_Cases = [Modiz_Healthy_to_Hyperparathyroidism, Modiz_Healthy_to_Osteoporosis, Modiz_Healthy_to_Postmenopausal_Osteoporosis, Modiz_Healthy_to_Hypercalcemia, Modiz_Healthy_to_Hypocalcemia, Modiz_Healthy_to_Glucocorticoid_Induced_Osteoporosis]
Reference_Load_Cases = [Modiz_Reference_Healthy_to_Hyperparathyroidism, Modiz_Reference_Healthy_to_Osteoporosis, Modiz_Reference_Healthy_to_Postmenopausal_Osteoporosis, Modiz_Reference_Healthy_to_Hypercalcemia, Modiz_Reference_Healthy_to_Hypocalcemia, Modiz_Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis]
solutions = {}
cellular_responsiveness_calibration_all = []
cellular_responsiveness_calibration_healthy = []
integrated_activity_calibration_all = []
integrated_activity_calibration_healthy = []
old_activation = []

for Disease_Load_Case in Disease_Load_Cases:
    disease_name = Disease_Load_Case.__name__
    solutions[disease_name] = {}
    solutions[disease_name]['cellular_responsiveness'] = {}
    solutions[disease_name]['integrated_activity'] = {}
    model = Modiz_Model(Disease_Load_Case(), model_type='cellular responsiveness', calibration_type='all')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['cellular_responsiveness']['calibration_type_all'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction,
                                                                                  'PTH_activation': model.parameters.calibration.cellular_responsiveness * model.parameters.disease_cellular_responsiveness}
    cellular_responsiveness_calibration_all.append(model.parameters.calibration.cellular_responsiveness * model.parameters.disease_cellular_responsiveness)

    model = Modiz_Model(Disease_Load_Case(), model_type='integrated activity', calibration_type='all')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['integrated_activity']['calibration_type_all'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction,
                                                                              'PTH_activation': model.parameters.calibration.integrated_activity * model.parameters.disease_integrated_activity}
    integrated_activity_calibration_all.append(model.parameters.calibration.integrated_activity * model.parameters.disease_integrated_activity)

    model = Modiz_Model(Disease_Load_Case(), model_type='cellular responsiveness', calibration_type='only for healthy state')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['cellular_responsiveness']['calibration_type_only_for_healthy_state'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction,
                                                                                                     'PTH_activation': model.parameters.calibration_only_for_healthy_state.cellular_responsiveness * model.parameters.disease_cellular_responsiveness}
    cellular_responsiveness_calibration_healthy.append(model.parameters.calibration_only_for_healthy_state.cellular_responsiveness * model.parameters.disease_cellular_responsiveness)

    model = Modiz_Model(Disease_Load_Case(), model_type='integrated activity', calibration_type='only for healthy state')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['integrated_activity']['calibration_type_only_for_healthy_state'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction,
                                                                                                 'PTH_activation': model.parameters.calibration_only_for_healthy_state.integrated_activity * model.parameters.disease_integrated_activity}
    integrated_activity_calibration_healthy.append(model.parameters.calibration_only_for_healthy_state.integrated_activity * model.parameters.disease_integrated_activity)

cellular_responsiveness_calibration_all.insert(0, model.parameters.calibration.cellular_responsiveness * model.parameters.healthy_cellular_responsiveness)
integrated_activity_calibration_all.insert(0, model.parameters.calibration.integrated_activity * model.parameters.healthy_integrated_activity)
cellular_responsiveness_calibration_healthy.insert(0, model.parameters.calibration_only_for_healthy_state.cellular_responsiveness * model.parameters.healthy_cellular_responsiveness)
integrated_activity_calibration_healthy.insert(0, model.parameters.calibration_only_for_healthy_state.integrated_activity * model.parameters.healthy_integrated_activity)

for Reference_Load_Case in Reference_Load_Cases:
    disease_name = Reference_Load_Case.__name__
    solutions[disease_name] = {}
    solutions[disease_name]['old_activation'] = {}
    model = Reference_Lemaire_Model(Reference_Load_Case())
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['old_activation']['calibration_type_all'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction, 'PTH_activation': model.load_case.PTH_elevation * model.calculate_PTH_activation_OB(None)}
    old_activation.append(model.calculate_PTH_activation_OB(t=50))
old_activation.insert(0, model.calculate_PTH_activation_OB(t=0))
plot_bone_volume_fractions(solutions)


#plot_PTH_activation_for_all_disease_states(cellular_responsiveness_calibration_all, integrated_activity_calibration_all, cellular_responsiveness_calibration_healthy, integrated_activity_calibration_healthy, old_activation)
# utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='cellular_responsiveness', calibration_type='calibration_type_only_for_healthy_state')
# utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='integrated_activity', calibration_type='calibration_type_only_for_healthy_state')

# plot_all_model_options(solutions, 'Healthy_to_Hyperparathyroidism', 'Reference_Healthy_to_Hyperparathyroidism')
# plot_all_model_options(solutions, 'Healthy_to_Glucocorticoid_Induced_Osteoporosis', 'Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis')
