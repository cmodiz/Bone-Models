from models import Modiz_Model
from models import Lemaire_Model
from models.modiz_model import identify_calibration_parameters, identify_calibration_parameters_only_for_healthy_state
from models.modiz_model import Reference_Lemaire_Model
from load_cases.modiz_load_cases import (Healthy_to_Hyperparathyroidism, Healthy_to_Osteoporosis, Healthy_to_Postmenopausal_Osteoporosis, Healthy_to_Hypercalcemia, Healthy_to_Hypocalcemia, Healthy_to_Glucocorticoid_Induced_Osteoporosis)
from load_cases.modiz_load_cases import Reference_Healthy_to_Hyperparathyroidism, Reference_Healthy_to_Osteoporosis, Reference_Healthy_to_Postmenopausal_Osteoporosis, Reference_Healthy_to_Hypercalcemia, Reference_Healthy_to_Hypocalcemia, Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis
import matplotlib.pyplot as plt
import utils.plots

run_calibration = False
if run_calibration:
    identify_calibration_parameters()
    identify_calibration_parameters_only_for_healthy_state()

tspan = [0, 140]
# Healthy to Hyperparathyroidism
Disease_Load_Cases = [Healthy_to_Hyperparathyroidism, Healthy_to_Osteoporosis, Healthy_to_Postmenopausal_Osteoporosis, Healthy_to_Hypercalcemia, Healthy_to_Hypocalcemia, Healthy_to_Glucocorticoid_Induced_Osteoporosis]
Reference_Load_Cases = [Reference_Healthy_to_Hyperparathyroidism, Reference_Healthy_to_Osteoporosis, Reference_Healthy_to_Postmenopausal_Osteoporosis, Reference_Healthy_to_Hypercalcemia, Reference_Healthy_to_Hypocalcemia, Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis]
solutions = {}

for Disease_Load_Case in Disease_Load_Cases:
    disease_name = Disease_Load_Case.__name__
    solutions[disease_name] = {}
    model = Modiz_Model(Disease_Load_Case(), model_type='cellular responsiveness', calibration_type='all')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['cellular_responsiveness']['calibration_type_all'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction}

    model = Modiz_Model(Disease_Load_Case(), model_type='integrated activity', calibration_type='all')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['integrated_activity']['calibration_type_all'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction}

    model = Modiz_Model(Disease_Load_Case(), model_type='cellular responsiveness', calibration_type='only for healthy state')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['cellular_responsiveness']['calibration_type_only_for_healthy_state'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction}

    model = Modiz_Model(Disease_Load_Case(), model_type='integrated activity', calibration_type='only for healthy state')
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['integrated_activity']['calibration_type_only_for_healthy_state'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction}


for Reference_Load_Case in Reference_Load_Cases:
    disease_name = Reference_Load_Case.__name__
    solutions[disease_name] = {}
    model = Reference_Lemaire_Model(Reference_Load_Case())
    solution = model.solve_bone_cell_population_model(tspan=tspan)
    bone_volume_fraction = model.calculate_bone_volume_fraction_change(solution.t, solution.y, [model.steady_state.OBp, model.steady_state.OBa, model.steady_state.OCa], 0.3)
    solutions[disease_name]['old_activation'] = {'t': solution.t, 'y': solution.y, 'bone_volume_fraction': bone_volume_fraction}

utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='cellular_responsiveness', calibration_type='all')
utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='integrated_activity', calibration_type='all')
utils.plots.plot_bone_volume_fractions(solutions, Reference_Load_Cases, model_type='old_activation', calibration_type='all')

utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='cellular_responsiveness', calibration_type='only for healthy state')
utils.plots.plot_bone_volume_fractions(solutions, Disease_Load_Cases, model_type='integrated_activity', calibration_type='only for healthy state')

utils.plots.plot_bone_cell_concentrations(solutions, 'Healthy_to_Hyperparathyroidism', 'Reference_Healthy_to_Hyperparathyroidism', calibration_type='all')


