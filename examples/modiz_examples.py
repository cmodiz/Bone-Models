from models import Modiz_Model
from models import Lemaire_Model
from models.modiz_model import identify_calibration_parameters
from models.modiz_model import Reference_Lemaire_Model
from load_cases.modiz_load_cases import (Healthy_to_Glucocorticoid_Induced_Osteoporosis, Healthy_to_Hyperparathyroidism,
                                         Reference_Healthy_to_Hyperparathyroidism)
import matplotlib.pyplot as plt
import utils.plots

run_calibration = False
if run_calibration:
    identify_calibration_parameters()

tspan = [0, 140]
# solve with cellular responsiveness model
model = Modiz_Model(Healthy_to_Hyperparathyroidism())
solution_hpt_cellular_responsiveness = model.solve_bone_cell_population_model(tspan=tspan)

model = Modiz_Model(Healthy_to_Glucocorticoid_Induced_Osteoporosis())
solution_gio_cellular_responsiveness = model.solve_bone_cell_population_model(tspan=tspan)
# utils.plots.plot_bone_cell_concentrations(solution.t, solution.y[0], solution.y[1], solution.y[2], r'$k_R \cdot \alpha_R$')

# solve with integrated activity model
model = Modiz_Model(Healthy_to_Hyperparathyroidism(), model_type='integrated activity')
solution_hpt_integrated_activity = model.solve_bone_cell_population_model(tspan=tspan)

model = Modiz_Model(Healthy_to_Glucocorticoid_Induced_Osteoporosis(), model_type='integrated activity')
solution_gio_integrated_activity = model.solve_bone_cell_population_model(tspan=tspan)

