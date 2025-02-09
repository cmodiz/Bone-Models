from bone_models.models import Scheiner_Model

tspan = [0, 200]
load_case = Scheiner_Load_Case()
model = Scheiner_Model(load_case)
# model.parameters.differentiation_rate.OCp = 10 * model.parameters.differentiation_rate.OCp
solution = model.solve_bone_cell_population_model(tspan=tspan)