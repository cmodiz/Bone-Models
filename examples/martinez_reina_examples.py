from bone_models.models import Martinez_Reina_Model
from bone_models.load_cases.martinez_reina_load_cases import Martinez_Reina_Load_Case
import matplotlib.pyplot as plt


tspan = [0, 1000 + 2*365]
load_case = Martinez_Reina_Load_Case()
model = Martinez_Reina_Model(load_case)
# model.parameters.differentiation_rate.OCp = 10 * model.parameters.differentiation_rate.OCp
solution = model.solve_bone_cell_population_model(tspan=tspan)

# plt.figure()
# plt.plot(solution.t, model.bone_apparent_density, label='apparent density')
# plt.xlabel('Time [days]')
# plt.ylabel('Apparent Density')
# plt.title('Apparent Density Over Time')

plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='OBp')
plt.plot(solution.t, solution.y[1], label='OBa')
plt.plot(solution.t, solution.y[2], label='OCa')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Time [days]')
plt.ylabel('Cell Concentrations [pM]')
plt.title('Bone Cell Population Over Time')
plt.legend()

plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[3])
plt.xlabel('Time [days]')
plt.ylabel('Vascular pores volume fraction [%]')

plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[4], label='fbm')
plt.xlabel('Time [days]')
plt.ylabel('Bone volume fraction [%]')

plt.show()