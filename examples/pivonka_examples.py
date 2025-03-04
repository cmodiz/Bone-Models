from bone_models.models import Pivonka_Model
from bone_models.load_cases.pivonka_load_cases import Pivonka_Load_Case_1
import matplotlib.pyplot as plt

tspan = [0, 140]
load_case = Pivonka_Load_Case_1()
model = Pivonka_Model(load_case)
# model.parameters.differentiation_rate.OCp = 10 * model.parameters.differentiation_rate.OCp
solution = model.solve_bone_cell_population_model(tspan=tspan)


plt.figure()
plt.plot(solution.t, solution.y[0], label='OBp', color='blue', linestyle='--')
plt.plot(solution.t, solution.y[1], label='OBa', color='blue')
plt.plot(solution.t, solution.y[2], label='OCa', color='red')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.grid(True)
plt.xlabel('Time [days]')
plt.ylabel('Cell Concentrations [pM]')
plt.title('Bone Cell Population Over Time')
plt.legend()
plt.show()
#
# bone_volume_fraction_change = model.calculate_bone_volume_fraction_change(solution, [model.steady_state.OBp,
#                                                                                      model.steady_state.OBa,
#                                                                                      model.steady_state.OCa], 1)
#
# plt.figure(figsize=(10, 6))
# plt.plot(solution.t, bone_volume_fraction_change, label='Bone Volume Fraction Change')
# plt.xlabel('Time')
# plt.ylabel('Bone Volume Fraction')
# plt.show()