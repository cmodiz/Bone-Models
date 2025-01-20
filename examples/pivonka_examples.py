from models import Pivonka_Model
from load_cases.pivonka_load_cases import Load_Case_1
import matplotlib.pyplot as plt

tspan = [0, 200]
load_case = Load_Case_1()
model = Pivonka_Model(load_case)
model.parameters.differentiation_rate.OCp = 10 * model.parameters.differentiation_rate.OCp
solution = model.solve_bone_cell_population_model(tspan=tspan)


plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='OBp')
plt.plot(solution.t, solution.y[1], label='OBa')
plt.plot(solution.t, solution.y[2], label='OCa')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Time')
plt.ylabel('Cell Population')
plt.title('Bone Cell Population Over Time')
plt.legend()
plt.show()

bone_volume_fraction_change = model.calculate_bone_volume_fraction_change(solution, [model.steady_state.OBp,
                                                                                     model.steady_state.OBa,
                                                                                     model.steady_state.OCa], 1)

plt.figure(figsize=(10, 6))
plt.plot(solution.t, bone_volume_fraction_change, label='Bone Volume Fraction Change')
plt.xlabel('Time')
plt.ylabel('Bone Volume Fraction')
plt.show()