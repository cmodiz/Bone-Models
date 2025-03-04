from bone_models.models import Lemaire_Model
from bone_models.load_cases.lemaire_load_cases import Lemaire_Load_Case_3
import matplotlib.pyplot as plt

load_case = Lemaire_Load_Case_3()
model = Lemaire_Model(load_case)

tspan = [0, 140]
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
# model.visualize.plot_OBp_OBa_OCa(solution)
