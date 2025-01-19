from models import Lemaire_Model
from load_cases.lemaire_load_cases import Load_Case_5
import matplotlib.pyplot as plt

load_case = Load_Case_5()
model = Lemaire_Model(load_case)

tspan = [0, 140]
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
# model.visualize.plot_OBp_OBa_OCa(solution)
