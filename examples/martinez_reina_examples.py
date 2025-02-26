from bone_models.models import Martinez_Reina_Model
from bone_models.load_cases.martinez_reina_load_cases import Martinez_Reina_Load_Case
import matplotlib.pyplot as plt
import numpy as np


tspan = [0, 4000]
load_case = Martinez_Reina_Load_Case()
model = Martinez_Reina_Model(load_case)
# model.parameters.differentiation_rate.OCp = 10 * model.parameters.differentiation_rate.OCp
solution = model.solve_bone_cell_population_model(tspan=tspan)
[time, OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction] = solution

plt.figure(figsize=(10, 6))
plt.plot(time, OBp, label='OBp')
plt.plot(time, OBa, label='OBa')
plt.plot(time, OCa, label='OCa')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Time [days]')
plt.ylabel('Cell Concentrations [pM]')
plt.title('Bone Cell Population Over Time')
plt.legend()

plt.figure(figsize=(10, 6))
plt.plot(time, bone_volume_fraction)
plt.xlabel('Time [days]')
plt.ylabel('Bone volume fraction [%]')

plt.figure(figsize=(10, 6))
plt.plot(time, vascular_pore_fraction, label='fbm')
plt.xlabel('Time [days]')
plt.ylabel('Vascular pore fraction [%]')

plt.show()

plt.figure()
plt.plot(np.arange(tspan[0], tspan[1], 1), model.bone_apparent_density, label='apparent density')
plt.xlabel('Time [days]')
plt.ylabel('Apparent Density')
plt.title('Apparent Density Over Time')

plt.figure()
plt.plot(np.arange(tspan[0], tspan[1], 1), model.bone_material_density, label='material density')
plt.xlabel('Time [days]')
plt.ylabel('Apparent Density')
plt.title('Apparent Density Over Time')
plt.show()