from bone_models.bone_cell_population_models.models.lerebours_model import Lerebours_Model
from bone_models.bone_cell_population_models.load_cases.lerebours_load_cases import Lerebours_Load_Case
import matplotlib.pyplot as plt
import numpy as np

# define the time span for the simulation
tspan = [0, 200]
# initialise the load case (underloading scenario)
load_case = Lerebours_Load_Case()
# define porosity between 0 and 1
porosity = 0.05
# initialise the model instance
model = Lerebours_Model(load_case, porosity)
# solve the model for the given porosity and time span
solution = model.solve_bone_cell_population_model(tspan, porosity)
# time solution contains OBp, OBa, OCp, OCa, vascular porosity, bone volume fraction

# 1. Bone Cell Population Dynamics
plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='OBp')
plt.plot(solution.t, solution.y[1], label='OBa')
plt.plot(solution.t, solution.y[3], label='OCa')
# Force scientific notation on Y axis
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Time [days]')
plt.ylabel('Cell Concentrations [pM]')
plt.title(f'Bone Cell Population Dynamics (Porosity={porosity})')
plt.legend()
plt.show()

# 2. Volume Fractions
fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(solution.t, solution.y[4], label='Vascular Porosity', color='tab:blue')
ax1.set_xlabel('Time [days]')
ax1.set_ylabel('Vascular Porosity', color='tab:blue')
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2 = ax1.twinx()
ax2.plot(solution.t, solution.y[5], label='Bone Volume Fraction', color='tab:orange')
ax2.set_ylabel('Bone Volume Fraction', color='tab:orange')
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title(f'Vascular and Bone Volume Fractions (Porosity={porosity})')
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
plt.show()

plt.figure()
plt.plot(np.linspace(5, 95, 100),
         [model.calculate_strain_energy_density(0, 0, 1 - x, x, t=0) for x in np.linspace(5, 95, 100)])
plt.xlabel('Bone Volume Fraction')
plt.ylabel('Strain Energy Density [GPa]')
plt.title('Strain Energy Density vs Porosity')

plt.figure()
plt.plot(np.linspace(0, 1000, 100),
         [model.calculate_strain_effect_on_OBp(0, 0, 5, 95, t) for t in np.linspace(0, 1000, 100)])
plt.xlabel('Time')
plt.ylabel('Strain effect on OBp')
plt.show()
