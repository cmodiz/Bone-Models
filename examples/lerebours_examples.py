from bone_models.bone_cell_population_models.models import Lerebours_Model
from bone_models.bone_cell_population_models.load_cases.lerebours_load_cases import Lerebours_Load_Case
import matplotlib.pyplot as plt
import numpy as np


tspan = [0, 4000]
load_case = Lerebours_Load_Case()

OBu_list = []
OBp_list = []
OBa_list = []
OCu_list = []
OCp_list = []
OCa_list = []
turnover_list = []

for porosity in np.linspace(0, 1, 50):
    model = Lerebours_Model(load_case, porosity)
    result = model.calculate_steady_state(porosity)
    OBu_list.append(model.steady_state.OBu)
    OBp_list.append(model.steady_state.OBp)
    OBa_list.append(model.steady_state.OBa)
    OCu_list.append(model.steady_state.OCu)
    OCp_list.append(model.steady_state.OCp)
    OCa_list.append(model.steady_state.OCa)
    turnover_list.append(model.calculate_turnover(porosity))

plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, 1, 50), OBu_list, label='OBu')
plt.plot(np.linspace(0, 1, 50), OCu_list, label='OCu')
plt.xlabel('Porosity')
plt.ylabel('Cell Concentrations [pM]')
plt.title('Bone Cell Population Over Porosity')
plt.legend()

plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, 1, 50), OBa_list, label='OBa')
plt.plot(np.linspace(0, 1, 50), OCa_list, label='OCa')
plt.xlabel('Porosity')
plt.ylabel('Cell Concentrations [pM]')
plt.title('Bone Cell Population Over Porosity')
plt.legend()

plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, 1, 50), turnover_list)
plt.xlabel('Porosity')
plt.ylabel('Turnover rate [day^-1]')
plt.title('Turnover Rate vs Porosity')
plt.show()
