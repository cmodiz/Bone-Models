from bone_models.models import Martonova_Model
from bone_models.load_cases import Martonova_Healthy, Martonova_Hyperparathyroidism_With_Drug
import matplotlib.pyplot as plt

load_case = Martonova_Healthy()
model = Martonova_Model(load_case)
cellular_activity, time = model.calculate_cellular_activity()
basal_activity, basal_integrated_activity, basal_cellular_responsiveness = model.calculate_activity_constants(cellular_activity, time)

print(f'Basal Activity: {basal_activity}')
print(f'Integrated Activity: {basal_integrated_activity}')
print(f'Cellular Responsiveness: {basal_cellular_responsiveness}')

plt.figure(figsize=(10, 6))
plt.plot(time, cellular_activity, label='Cellular Activity')
plt.xlabel('Time [min]')
plt.ylabel('Cellular Activity [-]')
plt.title('Cellular Activity Over Time')
plt.legend()
plt.show()

load_case = Martonova_Hyperparathyroidism_With_Drug()
model = Martonova_Model(load_case)
cellular_activity, time = model.calculate_cellular_activity()
basal_activity, basal_integrated_activity, basal_cellular_responsiveness = model.calculate_activity_constants(cellular_activity, time)

print(f'Basal Activity: {basal_activity}')
print(f'Integrated Activity: {basal_integrated_activity}')
print(f'Cellular Responsiveness: {basal_cellular_responsiveness}')

plt.figure(figsize=(10, 6))
plt.plot(time, cellular_activity, label='Cellular Activity')
plt.xlabel('Time [min]')
plt.ylabel('Cellular Activity [-]')
plt.title('Cellular Activity Over Time')
plt.legend()
plt.show()
