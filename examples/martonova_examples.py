from models import Martonova_Model
from load_cases.martonova_load_cases import Healthy

model = Martonova_Model(Healthy)
cellular_activity, time = model.calculate_cellular_activity()
basal_activity, integrated_activity_for_step_increase, basal_cellular_responsiveness = model.calculate_activity_constants(cellular_activity, time)