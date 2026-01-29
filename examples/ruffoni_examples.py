import time
from bone_models.bone_mineralisation_models.models.ruffoni_model import Ruffoni_Model

# Initialize the model with a BMDD
start_time = time.time()
model = Ruffoni_Model(simulation_time=1, number_of_grid_points=500, start='BMDD')
# Solve the model
BMDD_evolution, BV_evolution, time_points = model.solve_for_BMDD(save_interval=0.3)
end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")
# Plot results
model.plot_results(BMDD_evolution, BV_evolution, time_points)

# Initialize the model with a mineralization law
start_time = time.time()
model = Ruffoni_Model(simulation_time=1, number_of_grid_points=500, start='mineralization law')
# Solve the model
BMDD_evolution, BV_evolution, time_points = model.solve_for_BMDD(save_interval=0.3)
end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")
# Plot results
model.plot_results(BMDD_evolution, BV_evolution, time_points)
