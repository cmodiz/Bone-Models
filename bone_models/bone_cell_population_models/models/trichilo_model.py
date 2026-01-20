import numpy as np
from scipy.integrate import solve_ivp
from bone_models.bone_cell_population_models.parameters.trichilo_parameters import Trichilo_Parameters
from bone_models.bone_cell_population_models.models.pivonka_model import Pivonka_Model


class Trichilo_Model(Pivonka_Model):
    def __init__(self, parameters: Trichilo_Parameters):
        super().__init__(parameters)
        self.parameters = Trichilo_Parameters()
        # self.initial_guess_root = np.array([6.196390627918603e-004, 5.583931899482344e-004, 8.069635262731931e-004, self.parameters.bone_volume.vascular_pore_fraction, self.parameters.bone_volume.bone_fraction])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def bone_cell_population_model(self, x, t=None):
        pass
