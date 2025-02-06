import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from parameters.scheiner_parameters import Parameters
from models.pivonka_model import Pivonka_Model


class Scheiner_Model(Pivonka_Model):
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Parameters()
        self.initial_guess_root = np.array([6.196390627918603e-004, 5.583931899482344e-004, 8.069635262731931e-004])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def apply_mechanical_effects(self, OBp, OCa, t):
        if t is None:
            # no mechanical effects in steady-state
            return 0
        elif self.parameters.mechanics.update_OBp_proliferation_rate:
            strain_effect_on_OBp = self.calculate_strain_effect_on_OBp(t)
            self.parameters.proliferation_rate.OBp = ((self.parameters.differentiation_rate.OBu *
                                                       self.parameters.mechanics.fraction_of_OBu_differentiation_rate *
                                                       self.calculate_TGFb_activation_OBu(OCa, t)) /
                                                      (self.steady_state.OBp / strain_effect_on_OBp))
            return self.parameters.proliferation_rate.OBp * OBp * strain_effect_on_OBp
        else:
            return self.parameters.proliferation_rate.OBp * OBp * self.calculate_strain_effect_on_OBp(t)

    def calculate_strain_effect_on_OBp(self, t):
        if t <= self.load_case.start_time:
            return self.parameters.mechanics.strain_effect_on_OBp_steady_state
        else:
            return (self.parameters.mechanics.strain_effect_on_OBp_steady_state *
                    (1 + 1.25 * (self.calculate_strain_energy_density(t)
                                 - self.parameters.mechanics.strain_energy_density_steady_state)
                    / self.parameters.mechanics.strain_energy_density_steady_state))
