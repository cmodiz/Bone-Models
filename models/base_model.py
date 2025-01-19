import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from parameters.base_parameters import Parameters


class Base_Model:
    def __init__(self):
        self.parameters = Parameters()
        self.initial_guess_root = np.array([None, None, None])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def bone_cell_population_model(self, x, t=None):
        OBp, OBa, OCa = x
        dOBpdt = ((self.parameters.differentiation_rate.OBu * self.calculate_TGFb_activation_OBu(OCa) - (self.parameters.differentiation_rate.OBp * self.parameters.correction_factor.f0) * OBp * self.calculate_TGFb_repression_OBp(OCa) + self.calculate_external_injection_OBp(t)) + self.apply_mechanical_effects() + self.apply_medication_effects_OBp())
        dOBadt = ((self.parameters.differentiation_rate.OBp * self.parameters.correction_factor.f0) * OBp * self.calculate_TGFb_repression_OBp(OCa) -
                  self.parameters.apoptosis_rate.OBa * OBa * self.apply_medication_effects_OBa()
                  + self.calculate_external_injection_OBa(t))
        dOCadt = (self.parameters.differentiation_rate.OCp * self.calculate_RANKL_activation_OCp(OBp,OBa,t) -
                  self.parameters.apoptosis_rate.OCa * OCa * self.calculate_TGFb_activation_OCa(OCa) +
                  self.calculate_external_injection_OCa(t))
        dxdt = [dOBpdt, dOBadt, dOCadt]
        return dxdt

    def calculate_steady_state(self):
        print('Calculating steady state ...', end='')
        steady_state = root(self.bone_cell_population_model, self.initial_guess_root, tol=1e-30, method="lm", options={'xtol': 1e-30}) #tol=1e-5)
        self.steady_state.OBp = steady_state.x[0]
        self.steady_state.OBa = steady_state.x[1]
        self.steady_state.OCa = steady_state.x[2]
        print(f'done \n Steady state: {steady_state.x}')
        return steady_state.x

    def solve_bone_cell_population_model(self, tspan):
        x0 = self.calculate_steady_state()
        print('Solving bone cell population model ...', end='')
        solution = solve_ivp(lambda t, x: self.bone_cell_population_model(x, t), tspan, x0, rtol=1e-8, atol=1e-8)
        print('done')
        return solution

    def calculate_TGFb_activation_OBu(self, OCa):
        pass

    def calculate_TGFb_repression_OBp(self, OCa):
        pass

    def calculate_TGFb_activation_OCa(self, OCa):
        pass

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        pass

    def calculate_external_injection_OBp(self, t):
        pass

    def calculate_external_injection_OBa(self, t):
        pass

    def calculate_external_injection_OCa(self, t):
        pass

    def calculate_bone_volume_fraction_change(self, solution, steady_state, initial_bone_volume_fraction):
        bone_volume_fraction = [initial_bone_volume_fraction]
        for i in range(len(solution)):
            bone_volume_fraction.append(
                bone_volume_fraction[-1] +
                self.parameters.bone_volume.formation_rate * (solution[i][1] - steady_state[1]) -
                self.parameters.bone_volume.resorption_rate * (solution[i][2] - steady_state[2])
            )
        return bone_volume_fraction

    def apply_mechanical_effects(self):
        return 0

    def apply_medication_effects_OBp(self):
        return 0

    def apply_medication_effects_OBa(self):
        return 1

