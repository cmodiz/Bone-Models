import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from parameters.lemaire_parameters import Parameters
from models.base_model import Base_Model


class Lemaire_Model(Base_Model):
    def __init__(self, load_case):
        super().__init__()
        self.parameters = Parameters()
        self.parameters.differentiation_rate.OBp = self.parameters.differentiation_rate.OBp * self.parameters.correction_factor.f0
        self.initial_guess_root = np.array([0.7734e-3, 0.7282e-3, 0.9127e-3])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None
        self.load_case = load_case

    def calculate_TGFb_activation_OBu(self, OCa, t):
        TGFb_activation_OBu = ((OCa + self.parameters.correction_factor.f0 * self.parameters.binding_constant.TGFb_OC) /
                               (OCa + self.parameters.binding_constant.TGFb_OC))
        return TGFb_activation_OBu

    def calculate_TGFb_repression_OBp(self, OCa, t):
        TGFb_repression_OBp = 1 / self.calculate_TGFb_activation_OBu(OCa, t)
        return TGFb_repression_OBp

    def calculate_TGFb_activation_OCa(self, OCa, t):
        TGFb_activation_OCp = self.calculate_TGFb_activation_OBu(OCa, t)
        return TGFb_activation_OCp

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        PTH_effect = self.parameters.production_rate.max_RANKL_per_cell * self.calculate_PTH_activation_OB(t) * OBa
        OPG_concentration = self.calculate_OPG_concentration(OBp, t)
        kinetics_RANKL_RANK = (
                self.parameters.binding_constant.RANKL_RANK / self.parameters.unbinding_constant.RANKL_RANK)
        kinetics_RANKL_OPG = (self.parameters.binding_constant.RANKL_OPG / self.parameters.unbinding_constant.RANKL_OPG)
        temp = kinetics_RANKL_RANK * (PTH_effect / (
                1 + kinetics_RANKL_RANK * self.parameters.concentration.RANK + kinetics_RANKL_OPG * OPG_concentration))
        RANKL_activation_OCp = temp * (1 + self.calculate_external_injection_RANKL(t) / self.parameters.production_rate.intrinsic_RANKL)
        return RANKL_activation_OCp

    def calculate_PTH_activation_OB(self, t):
        PTH = self.calculate_PTH_concentration(t)
        PTH_kinetic = self.parameters.unbinding_constant.PTH_OB / self.parameters.binding_constant.PTH_OB
        PTH_activation_OB = PTH / (self.calculate_external_injection_PTH(t) / self.parameters.degradation_rate.PTH + PTH_kinetic)
        return PTH_activation_OB

    def calculate_PTH_concentration(self, t):
        PTH = ((self.parameters.production_rate.intrinsic_PTH + self.calculate_external_injection_PTH(t)) /
               self.parameters.degradation_rate.PTH)
        return PTH

    def calculate_OPG_concentration(self, OBp, t):
        OPG = (1 / self.parameters.degradation_rate.OPG) * (self.parameters.production_rate.min_OPG_per_cell * OBp /
                                                    self.calculate_PTH_activation_OB(t)
                                                    + self.calculate_external_injection_OPG(t))
        return OPG

    def calculate_external_injection_OBp(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OBp_injection

    def calculate_external_injection_OBa(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OBa_injection

    def calculate_external_injection_OCa(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OCa_injection

    def calculate_external_injection_PTH(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.PTH_injection

    def calculate_external_injection_OPG(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OPG_injection

    def calculate_external_injection_RANKL(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.RANKL_injection

    def calculate_bone_volume_fraction_change(self, solution, steady_state, initial_bone_volume_fraction):
        self.parameters.bone_volume.resorption_rate = self.parameters.bone_volume.formation_rate * steady_state[1]/steady_state[2]
        bone_volume_fraction = [initial_bone_volume_fraction]
        for i in range(len(solution)):
            bone_volume_fraction.append(
                bone_volume_fraction[-1] +
                self.parameters.bone_volume.formation_rate * (solution[i][1] - steady_state[1]) -
                self.parameters.bone_volume.resorption_rate * (solution[i][2] - steady_state[2])
            )
        return bone_volume_fraction
