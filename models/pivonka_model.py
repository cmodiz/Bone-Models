import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from parameters.pivonka_parameters import Parameters
from models.lemaire_model import Lemaire_Model


class Pivonka_Model(Lemaire_Model):
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Parameters()
        self.initial_guess_root = np.array([0.7734e-3, 0.7282e-3, 0.9127e-3])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def calculate_TGFb_activation_OBu(self, OCa, t):
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_activation_OBu = TGFb / (TGFb + self.parameters.activation_coefficient.TGFb_OBu)
        return TGFb_activation_OBu

    def calculate_TGFb_repression_OBp(self, OCa, t):
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_repression_OBp = self.parameters.repression_coefficient.TGFb_OBp / (
                    TGFb + self.parameters.repression_coefficient.TGFb_OBp)
        return TGFb_repression_OBp

    def calculate_TGFb_activation_OCa(self, OCa, t):
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_activation_OCp = TGFb / (TGFb + self.parameters.activation_coefficient.TGFb_OCa)
        return TGFb_activation_OCp

    def calculate_TGFb_concentration(self, OCa, t):
        TGFb = (self.parameters.bone_volume.stored_TGFb_content * self.parameters.bone_volume.resorption_rate * OCa +
                self.calculate_external_injection_TGFb(t)) / self.parameters.degradation_rate.TGFb
        return TGFb

    def calculate_external_injection_TGFb(self, t):
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.TGFb_injection

    def calculate_PTH_activation_OB(self, t):
        PTH = self.calculate_PTH_concentration(t)
        PTH_activation_OB = PTH / (PTH + self.parameters.activation_coefficient.PTH_OB)
        return PTH_activation_OB

    def calculate_PTH_repression_OB(self, t):
        PTH = self.calculate_PTH_concentration(t)
        PTH_repression_OB = self.parameters.repression_coefficient.PTH_OB / (
                    PTH + self.parameters.repression_coefficient.PTH_OB)
        return PTH_repression_OB

    def calculate_OPG_concentration(self, OBp, OBa, t):
        temp_PTH_OB = (self.parameters.production_rate.min_OPG_per_cell * OBp +
                       self.parameters.production_rate.min_OPG_per_cell * OBa) * self.calculate_PTH_repression_OB(t)
        OPG = (((temp_PTH_OB + self.calculate_external_injection_OPG(t)) * self.parameters.concentration.OPG_max) /
               (temp_PTH_OB + self.parameters.degradation_rate.OPG * self.parameters.concentration.OPG_max))
        return OPG

    def calculate_effective_carrying_capacity_RANKL(self, OBp, OBa, t):
        RANKL_eff = (self.parameters.concentration.max_RANKL_per_cell * OBp +
                     self.parameters.concentration.max_RANKL_per_cell * OBa) * self.calculate_PTH_activation_OB(t)
        return RANKL_eff

    def calculate_RANKL_concentration(self, OBp, OBa, t):
        RANKL_eff = self.calculate_effective_carrying_capacity_RANKL(OBp, OBa, t)
        RANKL_RANK_OPG = RANKL_eff / (1 + self.parameters.binding_constant.RANKL_OPG *
                                      self.calculate_OPG_concentration(OBp, OBa, t) +
                                      self.parameters.binding_constant.RANKL_RANK * self.parameters.concentration.RANK)
        RANKL = RANKL_RANK_OPG * ((self.parameters.production_rate.RANKL_rate_per_cell +
                                   self.calculate_external_injection_RANKL(t)) /
                                  self.parameters.production_rate.RANKL_rate_per_cell +
                                  self.parameters.degradation_rate.RANKL * RANKL_eff)
        return RANKL

    def calculate_RANKL_RANK_concentration(self, OBp, OBa, t):
        RANKL = self.calculate_RANKL_concentration(OBp, OBa, t)
        RANKL_RANK = self.parameters.activation_coefficient.RANKL_RANK * RANKL * self.parameters.concentration.RANK
        return RANKL_RANK

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        RANKL_RANK = self.calculate_RANKL_RANK_concentration(OBp, OBa, t)
        RANKL_activation_OCp = RANKL_RANK / (RANKL_RANK + self.parameters.activation_coefficient.RANKL_OCp)
        return RANKL_activation_OCp
