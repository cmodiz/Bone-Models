import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from bone_models.bone_cell_population_models.parameters.pivonka_parameters import Pivonka_Parameters
from bone_models.bone_cell_population_models.models.lemaire_model import Lemaire_Model


class Pivonka_Model(Lemaire_Model):
    """ This class implements the bone cell population model by Lemaire et al. (2004) as a subclass of the Base_Model class.

    .. note::
       **Source Publication**:
       Pivonka P., Zimak J., Smith D.W., Gardiner B.S., Dunstan C.R., Sims N.A., Martin J.T., Mundy G.R. (2008).
       *Model structure and control of bone remodeling: A theoretical study.*
       Bone, 43(2), 249-263.
       :doi:`10.1016/j.bone.2008.03.025`

    :param load_case: load case for the model
    :type load_case: object
    :param parameters: model parameters
    :type parameters: Pivonka_Parameters
    :param initial_guess_root: initial guess for the root-finding algorithm for steady-state
    :type initial_guess_root: numpy.ndarray
    :param steady_state: steady state values of the model
    :type steady_state: object
    :param t: time variable
    :type t: float
    :param tspan: time span for the ODE solver
    :type tspan: numpy.ndarray with start and end time
    :param x: state variables of the model
    :type x: list
    :param solution: solution of the ODE system
    :type solution: scipy.integrate._ivp.ivp.OdeResult
    :param initial_bone_volume_fraction: initial bone volume fraction
    :type initial_bone_volume_fraction: float
    :param dOBpdt: rate of change of precursor osteoblast cell concentration
    :type dOBpdt: float
    :param dOBadt: rate of change of active osteoblast cell concentration
    :type dOBadt: float
    :param dOCadt: rate of change of osteoclast cell concentration
    :type dOCadt: float
    :param dxdt: rate of change of state variables
    :type dxdt: list
    :param bone_volume_fraction: bone volume fraction over time
    :type bone_volume_fraction: list"""
    def __init__(self, load_case):
        super().__init__(load_case=load_case)
        self.parameters = Pivonka_Parameters()
        self.initial_guess_root = np.array([[6.196390627918603e-4, 5.583931899482344e-4, 8.069635262731931e-4]])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None
        self.is_load_case = False

    def calculate_TGFb_activation_OBu(self, OCa, t):
        """ Calculates the activation of uncommitted osteoblasts by TGF-beta based on the current TGF-beta concentration.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: activation of uncommitted osteoblasts by TGF-beta
        :rtype: float"""
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_activation_OBu = TGFb / (TGFb + self.parameters.activation_coefficient.TGFb_OBu)
        return TGFb_activation_OBu

    def calculate_TGFb_repression_OBp(self, OCa, t):
        """ Calculates the repression of precursor osteoblasts by TGF-beta based on the current TGF-beta concentration.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: repression of precursor osteoblasts by TGF-beta
        :rtype: float"""
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_repression_OBp = self.parameters.repression_coefficient.TGFb_OBp / (
                    TGFb + self.parameters.repression_coefficient.TGFb_OBp)
        return TGFb_repression_OBp

    def calculate_TGFb_activation_OCa(self, OCa, t):
        """ Calculates the activation of active osteoclasts by TGF-beta based on the current TGF-beta concentration.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: activation of active osteoclasts by TGF-beta
        :rtype: float"""
        TGFb = self.calculate_TGFb_concentration(OCa, t)
        TGFb_activation_OCp = TGFb / (TGFb + self.parameters.activation_coefficient.TGFb_OCa)
        return TGFb_activation_OCp

    def calculate_TGFb_concentration(self, OCa, t):
        """ Calculates the TGF-beta concentration based on the osteoclastic resorption, external injection and
        degradation rate.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: TGF-beta concentration
        :rtype: float"""
        TGFb = (self.parameters.bone_volume.stored_TGFb_content * self.parameters.bone_volume.resorption_rate * OCa +
                self.calculate_external_injection_TGFb(t)) / self.parameters.degradation_rate.TGFb
        return TGFb

    def calculate_external_injection_TGFb(self, t):
        """ Calculates the external injection of TGF-beta (used and determined in load case scenarios).

        :param t: time variable
        :type t: float
        :return: concentration of external injection of TGF-beta
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.TGFb_injection

    def calculate_PTH_activation_OB(self, t):
        """ Calculates the activation of osteoblasts by parathyroid hormone (PTH) based on the current PTH concentration.

        :param t: time variable
        :type t: float
        :return: activation of osteoblasts by PTH
        :rtype: float"""
        PTH = self.calculate_PTH_concentration(t)
        PTH_activation_OB = PTH / (PTH + self.parameters.activation_coefficient.PTH_OB)
        return PTH_activation_OB

    def calculate_PTH_repression_OB(self, t):
        """ Calculates the repression of osteoblasts by parathyroid hormone (PTH) based on the current PTH concentration.

        :param t: time variable
        :type t: float
        :return: repression of osteoblasts by PTH
        :rtype: float"""
        PTH = self.calculate_PTH_concentration(t)
        PTH_repression_OB = self.parameters.repression_coefficient.PTH_OB / (
                    PTH + self.parameters.repression_coefficient.PTH_OB)
        return PTH_repression_OB

    def calculate_OPG_concentration(self, OBp, OBa, t):
        """ Calculates the osteoprotegerin (OPG) concentration based on the osteoblasts' production, PTH repression,
        maximum concentration, external injection and degradation rate. As described in the source publication,
        OPG can be produced by either precursor or active osteoblasts - determined by boolean variables.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: OPG concentration
        :rtype: float"""
        temp_PTH_OB = (self.parameters.production_rate.bool_OBp_produce_OPG *
                       self.parameters.production_rate.min_OPG_per_cell * OBp +
                       self.parameters.production_rate.bool_OBa_produce_OPG *
                       self.parameters.production_rate.min_OPG_per_cell * OBa) * self.calculate_PTH_repression_OB(t)
        OPG = (((temp_PTH_OB + self.calculate_external_injection_OPG(t)) * self.parameters.concentration.OPG_max) /
               (temp_PTH_OB + self.parameters.degradation_rate.OPG * self.parameters.concentration.OPG_max))
        return OPG

    def calculate_effective_carrying_capacity_RANKL(self, OBp, OBa, t):
        """ Calculates the effective carrying capacity of RANKL based on the osteoblasts' production, PTH activation and
        maximum RANKL per cell. As described in the source publication, RANKL can be produced by either precursor or
        active osteoblasts - determined by boolean variables.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: effective carrying capacity of RANKL
        :rtype: float"""
        RANKL_eff = (self.parameters.production_rate.bool_OBp_produce_RANKL *
                     self.parameters.production_rate.max_RANKL_per_cell * OBp +
                     self.parameters.production_rate.bool_OBa_produce_RANKL *
                     self.parameters.production_rate.max_RANKL_per_cell * OBa) * self.calculate_PTH_activation_OB(t)
        return RANKL_eff

    def calculate_RANKL_concentration(self, OBp, OBa, t):
        """ Calculates the RANKL concentration based on the effective carrying capacity, RANKL-RANK-OPG binding,
        degradation rate, intrinsic RANKL production and external injection of RANKL.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: RANKL concentration
        :rtype: float"""
        RANKL_eff = self.calculate_effective_carrying_capacity_RANKL(OBp, OBa, t)
        RANKL_RANK_OPG = RANKL_eff / (1 + self.parameters.binding_constant.RANKL_OPG *
                                      self.calculate_OPG_concentration(OBp, OBa, t) +
                                      self.parameters.binding_constant.RANKL_RANK * self.parameters.concentration.RANK)
        RANKL = RANKL_RANK_OPG * ((self.parameters.production_rate.intrinsic_RANKL +
                                   self.calculate_external_injection_RANKL(t)) /
                                  (self.parameters.production_rate.intrinsic_RANKL +
                                  self.parameters.degradation_rate.RANKL * RANKL_eff))
        return RANKL

    def calculate_RANKL_RANK_concentration(self, OBp, OBa, t):
        """ Calculates the RANKL-RANK complex concentration based on the RANKL concentration, RANK concentration and
        binding constant.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: RANKL-RANK complex concentration
        :rtype: float"""
        RANKL = self.calculate_RANKL_concentration(OBp, OBa, t)
        RANKL_RANK = self.parameters.binding_constant.RANKL_RANK * RANKL * self.parameters.concentration.RANK
        return RANKL_RANK

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        """ Calculates the activation of precursor osteoclasts by RANKL based on the RANKL-RANK complex concentration
        and activation coefficient.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: activation of precursor osteoclasts by RANKL
        :rtype: float"""
        RANKL_RANK = self.calculate_RANKL_RANK_concentration(OBp, OBa, t)
        RANKL_activation_OCp = RANKL_RANK / (RANKL_RANK + self.parameters.activation_coefficient.RANKL_RANK)
        return RANKL_activation_OCp

    # def calculate_external_injection_OCa(self, t):
    #     if t is None or t < self.load_case.start_time:
    #         return 0
    #     elif t > self.load_case.end_time:
    #         self.parameters.differentiation_rate.OCp = self.parameters.differentiation_rate.OCp / self.load_case.differentiation_rate_OCp_multiplier
    #         return 0
    #     elif self.load_case.start_time <= t <= self.load_case.end_time:
    #         self.is_load_case = True
    #         self.parameters.differentiation_rate.OCp = self.parameters.differentiation_rate.OCp * self.load_case.differentiation_rate_OCp_multiplier
    #         return self.load_case.OCa_injection
