import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from ..parameters.lemaire_parameters import Lemaire_Parameters
from .base_model import Base_Model


class Lemaire_Model(Base_Model):
    """ This class implements the bone cell population model by Lemaire et al. (2004) as a subclass of the Base_Model class.

    .. note::
       **Source Publication**:
       Lemaire, V., Tobin, F. L., Greller, L. D., Cho, C. R., & Suva, L. J. (2004).
       *Modeling the interactions between osteoblast and osteoclast activities in bone remodeling.*
       Journal of Theoretical Biology, 229(3), 293-309.
       :doi:`10.1016/j.jtbi.2004.03.023`

    :param load_case: load case for the model
    :type load_case: object
    :param parameters: model parameters
    :type parameters: Lemaire_Parameters
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
        """ Constructor for the Lemaire_Model class.

        :param load_case: load case for the model
        :type load_case: object """
        super().__init__()
        self.parameters = Lemaire_Parameters()
        self.initial_guess_root = np.array([0.7734e-3, 0.7282e-3, 0.9127e-3])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None
        self.load_case = load_case

    def calculate_TGFb_activation_OBu(self, OCa, t):
        """ Calculate the activation of uncommitted osteoblasts by TGF-beta.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float

        :return: activation of uncommitted osteoblasts by TGF-beta
        :rtype: float"""
        TGFb_activation_OBu = ((OCa + self.parameters.correction_factor.f0 * self.parameters.binding_constant.TGFb_OC) /
                               (OCa + self.parameters.binding_constant.TGFb_OC))
        return TGFb_activation_OBu

    def calculate_TGFb_repression_OBp(self, OCa, t):
        """ Calculate the repression of precursor osteoblasts by TGF-beta.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float

        :return: repression of precursor osteoblasts by TGF-beta
        :rtype: float"""

        TGFb_repression_OBp = 1 / self.calculate_TGFb_activation_OBu(OCa, t)
        return TGFb_repression_OBp

    def calculate_TGFb_activation_OCa(self, OCa, t):
        """ Calculate the activation of active osteoclasts by TGF-beta.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float

        :return: activation of active osteoclasts by TGF-beta
        :rtype: float"""
        TGFb_activation_OCp = self.calculate_TGFb_activation_OBu(OCa, t)
        return TGFb_activation_OCp

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        """ Calculate the activation of precursor osteoclasts by RANKL.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float

        :return: activation of precursor osteoclasts by RANKL
        :rtype: float"""
        PTH_effect = self.parameters.production_rate.max_RANKL_per_cell * self.calculate_PTH_activation_OB(t) * OBa
        OPG_concentration = self.calculate_OPG_concentration(OBp, t)
        kinetics_RANKL_RANK = self.parameters.binding_constant.RANKL_RANK/self.parameters.unbinding_constant.RANKL_RANK
        kinetics_RANKL_OPG = self.parameters.binding_constant.RANKL_OPG/self.parameters.unbinding_constant.RANKL_OPG
        temp = kinetics_RANKL_RANK * (PTH_effect / (
                1 + kinetics_RANKL_RANK * self.parameters.concentration.RANK + kinetics_RANKL_OPG * OPG_concentration))
        RANKL_activation_OCp = temp * (1 + self.calculate_external_injection_RANKL(t) / self.parameters.production_rate.intrinsic_RANKL)
        return RANKL_activation_OCp

    def calculate_PTH_activation_OB(self, t):
        """ Calculate the activation of osteoblasts by PTH.

        :param t: time variable
        :type t: float

        :return: activation of osteoblasts by PTH
        :rtype: float"""
        PTH = self.calculate_PTH_concentration(t)
        PTH_kinetic = self.parameters.unbinding_constant.PTH_OB / self.parameters.binding_constant.PTH_OB
        PTH_activation_OB = PTH / (self.calculate_external_injection_PTH(t) / self.parameters.degradation_rate.PTH + PTH_kinetic)
        return PTH_activation_OB

    def calculate_PTH_concentration(self, t):
        """ Calculate the PTH concentration depending on intrinsic and external PTH and degradation rate.

        :param t: time variable
        :type t: float

        :return: PTH concentration
        :rtype: float"""
        PTH = ((self.parameters.production_rate.intrinsic_PTH + self.calculate_external_injection_PTH(t)) /
               self.parameters.degradation_rate.PTH)
        return PTH

    def calculate_OPG_concentration(self, OBp, t):
        """ Calculate the OPG concentration depending on intrinsic and external OPG and degradation rate.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param t: time variable
        :type t: float

        :return: OPG concentration
        :rtype: float"""
        OPG = (1 / self.parameters.degradation_rate.OPG) * (self.parameters.production_rate.min_OPG_per_cell * OBp /
                                                    self.calculate_PTH_activation_OB(t)
                                                    + self.calculate_external_injection_OPG(t))
        return OPG

    def calculate_external_injection_OBp(self, t):
        """ Calculate the external injection of precursor osteoblasts (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of precursor osteoblasts
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OBp_injection

    def calculate_external_injection_OBa(self, t):
        """ Calculate the external injection of active osteoblasts (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of active osteoblasts
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OBa_injection

    def calculate_external_injection_OCa(self, t):
        """ Calculate the external injection of active osteoclasts (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of active osteoclasts
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OCa_injection

    def calculate_external_injection_PTH(self, t):
        """ Calculate the external injection of PTH (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of PTH
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.PTH_injection

    def calculate_external_injection_OPG(self, t):
        """ Calculate the external injection of OPG (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of OPG
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.OPG_injection

    def calculate_external_injection_RANKL(self, t):
        """ Calculate the external injection of RANKL (used in load case scenarios).

        :param t: time variable
        :type t: float

        :return: external injection of RANKL
        :rtype: float"""
        if (t is None) or t < self.load_case.start_time or t > self.load_case.end_time:
            return 0
        else:
            return self.load_case.RANKL_injection

    def calculate_bone_volume_fraction_change(self, time, solution, steady_state, initial_bone_volume_fraction):
        """ Calculate the bone volume fraction over time numerically using the solution of the ODE system and the trapezoidal rule.

        :param time: time variable
        :type time: numpy.ndarray
        :param solution: solution of the ODE system
        :type solution: numpy.ndarray
        :param steady_state: steady state values of the model
        :type steady_state: numpy.ndarray
        :param initial_bone_volume_fraction: initial bone volume fraction
        :type initial_bone_volume_fraction: float

        :return: bone volume fraction over time
        :rtype: list"""
        self.parameters.bone_volume.resorption_rate = self.parameters.bone_volume.formation_rate * steady_state[1]/steady_state[2]
        bone_volume_fraction = [initial_bone_volume_fraction]
        dBV_dt = self.parameters.bone_volume.formation_rate * (solution[1][:]) - self.parameters.bone_volume.resorption_rate * (solution[2][:])  # Compute dBV/dt
        for i in range(1, len(time)):
            dt = time[i] - time[i - 1]
            trapezoid = (dBV_dt[i - 1] + dBV_dt[i]) / 2 * dt
            bone_volume_fraction.append(bone_volume_fraction[-1] + trapezoid)
        return bone_volume_fraction
