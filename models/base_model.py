import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from parameters.base_parameters import Parameters


class Base_Model:
    """This class implements the base model for bone cell population models. It essentially includes the ODE system that
    is shared across future models. It is not intended to be used as a model but only indicates the basis of the models.

    :param OBp: precursor osteoblast cell concentration
    :type OBp: float
    :param OBa: active osteoblast cell concentration
    :type OBa: float
    :param OCa: osteoclast cell concentration
    :type OCa: float
    :param parameters: Model parameters, see :class:`Parameters` for details.
    :type parameters: Parameters
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
    :type bone_volume_fraction: list
    """

    def __init__(self):
        """ Constructor for the Base_Model class. """
        self.parameters = Parameters()
        self.initial_guess_root = np.array([None, None, None])
        self.steady_state = type('', (), {})()
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCa = None

    def bone_cell_population_model(self, x, t=None):
        """ ODE system for the bone cell population model.

        :param x: state variables of the model
        :type x: list
        :param t: time variable
        :type t: float
        :return: rate of change of state variables
        :rtype: list
        """
        OBp, OBa, OCa = x
        dOBpdt = ((self.parameters.differentiation_rate.OBu * self.calculate_TGFb_activation_OBu(OCa,t) -
                   self.parameters.differentiation_rate.OBp * OBp *
                   self.calculate_TGFb_repression_OBp(OCa,t) + self.calculate_external_injection_OBp(t)) +
                  self.apply_mechanical_effects() + self.apply_medication_effects_OBp())
        dOBadt = (self.parameters.differentiation_rate.OBp * OBp * self.calculate_TGFb_repression_OBp(OCa,t) -
                  self.parameters.apoptosis_rate.OBa * OBa * self.apply_medication_effects_OBa()
                  + self.calculate_external_injection_OBa(t))
        dOCadt = (self.parameters.differentiation_rate.OCp * self.calculate_RANKL_activation_OCp(OBp,OBa,t) -
                  self.parameters.apoptosis_rate.OCa * OCa * self.calculate_TGFb_activation_OCa(OCa,t) +
                  self.calculate_external_injection_OCa(t))
        dxdt = [dOBpdt, dOBadt, dOCadt]
        return dxdt

    def calculate_steady_state(self):
        """ Calculate the steady state of the bone cell population model using root finding of the ODE system.
        :return: steady state values of the model
        :rtype: numpy.ndarray
        """
        print('Calculating steady state ...', end='')
        steady_state = root(self.bone_cell_population_model, self.initial_guess_root, tol=1e-30, method="lm", options={'xtol': 1e-30}) #tol=1e-5)
        self.steady_state.OBp = steady_state.x[0]
        self.steady_state.OBa = steady_state.x[1]
        self.steady_state.OCa = steady_state.x[2]
        print(f'done \n Steady state: {steady_state.x}')
        return steady_state.x

    def solve_bone_cell_population_model(self, tspan):
        """ Solve the bone cell population model using the ODE system.

        :param tspan: time span for the ODE solver
        :type tspan: numpy.ndarray with start and end time
        :return: solution of the ODE system
        :rtype: scipy.integrate._ivp.ivp.OdeResult
        """
        x0 = self.calculate_steady_state()
        print('Solving bone cell population model ...', end='')
        solution = solve_ivp(lambda t, x: self.bone_cell_population_model(x, t), tspan, x0, rtol=1e-8, atol=1e-8)
        print('done')
        return solution

    def calculate_TGFb_activation_OBu(self, OCa, t):
        """ Calculate the activation of uncommitted osteoblasts by TGFb.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: activation of uncommitted osteoblasts by TGFb
        :rtype: float
        """
        pass

    def calculate_TGFb_repression_OBp(self, OCa, t):
        """ Calculate the repression of precursor osteoblasts by TGFb.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: repression of precursor osteoblasts by TGFb"""
        pass

    def calculate_TGFb_activation_OCa(self, OCa, t):
        """ Calculate the activation of active osteoclasts by TGFb.

        :param OCa: osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: activation of active osteoclasts by TGFb
        :rtype: float"""
        pass

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
        pass

    def calculate_external_injection_OBp(self, t):
        """ Calculate the external injection of precursor osteoblasts (used in load case scenarios).

        :param t: time variable
        :type t: float
        :return: external injection of precursor osteoblasts
        :rtype: float"""
        pass

    def calculate_external_injection_OBa(self, t):
        """ Calculate the external injection of active osteoblasts (used in load case scenarios).

        :param t: time variable
        :type t: float
        :return: external injection of active osteoblasts
        :rtype: float"""
        pass

    def calculate_external_injection_OCa(self, t):
        """ Calculate the external injection of active osteoclasts (used in load case scenarios).

        :param t: time variable
        :type t: float
        :return: external injection of active osteoclasts
        :rtype: float"""
        pass

    def calculate_bone_volume_fraction_change(self, solution, steady_state, initial_bone_volume_fraction):
        """ Calculate the change in bone volume fraction over time.

        :param solution: solution of the ODE system
        :type solution: scipy.integrate._ivp.ivp.OdeResult
        :param steady_state: steady state values of the model
        :type steady_state: numpy.ndarray
        :param initial_bone_volume_fraction: initial bone volume fraction
        :type initial_bone_volume_fraction: float
        :return: bone volume fraction over time
        :rtype: list
        """
        bone_volume_fraction = [initial_bone_volume_fraction]
        for i in range(len(solution)):
            bone_volume_fraction.append(
                bone_volume_fraction[-1] +
                self.parameters.bone_volume.formation_rate * (solution[i][1] - steady_state[1]) -
                self.parameters.bone_volume.resorption_rate * (solution[i][2] - steady_state[2])
            )
        return bone_volume_fraction

    def apply_mechanical_effects(self):
        """ Apply mechanical effects to the bone cell population model.

        :return: mechanical effects acting on the bone cell population model
        :rtype: float
        """
        return 0

    def apply_medication_effects_OBp(self):
        """ Apply medication effects to precursor osteoblasts.

        :return: medication effects acting on precursor osteoblasts
        :rtype: float
        """
        return 0

    def apply_medication_effects_OBa(self):
        """ Apply medication effects to active osteoblasts.

        :return: medication effects acting on active osteoblasts
        :rtype: float
        """
        return 1

