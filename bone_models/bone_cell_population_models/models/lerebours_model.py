import numpy as np
from scipy.integrate import solve_ivp
import scipy as sc
from ..parameters.lerebours_parameters import Lerebours_Parameters
from .scheiner_model import Scheiner_Model


class Lerebours_Model(Scheiner_Model):
    """
    Implements the Lerebours mechanobiological model for bone cell population dynamics.

    This model extends the Scheiner framework by incorporating porosity-dependent
    feedback loops. It couples bone cell populations (osteoblasts/osteoclasts)
    to the bone volume fraction via the specific surface area, allowing for
    simulations of site-specific bone loss and remodeling.

    .. note::
       **Source Publication**:
       Lerebours, C., Buenzli, P. R., Scheiner, S., & Pivonka, P. (2016).
       *A multiscale mechanobiological model of bone remodeling predicts
       site-specific bone loss in the femur during osteoporosis and mechanical disuse.*
       Biomechanics and Modeling in Mechanobiology, 15(1), 43-67.
       :doi:`10.1007/s10237-015-0705-x`
    """
    def __init__(self, load_case, porosity, specific_surface_multiplier=1):
        """
        Initializes the model with local geometry and loading conditions.

        Args:
            load_case (Load_Case): Object containing mechanical loading parameters.
            porosity (float): Initial bone porosity (0 to 1), used to derive
                the specific surface area for cell recruitment.
            specific_surface_multiplier (float, optional): Adjusts the surface
                availability for cell attachment. Defaults to 1.
        """
        super().__init__(load_case)
        self.parameters = Lerebours_Parameters()
        self.initial_guess_root = np.array([0.0001, 0.0001, 0.001, 0.0001, porosity, 1 - porosity])
        self.steady_state = type('', (), {})()
        self.steady_state.OBu = None
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCu = None
        self.steady_state.OCp = None
        self.steady_state.OCa = None
        self.specific_surface_multiplier = specific_surface_multiplier
        self.update_mechanical_effects = True
        self.mechanical_effect_on_OBp = 0

    def bone_cell_population_model(self, x, t=None):
        """
        Defines the ODE system for bone cell populations and volume fractions.

        This function supports two modes:
        1. **Steady State** (t=None): Calculates residuals for uncommitted and precursor cells given fixed active cell densities.
        2. **Transient** (t=float): Calculates rates of change for precursors, active cells, and volume fractions.

        :param x: State vector. If t is None: [OBu, OBp, OCu, OCp, vascular_pore_fraction, bone_volume_fraction]. If t is float: [OBp, OBa, OCp, OCa, vascular_pore_fraction, bone_volume_fraction].
        :type x: list
        :param t: Time variable; if None, the model assumes a steady-state calculation.
        :type t: float or None
        :return: List of rates of change [dOBp/dt, dOBa/dt, dOCp/dt, dOCa/dt, dVPF/dt, dBVF/dt].
        :rtype: list
        """
        if t is None:
            # steady state for given OBa, OCa
            OCa = self.steady_state.OCa
            OBa = self.steady_state.OBa
            OBu, OBp, OCu, OCp, vascular_pore_fraction, bone_volume_fraction = x
        else:
            OCu = self.steady_state.OCu
            OBu = self.steady_state.OBu
            OBp, OBa, OCp, OCa, vascular_pore_fraction, bone_volume_fraction = x

        dOBpdt = (self.parameters.differentiation_rate.OBu * self.calculate_TGFb_activation_OBu(OCa, t) * OBu -
                  self.parameters.differentiation_rate.OBp * OBp *
                  self.calculate_TGFb_repression_OBp(OCa, t) + self.apply_mechanical_effects(OBp, OBa, OCa,
                                                                                             vascular_pore_fraction * 100,
                                                                                             bone_volume_fraction * 100,
                                                                                             t))

        dOBadt = (self.parameters.differentiation_rate.OBp * OBp * self.calculate_TGFb_repression_OBp(OCa, t) -
                  self.parameters.apoptosis_rate.OBa * OBa)

        dOCpdt = (self.parameters.differentiation_rate.OCu * self.calculate_RANKL_activation_OCu(OBp, OBa, t) *
                  self.calculate_MCSF_activation_OCu() * OCu - self.parameters.differentiation_rate.OCp *
                  self.calculate_RANKL_activation_OCp(OBp, OBa, t) * OCp)

        dOCadt = (self.parameters.differentiation_rate.OCp * self.calculate_RANKL_activation_OCp(OBp, OBa, t) * OCp -
                  self.parameters.apoptosis_rate.OCa * OCa * self.calculate_TGFb_activation_OCa(OCa, t))

        dvascular_pore_fractiondt = self.parameters.bone_volume.resorption_rate * OCa - self.parameters.bone_volume.formation_rate * OBa
        dbone_volume_fractiondt = self.parameters.bone_volume.formation_rate * OBa - self.parameters.bone_volume.resorption_rate * OCa
        dxdt = [dOBpdt, dOBadt, dOCpdt, dOCadt, dvascular_pore_fractiondt, dbone_volume_fractiondt]
        return dxdt

    def calculate_RANKL_activation_OCu(self, OBp, OBa, t):
        """
        Calculates the RANKL-mediated activation of uncommitted osteoclasts (OCu).
        This implementation assumes the activation rate for uncommitted cells is
        identical to that of osteoclast precursors (OCp).

        :param OBp: Concentration of precursor osteoblasts.
        :type OBp: float
        :param OBa: Concentration of active osteoblasts.
        :type OBa: float
        :param t: Time variable.
        :type t: float
        :return: RANKL activation factor for OCu.
        :rtype: float
        """
        RANKL_activation_OCu = self.calculate_RANKL_activation_OCp(OBp, OBa, t)
        return RANKL_activation_OCu

    def calculate_RANKL_activation_OCp(self, OBp, OBa, t):
        """ Calculates the RANKL activation for osteoclast precursor cells (OCp) based on the RANKL concentration.
        It calls the function calculate_RANKL_concentration.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: RANKL activation for osteoclast precursor cells (OCp)
        :rtype: float"""
        RANKL = self.calculate_RANKL_concentration(OBp, OBa, t)
        RANKL_activation_OCp = RANKL / (RANKL + self.parameters.activation_coefficient.RANKL_RANK)
        return RANKL_activation_OCp

    def calculate_RANKL_concentration(self, OBp, OBa, t):
        """ Calculates the RANKL concentration based on the effective carrying capacity, RANKL-RANK-OPG binding,
        degradation rate, intrinsic RANKL production and external injection of RANKL. An additional RANKL production is
        added due to mechanical effects.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: RANKL concentration
        :rtype: float """
        RANKL_eff = self.calculate_effective_carrying_capacity_RANKL(OBp, OBa, t)
        RANKL_RANK_OPG = RANKL_eff / (1 + self.parameters.binding_constant.RANKL_OPG *
                                      self.calculate_OPG_concentration(OBp, OBa, t) +
                                      self.parameters.binding_constant.RANKL_RANK * self.parameters.concentration.RANK)
        RANKL = RANKL_RANK_OPG * ((self.parameters.production_rate.intrinsic_RANKL * OBp +
                                   self.calculate_external_injection_RANKL(
                                       t) + self.parameters.mechanics.RANKL_production) /
                                  (self.parameters.production_rate.intrinsic_RANKL * OBp +
                                   self.parameters.degradation_rate.RANKL * RANKL_eff))
        return RANKL

    def calculate_PTH_concentration(self, t):
        """ Calculates the PTH concentration based on intrinsic PTH production and external injection of PTH.

        :param t: time variable
        :type t: float
        :return: PTH concentration
        :rtype: float"""
        PTH = self.parameters.production_rate.intrinsic_PTH + self.calculate_external_injection_PTH(t)
        return PTH

    def calculate_MCSF_activation_OCu(self):
        """ Calculates the MCSF activation for uncommitted osteoclast cells (OCu) based on the MCSF concentration
        (constant parameter) and activation coefficient (Hill function).

        :return: MCSF activation for uncommitted osteoclast cells (OCu)
        :rtype: float"""
        MCSF_activation_OCu = self.parameters.concentration.MCSF / (
                self.parameters.concentration.MCSF + self.parameters.activation_coefficient.MCSF_OCu)
        return MCSF_activation_OCu

    def calculate_TGFb_concentration(self, OCa, t):
        """ Calculates the TGFb concentration based on the stored TGFb content in bone volume, active osteoclasts,
        resorption rate, degradation rate and calibration parameter.

        :param OCa: active osteoclast cell concentration
        :type OCa: float
        :param t: time variable
        :type t: float
        :return: TGFb concentration
        :rtype: float """
        TGFb_concentration = ((self.parameters.bone_volume.stored_TGFb_content * OCa *
                               self.parameters.bone_volume.resorption_rate * (
                                      1 / self.parameters.calibration.OCa)) /
                              self.parameters.degradation_rate.TGFb)
        return TGFb_concentration

    def calculate_OPG_concentration(self, OBp, OBa, t):
        """ Calculates the OPG concentration based on the effective carrying capacity of OPG, OPG production,
        calibration parameter, PTH repression and external injection of OPG.

        :param OBp: precursor osteoblast cell concentration
        :type OBp: float
        :param OBa: active osteoblast cell concentration
        :type OBa: float
        :param t: time variable
        :type t: float
        :return: OPG concentration
        :rtype: float"""
        temp_PTH_OB = ((self.parameters.production_rate.bool_OBp_produce_OPG *
                        self.parameters.production_rate.min_OPG_per_cell * OBp +
                        self.parameters.production_rate.bool_OBa_produce_OPG *
                        self.parameters.production_rate.min_OPG_per_cell * OBa) * (
                               1 / self.parameters.calibration.OBa) *
                       self.calculate_PTH_repression_OB(t))
        OPG = (((temp_PTH_OB + self.calculate_external_injection_OPG(t)) * self.parameters.concentration.OPG_max) /
               (temp_PTH_OB + self.parameters.degradation_rate.OPG * self.parameters.concentration.OPG_max))
        return OPG

    def specific_surface(self, porosity):
        """ This function calculates the specific surface of bone based on the porosity. The specific surface multiplier
        is used for the boundaries in Modiz et al. (irrelevant for this model), default is 1.

        :param porosity: Porosity of the bone
        :type porosity: float
        :return: Specific surface of bone
        :rtype: float"""
        specific_surface = self.specific_surface_multiplier * (
                32.2 * (1 - porosity) - 93.9 * (1 - porosity) ** 2 + 134 * (1 - porosity) ** 3 - 101 * (
                1 - porosity) ** 4 + 28.8 * (1 - porosity) ** 5)
        return specific_surface

    def calculate_turnover(self, porosity):
        """ Calculates the turnover based on porosity, calibration factor and specific surface.
        The calibration factor was identified by fitting the SV curve to turnover datapoints.

        :param porosity: Porosity of the bone
        :type porosity: float
        :return: Turnover rate
        :rtype: float """
        # calibration_factor = 3.996636532576335e-05
        # calibration_factor = 0.255
        if porosity == 0 or porosity == 1:
            turnover = 0
        else:
            turnover = self.parameters.calibration.turnover * self.specific_surface(porosity)
        return turnover

    def calculate_steady_state(self, porosity):
        """ Calculates the steady state for the bone cell population model based on the porosity. It determines the
        active osteoclasts and osteoblasts based on the turnover, and then solves the bone cell population model for
        steady state of uncommitted and precursor cells. The steady state concentrations are stored as parameters
        (uncommitted stays constant for all further calculation).

        :param porosity: Porosity of the bone
        :type porosity: float
        :return: None
        """
        turnover = self.calculate_turnover(porosity)
        self.steady_state.OCa = turnover / self.parameters.bone_volume.resorption_rate
        self.steady_state.OBa = turnover / self.parameters.bone_volume.formation_rate

        cells_steady_state = sc.optimize.root(self.bone_cell_population_model, self.initial_guess_root, tol=1e-15,
                                              options={'xtol': 1e-15}, method='lm')
        self.steady_state.OBu = cells_steady_state.x[0]
        self.steady_state.OBp = cells_steady_state.x[1]
        self.steady_state.OCu = cells_steady_state.x[2]
        self.steady_state.OCp = cells_steady_state.x[3]
        print(f"Steady state calculated for porosity {porosity:.2f}: OBu={self.steady_state.OBu:}, "
              f"OBp={self.steady_state.OBp:}, OBa={self.steady_state.OBa:}, "
              f"OCu={self.steady_state.OCu:}, OCp={self.steady_state.OCp:}, "
              f"OCa={self.steady_state.OCa:}, Turnover in % per day ={turnover:}")
        pass

    def solve_bone_cell_population_model(self, tspan, porosity, initial_conditions=None):
        """ Solve the bone cell population model and volume fractions using the ODE system over a given time interval.
        The initial conditions are set to the steady-state values which are calculated if initial conditions are None.
        This function is overwritten from the source model to add vascular pore volume fraction and bone volume fraction,
        that are necessary to solve in every time step.

        :param tspan: time span for the ODE solver
        :type tspan: numpy.ndarray with start and end time
        :param porosity: Porosity of the bone ([0,1]), used to calculate specific surface and steady state
        :type porosity: float
        :param initial_conditions: Initial conditions for the ODE system, if None, steady state is calculated
        :type initial_conditions: list or None
        :return: solution of the ODE system
        :rtype: scipy.integrate._ivp.ivp.OdeResult
        """
        if initial_conditions is None:
            self.calculate_steady_state(porosity)
            x0 = np.array(
                [self.steady_state.OBp, self.steady_state.OBa, self.steady_state.OCp, self.steady_state.OCa, porosity,
                 1 - porosity])
        else:
            x0 = initial_conditions
        solution = solve_ivp(lambda t, x: self.bone_cell_population_model(x, t), tspan, x0, rtol=1e-8, atol=1e-10,
                             method='BDF', max_step=1)
        if not solution.success:
            print(f"Integration failed: {solution.message}")
        return solution

    def apply_mechanical_effects(self, OBp, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t):
        """ Computes the mechanical stimulus-driven change in osteoblast precursor (OBp) proliferation.

        This method applies a piecewise mechanotransduction function to the baseline
        proliferation rate. The effect is determined by the strain energy density and
        the biomechanical transduction strength defined in the model parameters.

        :param OBp: Precursor osteoblast concentration.
        :type OBp: float
        :param OBa: Active osteoblast concentration.
        :type OBa: float
        :param OCa: Active osteoclast concentration.
        :type OCa: float
        :param vascular_pore_fraction: Vascular pore volume fraction.
        :type vascular_pore_fraction: float
        :param bone_volume_fraction: Bone volume fraction.
        :type bone_volume_fraction: float
        :param t: Time variable.
        :type t: float
        :return: Rate of change of OBp concentration due to mechanical effects.
        :rtype: float
        """
        strain_effect_on_OBp = self.calculate_strain_effect_on_OBp(OBa, OCa, vascular_pore_fraction,
                                                                      bone_volume_fraction, t)
        self.strain_effect_on_OBp = strain_effect_on_OBp
        if strain_effect_on_OBp <= 0:
            return self.parameters.proliferation_rate.OBp * OBp
        elif 0 < strain_effect_on_OBp < 1 / self.parameters.mechanics.biomech_transduction_strength:
            return self.parameters.proliferation_rate.OBp * (
                    1 + self.parameters.mechanics.biomech_transduction_strength * strain_effect_on_OBp) * OBp
        else:
            return 2 * self.parameters.proliferation_rate.OBp * OBp

    def calculate_strain_effect_on_OBp(self, OBa, OCa, vascular_pore_fraction, bone_volume_fraction, t):
        """
        Calculates the normalized mechanical stimulus based on strain energy density (SED).
        In steady-state (t=None), this method initializes the baseline SED. During loading
        scenarios, it computes the relative deviation from this baseline to determine
        the mechanical effect on OBp proliferation and adapt RANKL production rates.

        :param OBa: Active osteoblast concentration.
        :type OBa: float
        :param OCa: Active osteoclast concentration.
        :type OCa: float
        :param vascular_pore_fraction: Vascular pore volume fraction.
        :type vascular_pore_fraction: float
        :param bone_volume_fraction: Bone volume fraction.
        :type bone_volume_fraction: float
        :param t: Time variable; if None or before start_time, returns 0 (steady-state).
        :type t: float or None
        :return: Normalized strain effect (relative change in SED).
        :rtype: float
        """
        if t is None:
            if self.parameters.mechanics.strain_energy_density_steady_state is None:
                self.parameters.mechanics.strain_energy_density_steady_state = self.calculate_strain_energy_density(OBa,
                                                                                                                    OCa,
                                                                                                                    vascular_pore_fraction,
                                                                                                                    bone_volume_fraction,
                                                                                                                    t=0)
            return 0  # no effect in steady-state
        elif t <= self.load_case.start_time:
            return 0  # no effect in steady-state
        else:
            strain_energy_density = self.calculate_strain_energy_density(OBa, OCa, vascular_pore_fraction,
                                                                         bone_volume_fraction, t)
            self.strain_energy_density = strain_energy_density
            strain_effect_on_OBp = (
                        (strain_energy_density - self.parameters.mechanics.strain_energy_density_steady_state) /
                        (self.parameters.mechanics.strain_energy_density_steady_state + self.parameters.mechanics.correction_factor))
            if strain_effect_on_OBp > 0:
                self.parameters.mechanics.RANKL_production = 0
            elif strain_effect_on_OBp <= 0:
                self.parameters.mechanics.RANKL_production = - self.parameters.mechanics.biomech_transduction_strength_RANKL * strain_effect_on_OBp
            return strain_effect_on_OBp

    def set_macroscopic_stress_tensor(self, stress_xx, stress_xy, stress_xz, steady_state=False):
        """
        Sets the macroscopic stress tensor and updates the model or load case parameters.

        .. note::
           **Coordinate System Convention**: To maintain compatibility between the Scheiner
           and Lerebours frameworks, the input ``stress_xx`` is mapped to the internal
           z-direction (tensor index [2,2]). This discrepancy is essential for the
           correct operation of the spatial model component.

        :param stress_xx: Axial stress (mapped to the internal z-direction).
        :type stress_xx: float
        :param stress_xy: Shear stress component in the xy plane.
        :type stress_xy: float
        :param stress_xz: Shear stress component in the xz plane.
        :type stress_xz: float
        :param steady_state: If True, updates the normal loading parameters; otherwise, updates the current load case.
        :type steady_state: bool
        :return: None
        """
        stress_tensor = np.zeros((3, 3))
        stress_tensor[2, 2] = stress_xx
        stress_tensor[0, 2] = stress_xz
        stress_tensor[2, 0] = stress_xz
        stress_tensor[1, 2] = stress_xy
        stress_tensor[2, 1] = stress_xy
        if steady_state:
            self.parameters.mechanics.stress_tensor_normal_loading = stress_tensor
        else:
            self.load_case.stress_tensor = stress_tensor
        pass
