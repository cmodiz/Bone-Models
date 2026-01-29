import numpy as np
from scipy.integrate import solve_ivp
import scipy as sc
from bone_models.bone_cell_population_models.parameters.lerebours_parameters import Lerebours_Parameters
from bone_models.bone_cell_population_models.models.scheiner_model import Scheiner_Model


class Lerebours_Model(Scheiner_Model):
    def __init__(self, load_case, porosity):
        super().__init__(load_case)
        self.parameters = Lerebours_Parameters()
        self.initial_guess_root = np.array([0.001, 0.001, 0.001, 0.001, porosity,  1-porosity])
        self.steady_state = type('', (), {})()
        self.steady_state.OBu = None
        self.steady_state.OBp = None
        self.steady_state.OBa = None
        self.steady_state.OCu = None
        self.steady_state.OCp = None
        self.steady_state.OCa = None

    def bone_cell_population_model(self, x, t=None):
        if t is None:
            OBu, OBp, OCu, OCp, vascular_pore_fraction, bone_volume_fraction = x
            OCa = self.steady_state.OCa
            OBa = self.steady_state.OBa
        else:
            OBp, OBa, OCp, OCa, vascular_pore_fraction, bone_volume_fraction = x
            OCu = self.steady_state.OCu
            OBu = self.steady_state.OBu

        dOBpdt = (self.parameters.differentiation_rate.OBu * self.calculate_TGFb_activation_OBu(OCa, t) * OBu -
                  self.parameters.differentiation_rate.OBp * OBp *
                  self.calculate_TGFb_repression_OBp(OCa, t) + self.apply_mechanical_effects(OBp, OBa, OCa,
                                                                                             vascular_pore_fraction,
                                                                                             bone_volume_fraction, t))

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
        RANKL_activation_OCu = self.calculate_RANKL_activation_OCp(OBp, OBa, t)
        return RANKL_activation_OCu

    def calculate_MCSF_activation_OCu(self):
        MCSF_activation_OCu = self.parameters.concentration.MCSF / (
                    self.parameters.concentration.MCSF + self.parameters.activation_coefficient.MCSF_OCu)
        return MCSF_activation_OCu

    def calculate_TGFb_concentration(self, OCa, t):
        TGFb_concentration = (self.parameters.bone_volume.stored_TGFb_content * OCa *
                              self.parameters.bone_volume.resorption_rate * (1/self.parameters.calibration.OCa) /
                              self.parameters.degradation_rate.TGFb)
        return TGFb_concentration

    def calculate_OPG_concentration(self, OBp, OBa, t):
        temp_PTH_OB = (self.parameters.production_rate.bool_OBp_produce_OPG *
                       self.parameters.production_rate.min_OPG_per_cell * OBp +
                       self.parameters.production_rate.bool_OBa_produce_OPG *
                       self.parameters.production_rate.min_OPG_per_cell * OBa) * self.calculate_PTH_repression_OB(t)
        OPG = (((temp_PTH_OB * self.parameters.calibration.OBa + self.calculate_external_injection_OPG(t)) * self.parameters.concentration.OPG_max) /
               (temp_PTH_OB + self.parameters.degradation_rate.OPG * self.parameters.concentration.OPG_max))
        return OPG

    def specific_surface(self, porosity):
        """ This function calculates the specific surface of bone based on the porosity. """
        specific_surface = 32.2 * (1 - porosity) - 93.9 * (1 - porosity) ** 2 + 134 * (1 - porosity) ** 3 - 101 * (
                    1 - porosity) ** 4 + 28.8 * (1 - porosity) ** 5
        return specific_surface

    def calculate_turnover(self, porosity):
        # calibration_factor = 3.996636532576335e-05
        # calibration_factor = 0.255
        turnover = self.parameters.calibration.turnover * self.specific_surface(porosity)
        return turnover

    def calculate_steady_state(self, porosity):
        activation_frequency = self.calculate_turnover(porosity)
        self.steady_state.OBa = activation_frequency / self.parameters.bone_volume.resorption_rate
        self.steady_state.OCa = activation_frequency / self.parameters.bone_volume.formation_rate

        cells_steady_state = sc.optimize.root(self.bone_cell_population_model, self.initial_guess_root, tol=1e-15,
                                              options={'xtol': 1e-15}, method='lm')
        self.steady_state.OBu = cells_steady_state.x[0]
        self.steady_state.OBp = cells_steady_state.x[1]
        self.steady_state.OCu = cells_steady_state.x[2]
        self.steady_state.OCp = cells_steady_state.x[3]
        pass
