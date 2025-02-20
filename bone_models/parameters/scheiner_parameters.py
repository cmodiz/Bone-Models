import numpy as np

class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types """
    def __init__(self):
        # -> D_OB_u
        self.OBu = 7.00e-4  # corrected differentiation rate of osteoblast progenitors [pM/day]
        # -> D_OB_p
        self.OBp = 0.165696312976030  # differentiation rate of preosteoblasts [pM/day]
        self.OCu = 4.200000000000000e-003  # differentiation rate of uncommitted osteoclast [pM/day]
        # -> D_OC_p
        self.OCp = 2.100000000000000e+000 * 0.001 # differentiation rate of preosteoclasts [pM/day]


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types. """
    def __init__(self):
        # -> A_OB_a
        self.OBa = 0.211072625806496  # apoptosis rate of active osteoblast [1/day]
        # -> A_OC_a
        self.OCa = 5.64874468409633   # apoptosis rate of active osteoclasts [pM/day]


class proliferation_rate:
    def __init__(self):
        # self.OBp_fraction = 0.1
        self.OBp = 0


class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand binding. """
    def __init__(self):
        # Activation coefficients related to TGF-beta binding on OBu and OCa [pM]
        # -> K^{TGF-beta}_{act,OBu}
        self.TGFb_OBu = 0.000563278809675429
        # -> K^{TGF-beta}_{act,OCa}
        self.TGFb_OCa = 0.000563278809675429
        # Activation coefficient related to RANKL binding to RANK [pM]
        # -> K_{d, RANKL-RANK}
        self.RANKL_OCp = 5.67971833061048
        # Activation coefficient for RANKL production related to PTH binding to osteoblasts [pM]
        # -> K^{PTH}_{act,OB}
        self.PTH_OB = 150
        # Activation coefficient related to MCSF binding on OCu [pM]
        self.MCSF_OCu = None
        # Activation coefficient related to RANKL binding on RANK [pM]
        # -> K_{d, [RANKL-RANK]}
        self.RANKL_RANK = 5.6797


class repression_coefficient:
    """ This class defines the repression coefficients of respective receptor-ligand binding. """
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        # -> K_{rep, OBp}^{TGFb}
        self.TGFb_OBp = 0.000175426051821094
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        # -> K_{rep, OBp}^{PTH}
        self.PTH_OB = 0.222581427709954


class degradation_rate:
    """ This class defines the degradation rates of the different factors. """
    def __init__(self):
        # Degradation rate of PTH [1/day]
        # -> D^tilde_{PTH}
        self.PTH = 86
        # Degradation rate of OPG [1/day]
        # -> D^tilde_{OPG}
        self.OPG = 3.50e-1
        # Degradation rate of RANKL [1/day]
        # -> D^tilde_{RANKL}
        self.RANKL = 1.0132471014805027e+1
        # Degradation rate of TGFb [1/day]
        # -> D^tilde_{TGFb}
        self.TGFb = 1.00e+0


class concentration:
    """ This class defines fixed concentrations. """
    def __init__(self):
        # -> C^max_OPG
        self.OPG_max = 2.00e+8  # Maximum concentration of OPG [pM]
        # self.MCSF = None
        # -> RANK
        self.RANK = 1.00e+1  # [pM] fixed concentration of RANK
        self.OCp = 0.001  # [pM] fixed concentration of preosteoclasts


class binding_constant:
    """ This class defines the binding constants of RANK RANKL and OPG. """
    def __init__(self):
        # Association binding constant for RANKL-OPG [(pM day)^{-1}]
        # -> K_{a, RANKL-OPG}
        self.RANKL_OPG = 1.00e-3
        # Association binding constant for RANKL-RANK [(pM day)^{-1}]
        # -> K_{a, RANKL-RANK}
        self.RANKL_RANK = 3.411764705882353e-002
        # dissociation binding coefficient of TGFb with its receptor
        # [pM] value of OC to get half differentiation flux
        # -> C^s
        # self.TGFb_OC = 5.00e-3
        # [(pM day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_5
        # self.PTH_OB = 2.00e-2


# class unbinding_constant:
#     """ This class defines the unbinding constants of RANK RANKL and OPG. """
#     def __init__(self):
        # Association binding constant for RANKL-OPG [1/day]
        # -> k_2
        # self.RANKL_OPG = 1.00e+1
        # Association binding constant for RANKL-RANK [1/pM]
        # -> k_4
        # self.RANKL_RANK = 1.70e-2
        # dissociation binding coefficient of TGFb with its receptor
        # self.TGFb_OC = None
        # [(day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_6
        # self.PTH_OB = 3.00e+0


class production_rate:
    """ This class defines the intrinsic/ endogenous production rates of the different factors."""
    def __init__(self):
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        # -> beta_PTH
        self.intrinsic_PTH = 250
        # Intrinsic production rate of RANKL [pM/day]
        # -> beta_RANKL
        # Note: this value is e+4 in the paper but e+2 in the code
        self.intrinsic_RANKL = 1.684195714712206e+2
        # Minimal rate of OPG production per cell
        # -> p_OB^{OPG}
        self.min_OPG_per_cell = 1.624900337835679e+008
        # Boolean variables determining which cells produce OPG
        self.bool_OBp_produce_OPG = 0  # 0=no
        self.bool_OBa_produce_OPG = 1  # 1=yes
        # Constant describing how much RANKL is produced per cell [pM/pM]
        # self.RANKL_rate_per_cell = 2.703476379131062e+006
        # Production rate of RANKL per cell [pM/pM]
        # -> N_{RANKL}^OB
        self.max_RANKL_per_cell = 2.703476379131062e+006
        # Production rate of RANK per cell [pM/pM]
        self.max_RANK_per_cell = 1.000e+004
        # Boolean variables determining which cells produce RANKL
        self.bool_OBp_produce_RANKL = 1  # 1=yes
        self.bool_OBa_produce_RANKL = 0  # 0=no


class correction_factor:
    """ This class defines the correction factors. """
    def __init__(self):
        # -> f_0
        self.f0 = 5.00e-2  # correction factor for OBp differentiation rate and TGFb activation function


class bone_volume:
    """ This class defines the parameters relevant for bone volume of the bone model. """
    def __init__(self):
        # -> k_form
        self.formation_rate = 40
        # -> k_res
        self.resorption_rate = 200
        # -> alpha
        self.stored_TGFb_content = 1.0 # proportionality constant expressing the TGF-β content stored in bone volume
        self.vascular_pore_fraction = 5  # fraction of vascular pores in bone volume in percentage
        self.bone_fraction = 95  # fraction of bone matrix in bone volume in percentage


class mechanics:
    def __init__(self):
        self.strain_effect_on_OBp_steady_state = 0.5
        self.strain_effect_on_OBp = np.array(self.strain_effect_on_OBp_steady_state)
        self.strain_energy_density_steady_state = None
        self.update_OBp_proliferation_rate = True
        self.fraction_of_OBu_differentiation_rate = 0.1
        self.RANKL_production = 0
        # self.OBu_proliferation_rate = 0
        self.bulk_modulus_water = 2.3  # [GPa]
        self.shear_modulus_water = 0  # [GPa]
        self.volumetric_part_of_unit_tensor = np.array([[1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0]]) * (1 / 3)
        self.unit_tensor_as_matrix = np.array([[1, 0, 0, 0, 0, 0],
                                               [0, 1, 0, 0, 0, 0],
                                               [0, 0, 1, 0, 0, 0],
                                               [0, 0, 0, 1, 0, 0],
                                               [0, 0, 0, 0, 1, 0],
                                               [0, 0, 0, 0, 0, 1]])
        self.deviatoric_part_of_unit_tensor = self.unit_tensor_as_matrix - self.volumetric_part_of_unit_tensor
        self.stiffness_tensor_vascular_pores = 3 * self.bulk_modulus_water * self.volumetric_part_of_unit_tensor + 2 * self.shear_modulus_water  * self.deviatoric_part_of_unit_tensor
        self.stiffness_tensor_bone_matrix = np.array([[18.5, 10.3, 10.4, 0, 0, 0],
                                                    [10.3, 20.8, 11.0, 0, 0, 0],
                                                    [10.4, 11.0, 28.4, 0, 0, 0],
                                                    [0, 0, 0, 12.9, 0, 0],
                                                    [0, 0, 0, 0, 11.5, 0],
                                                    [0, 0, 0, 0, 0, 9.3]])
        self.step_size_for_Hill_tensor_integration = 2 * np.pi / 50
        self.hill_tensor_cylindrical_inclusion = None
        self.stress_tensor_normal_loading = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -30]]) * (10 ** -3)  # [GPa]


class Scheiner_Parameters:
    """ This class defines the parameters of the bone model. """
    def __init__(self):
        self.differentiation_rate = differentiation_rate()
        self.apoptosis_rate = apoptosis_rate()
        self.activation_coefficient = activation_coefficient()
        self.repression_coefficient = repression_coefficient()
        self.correction_factor = correction_factor()
        self.degradation_rate = degradation_rate()
        self.concentration = concentration()
        self.binding_constant = binding_constant()
        # self.unbinding_constant = unbinding_constant()
        self.production_rate = production_rate()
        self.mechanics = mechanics()
        self.proliferation_rate = proliferation_rate()
        # self.capacity = capacity()
        self.bone_volume = bone_volume()
        # self.differentiation_rate.OBp = self.differentiation_rate.OBp * self.correction_factor.f0