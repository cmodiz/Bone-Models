class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types """
    def __init__(self):
        # -> D_OB_u
        self.OBu = 7.00e-4  # corrected differentiation rate of osteoblast progenitors [pM/day]
        # -> D_OB_p
        # self.OBp = 2.674077909527713e-1   # differentiation rate of preosteoblasts [pM/day]
        self.OBp = 5.348   # differentiation rate of preosteoblasts [pM/day]
        self.OCu = None  # differentiation rate of uncommitted osteoclast [pM/day]
        # -> D_OC_p
        self.OCp = 2.10e-3  # differentiation rate of preosteoclasts [pM/day]


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types. """
    def __init__(self):
        # -> A_OB_a
        self.OBa = 1.89e-1  # apoptosis rate of active osteoblast [1/day]
        # -> A_OC_a
        self.OCa = 7.00e-1  # apoptosis rate of active osteoclasts [pM/day]


class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand binding. """
    def __init__(self):
        # Activation coefficients related to TGF-beta binding on OBu and OCa [pM]
        # -> K_{D1, TGF-beta}
        self.TGFb_OBu = 4.545454545454545e-3
        # -> K_{D3, TGF-beta}
        self.TGFb_OCa = 4.545454545454545e-3
        # Activation coefficient related to RANKL binding to RANK [pM]
        # -> K_{D6, PTH}, K_{D7, PTH}
        # self.RANKL_OCp = 4.457452802710724
        # Activation coefficient for RANKL production related to PTH binding to osteoblasts [pM]
        # -> K_{D4, PTH}, K_{D5, PTH}
        self.PTH_OB = 1.5e+2
        # Activation coefficient related to MCSF binding on OCu [pM]
        self.MCSF_OCu = None
        # Activation coefficient related to RANKL binding on RANK [pM]
        # -> K_{D8, RANKL}
        self.RANKL_RANK = 1.306e+1


class repression_coefficient:
    """ This class defines the repression coefficients of respective receptor-ligand binding. """
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        # -> K_{D2, TGF-beta}
        self.TGFb_OBp = 1.415624253823446e-3
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        # -> K_{D6, PTH}, K_{D7, PTH}
        self.PTH_OB = 2.225814277099542e-1


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
        # -> OPG_max
        self.OPG_max = 2.00e+8  # Maximum concentration of OPG [pM]
        self.MCSF = None
        # -> RANK
        self.RANK = 1.00e+1  # [pM] fixed concentration of RANK


class binding_constant:
    """ This class defines the binding constants of RANK RANKL and OPG. """
    def __init__(self):
        # Association binding constant for RANKL-OPG [(pM day)^{-1}]
        # -> K_A1,RANKL
        self.RANKL_OPG = 1.00e-3
        # Association binding constant for RANKL-RANK [(pM day)^{-1}]
        # -> K_A2,RANKL
        self.RANKL_RANK = 3.411764705882353e-2
        # dissociation binding coefficient of TGFb with its receptor
        # [pM] value of OC to get half differentiation flux
        # -> C^s
        self.TGFb_OC = 5.00e-3
        # [(pM day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_5
        self.PTH_OB = 2.00e-2


class unbinding_constant:
    """ This class defines the unbinding constants of RANK RANKL and OPG. """
    def __init__(self):
        # Association binding constant for RANKL-OPG [1/day]
        # -> k_2
        self.RANKL_OPG = 1.00e+1
        # Association binding constant for RANKL-RANK [1/pM]
        # -> k_4
        self.RANKL_RANK = 1.70e-2
        # dissociation binding coefficient of TGFb with its receptor
        self.TGFb_OC = None
        # [(day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_6
        self.PTH_OB = 3.00e+0


class production_rate:
    """ This class defines the intrinsic/ endogenous production rates of the different factors."""
    def __init__(self):
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        # -> beta_PTH
        self.intrinsic_PTH = 250
        # Intrinsic production rate of RANKL [pM/day]
        # -> beta_RANKL
        self.intrinsic_RANKL = 1.684195714712206e+004
        # Minimal rate of OPG production per cell
        # -> beta_1,OPG, beta_2,OPG
        self.min_OPG_per_cell = 1.464e+8
        # Boolean variables determining which cells produce OPG
        self.bool_OBp_produce_OPG = 0  # 0=no
        self.bool_OBa_produce_OPG = 1  # 1=yes
        # Constant describing how much RANKL is produced per cell [pM/pM]
        self.RANKL_rate_per_cell = None
        # Production rate of RANKL per cell [pM/pM]
        # -> R_1^RANKL, R_2^RANKL
        self.max_RANKL_per_cell = 3e+006
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
        self.formation_rate = 1.571
        # -> k_res
        self.resorption_rate = 200.0
        # -> alpha
        self.stored_TGFb_content = 1.0  # proportionality constant expressing the TGF-β content stored in bone volume


class Parameters:
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
        self.unbinding_constant = unbinding_constant()
        self.production_rate = production_rate()
        # self.capacity = capacity()
        self.bone_volume = bone_volume()