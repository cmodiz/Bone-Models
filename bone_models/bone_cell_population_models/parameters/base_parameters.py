class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types. """
    def __init__(self):
        self.OBu = None  # differentiation rate of osteoblast progenitors [pM/day]
        self.OBp = None  # differentiation rate of preosteoblasts [pM/day]
        self.OCu = None  # differentiation rate of uncommitted osteoclast [pM/day]
        self.OCp = None  # differentiation rate of preosteoclasts [pM/day]


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types. """
    def __init__(self):
        self.OBa = None  # apoptosis rate of active osteoblast [1/day]
        self.OCa = None  # apoptosis rate of active osteoclasts [pM/day]


class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand binding. """
    def __int__(self):
        # Activation coefficients related to TGF-beta binding on OBu and OCa [pM]
        self.TGFb_OBu = None
        self.TGFb_OCa = None
        # Activation coefficient related to RANKL binding to RANK [pM]
        self.RANKL_OCp = None
        # Activation coefficient for RANKL production related to PTH binding to osteoblasts [pM]
        self.PTH_OB = None
        # Activation coefficient related to MCSF binding on OCu [pM]
        self.MCSF_OCu = None


class repression_coefficient:
    """ This class defines the repression coefficients of respective receptor-ligand binding. """
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        self.TGFb_OBp = None
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        self.PTH_OB = None


class external_dosage:
    """ This class defines the external dosages of the different cell types.
    They are 0 but can be altered in a load case scenario to simulate external injections."""
    def __init__(self):
        # External RANKL-injection [pM/day]
        self.RANKL = None
        # External OPG injection [pM/day]
        self.OPG = None
        # External PTH injection [pM/day]
        self.PTH = None
        # External OBa injection [pM/day]
        self.OBa = None
        # External OBp injection [pM/day]
        self.OBp = None
        # External OCY injection [pM/day]
        self.OCY = None
        # External OCp injection [pM/day]
        self.OCp = None
        # External OCa injection [pM/day]
        self.OCa = None


class degradation_rate:
    """ This class defines the degradation rates of the different factors. """
    def __init__(self):
        # Degradation rate of PTH [1/day]
        self.PTH = None
        # Degradation rate of OPG [1/day]
        self.OPG = None
        # Degradation rate of RANKL [1/day]
        self.RANKL = None


class concentration:
    """ This class defines fixed concentrations. """
    def __init__(self):
        self.OPG_max = None  # Maximum concentration of OPG [pM]
        self.MCSF = None


class binding_constant:
    """ This class defines the binding constants of RANK RANKL and OPG. """
    def __init__(self):
        # Association binding constant for RANKL-OPG [1/pM]
        self.RANKL_OPG = None
        # Association binding constant for RANKL-RANK [1/pM]
        self.RANKL_RANK = None


class production_rate:
    """ This class defines the intrinsic/ endogenous production rates of the different factors."""
    def __init__(self):
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        self.intrinsic_PTH = None
        # Intrinsic production rate of RANKL [pM/day]
        self.intrinsic_RANKL = None
        # Rate of OPG production per cell
        self.OPGrate_per_cell = None
        # Boolean variables determining which cells produce OPG
        self.bool_OPGprod_OBp = None  # 0=no
        self.bool_OPGprod_OBa = None  # 1=yes
        # Constant describing how much RANK is produced per cell [pM/pM]
        self.RANKrate_per_cell = None
        # Production rate of RANKL per cell [pM/pM]
        self.RANKLrate_per_cell = None
        # Boolean variables determining which cells produce RANKL
        self.bool_RANKLprod_OBp = None  # 1=yes
        self.bool_RANKLprod_OBa = None  # 0=no


class correction_factor:
    """ This class defines the correction factors. """
    def __init__(self):
        # -> f_0
        self.f0 = None  # correction factor for OBp differentiation rate and TGFb activation function


class bone_volume:
    """ This class defines the parameters relevant for bone volume of the bone model. """
    def __init__(self):
        self.formation_rate = None
        self.resorption_rate = None


class Base_Parameters:
    """ This class defines the parameters of the bone model. """
    def __init__(self):
        self.differentiation_rate = differentiation_rate()
        self.apoptosis_rate = apoptosis_rate()
        self.activation_coefficient = activation_coefficient()
        self.repression_coefficient = repression_coefficient()
        self.external_dosage = external_dosage()
        self.degradation_rate = degradation_rate()
        self.concentration = concentration()
        self.binding_constant = binding_constant()
        self.production_rate = production_rate()
        self.correction_factor = correction_factor()
        self.bone_volume = bone_volume()