class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OBu             | D_R       | pM/day          |
    | OBp             | d_B       | pM/day              |
    | OCp             | D_C       | pM/day              |
    +------------------+-----------+------------------+

     :param OBu: differentiation rate of uncommitted osteoblasts
     :type OBu: float
     :param OBp: differentiation rate of precursor osteoblasts
     :type OBp: float
     :param OCu: differentiation rate of uncommitted osteoclast
     :type OCu: float
     :param OCp: differentiation rate of precursor osteoclasts
     :type OCp: float
     """
    def __init__(self):
        """ Constructor method. """
        self.OBu = 7.00e-4
        self.OBp = 7.00e-1
        self.OCu = None
        self.OCp = 2.10e-3


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OBa             | k_B       | 1/day           |
    | OCa             | D_A       | pM/day          |
    +------------------+-----------+------------------+

    :param OBa: apoptosis rate of active osteoblasts
    :type OBa: float
    :param OCa: apoptosis rate of active osteoclasts
    :type OCa: float """
    def __init__(self):
        self.OBa = 1.89e-1
        self.OCa = 7.00e-1


class activation_coefficient:
    """ This class defines the activation coefficients of the different cell types in the Lemaire bone cell population model. """
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
    """ This class defines the repression coefficients of the different cell types in the Lemaire bone cell population model. """
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        self.TGFb_OBp = None
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        self.PTH_OB = None


class degradation_rate:
    """ This class defines the degradation rates of the different factors in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | PTH             | k_P       | 1/day           |
    | OPG             | k_O       | 1/day           |
    +------------------+-----------+------------------+

    :param PTH: degradation rate of PTH
    :type PTH: float
    :param OPG: degradation rate of OPG
    :type OPG: float
    :param RANKL: degradation rate of RANKL
    :type RANKL: float """
    def __init__(self):
        self.PTH = 86
        self.OPG = 3.50e-1
        self.RANKL = None


class concentration:
    """ This class defines the fixed concentrations of the different factors in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OPG_max         | OPG_max   | pM              |
    | RANK            | K        | pM              |
    +------------------+-----------+------------------+

    :param OPG_max: Maximum concentration of OPG
    :type OPG_max: float
    :param RANK: fixed concentration of RANK
    :type RANK: float """
    def __init__(self):
        """ Constructor method. """
        self.OPG_max = 2.00e+8  # Maximum concentration of OPG [pM]
        self.MCSF = None
        # -> K
        self.RANK = 1.00e+1  # [pM] fixed concentration of RANK


class binding_constant:
    """ This class defines the binding constants of the possible receptor-ligand interactions in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | RANKL_OPG       | k_1       | (pM day)^{-1}   |
    | RANKL_RANK      | k_3       | (pM day)^{-1}   |
    | TGFb_OC         | C^s       | pM              |
    | PTH_OB          | k_5       | (pM day)^{-1}   |
    +------------------+-----------+------------------+

    :param RANKL_OPG: Association binding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: Association binding constant for RANKL-RANK
    :type RANKL_RANK: float
    :param TGFb_OC: dissociation binding coefficient of TGFb with its receptor
    :type TGFb_OC: float
    :param PTH_OB: rate of PTH binding with its receptor on OB
    :type PTH_OB: float """
    def __init__(self):
        """ Constructor method. """
        # Association binding constant for RANKL-OPG [(pM day)^{-1}]
        # -> k_1
        self.RANKL_OPG = 1.00e-2
        # Association binding constant for RANKL-RANK [(pM day)^{-1}]
        # -> k_3
        self.RANKL_RANK = 5.80e-4
        # dissociation binding coefficient of TGFb with its receptor
        # [pM] value of OC to get half differentiation flux
        # -> C^s
        self.TGFb_OC = 5.00e-3
        # [(pM day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_5
        self.PTH_OB = 2.00e-2


class unbinding_constant:
    """ This class defines the unbinding constants of the possible receptor-ligand interactions in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | RANKL_OPG       | k_2       | day^{-1}        |
    | RANKL_RANK      | k_4       | day^{-1}        |
    | PTH_OB          | k_6       | day^{-1}        |
    +------------------+-----------+------------------+

    :param RANKL_OPG: Unbinding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: Unbinding constant for RANKL-RANK
    :type RANKL_RANK: float
    :param PTH_OB: rate of PTH unbinding with its receptor on OB
    :type PTH_OB: float """
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
    """ This class defines the intrinsic/ endogenous production rates of the different factors in the Lemaire bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | intrinsic_PTH    | S_P       | pM/day          |
    | intrinsic_RANKL  | r_L       | pM/day          |
    | min_OPG_per_cell | K^P_O     | pM/day          |
    | max_RANKL_per_cell| K^P_L     | pM/pM           |
    +------------------+-----------+------------------+

    :param intrinsic_PTH: Intrinsic production rate of PTH
    :type intrinsic_PTH: float
    :param intrinsic_RANKL: Intrinsic production rate of RANKL
    :type intrinsic_RANKL: float
    :param min_OPG_per_cell: Minimal rate of OPG production per cell
    :type min_OPG_per_cell: float
    :param max_RANKL_per_cell: Production rate of RANKL per cell
    :type max_RANKL_per_cell: float"""
    def __init__(self):
        """ Constructor method. """
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        # -> S_P
        self.intrinsic_PTH = 250
        # Intrinsic production rate of RANKL [pM/day]
        # -> r_L
        self.intrinsic_RANKL = 1.0e+3
        # Minimal rate of OPG production per cell
        # -> K^P_O
        self.min_OPG_per_cell = 2.00e+5
        # Boolean variables determining which cells produce OPG
        self.bool_OPGprod_OBp = None  # 0=no
        self.bool_OPGprod_OBa = None  # 1=yes
        # Constant describing how much RANK is produced per cell [pM/pM]
        self.RANKrate_per_cell = None
        # Production rate of RANKL per cell [pM/pM]
        # -> K^P_L
        self.max_RANKL_per_cell = 3.00e+6
        # Boolean variables determining which cells produce RANKL
        self.bool_RANKLprod_OBp = None  # 1=yes
        self.bool_RANKLprod_OBa = None  # 0=no


class correction_factor:
    """ This class defines the correction factors for the Lemaire bone cell population model. They are used to adjust the model parameters.

    :param f0: correction factor for OBp differentiation rate and TGFb activation function
    :type f0: float """
    def __init__(self):
        """ Constructor method. """
        # -> f_0
        self.f0 = 5.00e-2  # correction factor for OBp differentiation rate and TGFb activation function


class bone_volume:
    """ This class defines the parameters relevant for bone volume of the bone model.

    :param formation_rate: rate of bone formation
    :type formation_rate: float
    :param resorption_rate: rate of bone resorption
    :type resorption_rate: float """
    def __init__(self):
        self.formation_rate = 1.571/ 100
        self.resorption_rate = None


class Parameters:
    """ This class defines the parameters of the Lemaire bone cell population model.

    :param differentiation_rate: differentiation rates of the different cell types, see :class:`differentiation_rate` for details
    :type differentiation_rate: differentiation_rate
    :param apoptosis_rate: apoptosis rates of the different cell types, see :class:`apoptosis_rate` for details
    :type apoptosis_rate: apoptosis_rate
    :param activation_coefficient: activation coefficients of the different cell types, see :class:`activation
    _coefficient` for details
    :type activation_coefficient: activation_coefficient
    :param repression_coefficient: repression coefficients of the different cell types, see :class:`repression
    _coefficient` for details
    :type repression_coefficient: repression_coefficient
    :param correction_factor: correction factors for the model, see :class:`correction
    _factor` for details
    :type correction_factor: correction_factor
    :param degradation_rate: degradation rates of the different factors, see :class:`degradation_rate` for details
    :type degradation_rate: degradation_rate
    :param concentration: fixed concentrations of the different factors, see :class:`concentration` for details
    :type concentration: concentration
    :param binding_constant: binding constants of the receptor-ligand interactions, see :class:`binding_constant` for details
    :type binding_constant: binding_constant
    :param unbinding_constant: unbinding constants of the receptor-ligand interactions, see :class:`unbinding_constant` for details
    :type unbinding_constant: unbinding_constant
    :param production_rate: intrinsic production rates of the different factors, see :class:`production
    _rate` for details
    :type production_rate: production_rate
    :param bone_volume: parameters relevant for bone volume, see :class:`bone_volume` for details
    :type bone_volume: bone_volume
    :param differentiation_rate.OBp: corrected differentiation rate of precursor osteoblasts
    :type differentiation_rate.OBp: float """
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

        self.differentiation_rate.OBp = self.differentiation_rate.OBp * self.correction_factor.f0
