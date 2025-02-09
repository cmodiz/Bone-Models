class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OBu             | D_OB_u       | pM/day        |
    | OBp             | D_OB_p      | pM/day         |
    | OCp             | D_OC_p      | pM/day         |
    +------------------+-----------+------------------+

     :param OBu: differentiation rate of uncommitted osteoblasts
     :type OBu: float
     :param OBp: differentiation rate of precursor osteoblasts
     :type OBp: float
     :param OCp: differentiation rate of precursor osteoclasts
     :type OCp: float
    """
    def __init__(self):
        # -> D_OB_u
        self.OBu = 7.00e-4  # corrected differentiation rate of osteoblast progenitors [pM/day]
        # -> D_OB_p
        # self.OBp = 2.674077909527713e-1   # differentiation rate of preosteoblasts [pM/day]
        self.OBp = 5.348   # differentiation rate of preosteoblasts [pM/day]
        self.OCu = None  # differentiation rate of uncommitted osteoclast [pM/day]
        # -> D_OC_p
        self.OCp = 2.1e-3  # differentiation rate of preosteoclasts [pM/day]


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OBa             | A_OB_a     | 1/day          |
    | OCa             | A_OC_a     | pM/day         |
    +------------------+-----------+------------------+

    :param OBa: apoptosis rate of active osteoblasts
    :type OBa: float
    :param OCa: apoptosis rate of active osteoclasts
    :type OCa: float """
    def __init__(self):
        # -> A_OB_a
        self.OBa = 1.89e-1  # apoptosis rate of active osteoblast [1/day]
        # -> A_OC_a
        self.OCa = 7.00e-1  # apoptosis rate of active osteoclasts [pM/day]


class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand bindings in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | TGFb_OBu         | K_{D1, TGF-beta} | pM           |
    | TGFb_OCa         | K_{D3, TGF-beta} | pM           |
    | PTH_OB           | K_{D4, PTH}, K_{D5, PTH} | pM |
    | RANKL_RANK       | K_{D8, RANKL} | pM           |
    +------------------+-----------+------------------+

    :param TGFb_OBu: activation coefficient related to TGF-beta binding on OBu
    :type TGFb_OBu: float
    :param TGFb_OCa: activation coefficient related to TGF-beta binding on OCa
    :type TGFb_OCa: float
    :param PTH_OB: activation coefficient related to PTH binding to osteoblasts
    :type PTH_OB: float
    :param RANKL_RANK: activation coefficient related to RANKL binding on RANK
    :type RANKL_RANK: float """
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
    """ This class defines the repression coefficients of respective receptor-ligand binding in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | TGFb_OBp         | K_{D2, TGF-beta} | pM           |
    | PTH_OB           | K_{D6, PTH}, K_{D7, PTH} | pM |
    +------------------+-----------+------------------+

    :param TGFb_OBp: repression coefficient related to TGF-beta binding on OBp
    :type TGFb_OBp: float
    :param PTH_OB: repression coefficient for OPG production related to PTH binding on osteoblasts
    :type PTH_OB: float """
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        # -> K_{D2, TGF-beta}
        self.TGFb_OBp = 1.415624253823446e-3
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        # -> K_{D6, PTH}, K_{D7, PTH}
        self.PTH_OB = 2.225814277099542e-1


class degradation_rate:
    """ This class defines the degradation rates of the different factors in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | PTH              | D^tilde_{PTH} | 1/day        |
    | OPG              | D^tilde_{OPG} | 1/day        |
    | RANKL            | D^tilde_{RANKL} | 1/day      |
    | TGFb             | D^tilde_{TGFb} | 1/day       |
    +------------------+-----------+------------------+

    :param PTH: degradation rate of PTH
    :type PTH: float
    :param OPG: degradation rate of OPG
    :type OPG: float
    :param RANKL: degradation rate of RANKL
    :type RANKL: float
    :param TGFb: degradation rate of TGFb
    :type TGFb: float """
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
    """ This class defines fixed concentrations in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | OPG_max          | C^max_OPG | pM              |
    | RANK             | C^RANK    | pM              |
    +------------------+-----------+------------------+

    :param OPG_max: maximum concentration of OPG
    :type OPG_max: float
    :param RANK: fixed concentration of RANK
    :type RANK: float """
    def __init__(self):
        # -> OPG_max
        self.OPG_max = 2.00e+8  # Maximum concentration of OPG [pM]
        self.MCSF = None
        # -> RANK
        self.RANK = 1.00e+1  # [pM] fixed concentration of RANK


class binding_constant:
    """ This class defines the binding constants of receptor-ligand bindings in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | RANKL_OPG        | K_A1,RANKL | (pM day)^{-1}  |
    | RANKL_RANK       | K_A2,RANKL | (pM day)^{-1}  |
    | TGFb_OC          | C^s        | pM              |
    | PTH_OB           | k_5        | (pM day)^{-1}  |
    +------------------+-----------+------------------+

    :param RANKL_OPG: association binding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: association binding constant for RANKL-RANK
    :type RANKL_RANK: float
    :param TGFb_OC: dissociation binding coefficient of TGFb with its receptor
    :type TGFb_OC: float
    :param PTH_OB: rate of PTH binding with its receptor on OB
    :type PTH_OB: float """
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
    """ This class defines the unbinding constants of receptor-ligand bindings in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | RANKL_OPG        | k_2       | 1/day           |
    | RANKL_RANK       | k_4       | 1/pM            |
    | PTH_OB           | k_6       | 1/day           |
    +------------------+-----------+------------------+

    :param RANKL_OPG: unbinding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: unbinding constant for RANKL-RANK
    :type RANKL_RANK: float
    :param PTH_OB: unbinding constant for PTH binding with its receptor on OB
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
    """ This class defines the intrinsic/ endogenous production rates of the different factors in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | intrinsic_PTH    | beta_PTH  | pM/day          |
    | intrinsic_RANKL  | beta_RANKL | pM/day          |
    | min_OPG_per_cell | beta_1,OPG, beta_2,OPG | pM/day |
    | bool_OBp_produce_OPG | -     | -              |
    | bool_OBa_produce_OPG | -     | -              |
    | max_RANKL_per_cell | R_1^RANKL, R_2^RANKL | pM/pM |
    | bool_OBp_produce_RANKL | -   | -              |
    | bool_OBa_produce_RANKL | -   | -              |
    +------------------+-----------+------------------+

    :param intrinsic_PTH: intrinsic production rate of PTH
    :type intrinsic_PTH: float
    :param intrinsic_RANKL: intrinsic production rate of RANKL
    :type intrinsic_RANKL: float
    :param min_OPG_per_cell: minimal rate of OPG production per cell
    :type min_OPG_per_cell: float
    :param bool_OBp_produce_OPG: boolean variable determining if OBp produce OPG
    :type bool_OBp_produce_OPG: int
    :param bool_OBa_produce_OPG: boolean variable determining if OBa produce OPG
    :type bool_OBa_produce_OPG: int
    :param max_RANKL_per_cell: production rate of RANKL per cell
    :type max_RANKL_per_cell: float
    :param bool_OBp_produce_RANKL: boolean variable determining if OBp produce RANKL
    :type bool_OBp_produce_RANKL: int
    :param bool_OBa_produce_RANKL: boolean variable determining if OBa produce RANKL
    :type bool_OBa_produce_RANKL: int """
    def __init__(self):
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        # -> beta_PTH
        self.intrinsic_PTH = 250
        # Intrinsic production rate of RANKL [pM/day]
        # -> beta_RANKL
        # Note: this value is e+4 in the paper but e+2 in the code
        self.intrinsic_RANKL = 1.684195714712206e+2
        # Minimal rate of OPG production per cell
        # -> beta_1,OPG, beta_2,OPG
        # different value in paper and code
        self.min_OPG_per_cell = 1.624900337835679e+008
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
    """ This class defines the correction factors in the Pivonka bone cell population model.
    Remark: f0 is implemented in the model (as in Lemaire 2004), but not explicitly given in the publication.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | f0               | f_0       | -               |
    +------------------+-----------+------------------+

    :param f0: correction factor for OBp differentiation rate and TGFb activation function
    :type f0: float """
    def __init__(self):
        # -> f_0
        self.f0 = 5.00e-2  # correction factor for OBp differentiation rate and TGFb activation function


class bone_volume:
    """ This class defines the parameters relevant for bone volume of the bone model in the Pivonka bone cell population model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units           |
    +==================+===========+==================+
    | formation_rate   | k_form    | 1/day           |
    | resorption_rate  | k_res     | 1/day           |
    | stored_TGFb_content | alpha  | -               |
    +------------------+-----------+------------------+

    :param formation_rate: formation rate of bone volume
    :type formation_rate: float
    :param resorption_rate: resorption rate of bone volume
    :type resorption_rate: float
    :param stored_TGFb_content: proportionality constant expressing the TGF-β content stored in bone volume
    :type stored_TGFb_content: float """
    def __init__(self):
        # -> k_form
        self.formation_rate = 1.571
        # -> k_res
        self.resorption_rate = 1
        # -> alpha
        self.stored_TGFb_content = 1.0  # proportionality constant expressing the TGF-β content stored in bone volume


class Pivonka_Parameters:
    """ This class defines the parameters of the Pivonka bone cell population model.

    :param differentiation_rate: differentiation rates of the different cell types, see :class:`differentiation_rate` for details
    :type differentiation_rate: differentiation_rate
    :param apoptosis_rate: apoptosis rates of the different cell types, see :class:`apoptosis_rate` for details
    :type apoptosis_rate: apoptosis_rate
    :param activation_coefficient: activation coefficients of respective receptor-ligand bindings, see :class:`activation
    _coefficient` for details
    :type activation_coefficient: activation_coefficient
    :param repression_coefficient: repression coefficients of respective receptor-ligand binding, see :class:`repression
    _coefficient` for details
    :type repression_coefficient: repression_coefficient
    :param correction_factor: correction factors, see :class:`correction_factor` for details
    :type correction_factor: correction_factor
    :param degradation_rate: degradation rates of the different factors, see :class:`degradation_rate` for details
    :type degradation_rate: degradation_rate
    :param concentration: fixed concentrations, see :class:`concentration` for details
    :type concentration: concentration
    :param binding_constant: binding constants of receptor-ligand bindings, see :class:`binding_constant` for details
    :type binding_constant: binding_constant
    :param unbinding_constant: unbinding constants of receptor-ligand bindings, see :class:`unbinding_constant` for details
    :type unbinding_constant: unbinding_constant
    :param production_rate: intrinsic/ endogenous production rates of the different factors, see :class:`production_rate` for details
    :type production_rate: production_rate
    :param bone_volume: parameters relevant for bone volume, see :class:`bone_volume` for details
    :type bone_volume: bone_volume """
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
        # this is correct but wrong in the paper (f0 is 0.05 in the Lemaire the paper)
        self.differentiation_rate.OBp = self.differentiation_rate.OBp * self.correction_factor.f0