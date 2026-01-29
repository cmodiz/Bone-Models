import numpy as np


class differentiation_rate:
    """ This class defines the differentiation rates of the different cell types.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +--------------+--------------------------+----------+
    |Parameter Name| Symbol                   | Units    |
    +==============+==========================+==========+
    | OBu          |:math:`D_{OB_u}^{Pivonka}`|  1/day   |
    +--------------+--------------------------+----------+
    | OBp          | :math:`D_{OB_p}`         |  1/day   |
    +--------------+--------------------------+----------+
    | OCp          | :math:`D_{OC_p}`         |  1/day   |
    +--------------+--------------------------+----------+

    :param OBu: differentiation rate of uncommitted osteoblasts
    :type OBu: float
    :param OBp: differentiation rate of precursor osteoblasts
    :type OBp: float
    :param OCp: differentiation rate of precursor osteoclasts
    :type OCp: float"""
    def __init__(self):
        # -> D_OB_u
        self.OBu = 0.7  # corrected differentiation rate of osteoblast progenitors [pM/day]
        # -> D_OB_p
        self.OBp = 0.165696312976030  # differentiation rate of preosteoblasts [pM/day]
        # self.OCu = 4.200000000000000e-003  # differentiation rate of uncommitted osteoclast [pM/day]
        # -> D_OC_p
        self.OCp = 2.1  # differentiation rate of preosteoclasts [pM/day]
        # -> D_OC_u
        self.OCu = 0.42   # differentiation rate of uncommitted osteoclast [1/day]


class apoptosis_rate:
    """ This class defines the apoptosis rates of the different cell types.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+------------------+---------+
    | Parameter Name   | Symbol           | Units   |
    +==================+==================+=========+
    | OBa              |:math:`A_{OB_a}`  | 1/day   |
    +------------------+------------------+---------+
    | OCa              | :math:`A_{OC_a}` |  1/day  |
    +------------------+------------------+---------+

    :param OBa: apoptosis rate of active osteoblasts
    :type OBa: float
    :param OCa: apoptosis rate of active osteoclasts
    :type OCa: float """
    def __init__(self):
        # -> A_OB_a
        self.OBa = 0.211072625806496  # apoptosis rate of active osteoblast [1/day]
        # -> A_OC_a
        self.OCa = 5.64874468409633   # apoptosis rate of active osteoclasts [pM/day]


class proliferation_rate:
    """ This class defines the proliferation rates. The prolifaration rate of OBp is depends the mechanics effect and is
    thus computed in the model (Eq. (16) in the paper).

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+------------------+---------+
    | Parameter Name   | Symbol           | Units   |
    +==================+==================+=========+
    | OBp              |:math:`P_{OB_p}`  | 1/day   |
    +------------------+------------------+---------+

    :param OBp: proliferation rate of precursor osteoblasts
    :type OBp: float
    """
    def __init__(self):
        # self.OBp_fraction = 0.1
        self.OBp = 3.5e-3  # proliferation rate of precursor osteoblasts [1/day]


class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand binding.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | TGFb_OBu         |:math:`K^{TGFb}_{act,OBu}` |  pM     |
    +------------------+---------------------------+---------+
    | TGFb_OCa         |:math:`K^{TGFb}_{act,OCa}` |  pM     |
    +------------------+---------------------------+---------+
    | PTH_OB           | :math:`K^{PTH}_{act,OB}`  |  pM     |
    +------------------+---------------------------+---------+
    | RANKL_RANK       |:math:`K_{d, [RANKL-RANK]}`|  pM     |
    +------------------+---------------------------+---------+

    :param TGFb_OBu: activation coefficient related to TGF-beta binding on OBu
    :type TGFb_OBu: float
    :param TGFb_OCa: activation coefficient related to TGF-beta binding on OCa
    :type TGFb_OCa: float
    :param PTH_OB: activation coefficient related to PTH binding to osteoblasts
    :type PTH_OB: float
    :param RANKL_RANK: activation coefficient related to RANKL binding on RANK
    :type RANKL_RANK: float"""
    def __init__(self):
        # Activation coefficients related to TGF-beta binding on OBu and OCa [pM]
        # -> K^{TGF-beta}_{act,OBu}
        self.TGFb_OBu = 0.000563278809675429
        # -> K^{TGF-beta}_{act,OCa}
        self.TGFb_OCa = 0.000563278809675429
        # Activation coefficient related to RANKL binding to RANK [pM]
        # -> K_{d, [RANKL-RANK]}
        # self.RANKL_OCp = 5.67971833061048
        # Activation coefficient for RANKL production related to PTH binding to osteoblasts [pM]
        # -> K^{PTH}_{act,OB}
        self.PTH_OB = 150
        # Activation coefficient related to MCSF binding on OCu [pM]
        # self.MCSF_OCu = None
        # Activation coefficient related to RANKL binding on RANK [pM]
        # -> K_{d, [RANKL-RANK]}
        self.RANKL_RANK = 16.65
        self.MCSF_OCu = 0.001


class repression_coefficient:
    """ This class defines the repression coefficients of respective receptor-ligand binding.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | TGFb_OBp         |:math:`K^{TGFb}_{rep,OBp}` |  pM     |
    +------------------+---------------------------+---------+
    | PTH_OB           | :math:`K^{PTH}_{rep,OB}`  |  pM     |
    +------------------+---------------------------+---------+

    :param TGFb_OBp: repression coefficient related to TGF-beta binding on OBp
    :type TGFb_OBp: float
    :param PTH_OB: repression coefficient for OPG production related to PTH binding on osteoblasts
    :type PTH_OB: float"""
    def __init__(self):
        # Repression coefficient related to TGF-beta binding on OBp [pM]
        # -> K_{rep, OBp}^{TGFb}
        self.TGFb_OBp = 0.000175426051821094
        # Repression coefficient for OPG production related to PTH binding on osteoblasts [pM]
        # -> K_{rep, OBp}^{PTH}
        self.PTH_OB = 0.222581427709954


class degradation_rate:
    r""" This class defines the degradation rates of the different factors.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+------------------------------+---------+
    | Parameter Name   | Symbol                       | Units   |
    +==================+==============================+=========+
    | PTH              |:math:`\check{D}_{PTH}`       |  1/day  |
    +------------------+------------------------------+---------+
    | OPG              |:math:`\check{D}_{OPG}`       |  1/day  |
    +------------------+------------------------------+---------+
    | RANKL            |:math:`\check{D}_{RANKL}`     |  1/day  |
    +------------------+------------------------------+---------+
    | TGFb             |:math:`\check{D}_{TGF-\beta}` |  1/day  |
    +------------------+------------------------------+---------+

    :param PTH: degradation rate of PTH
    :type PTH: float
    :param OPG: degradation rate of OPG
    :type OPG: float
    :param RANKL: degradation rate of RANKL
    :type RANKL: float
    :param TGFb: degradation rate of TGF-beta
    :type TGFb: float"""
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
        self.TGFb = 2


class concentration:
    """ This class defines fixed concentrations.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | OPG_max          |:math:`C^{max}_{OPG}`      |  pM     |
    +------------------+---------------------------+---------+
    | RANK             |-                          |  pM     |
    +------------------+---------------------------+---------+
    | OCp              |-                          |  pM     |
    +------------------+---------------------------+---------+"""
    def __init__(self):
        # -> C^max_OPG
        self.OPG_max = 2.00e+8  # Maximum concentration of OPG [pM]
        self.MCSF = 0.001
        # -> RANK
        self.RANK = 1.00e+1  # [pM] fixed concentration of RANK
        # self.OCp = 0.001  # [pM] fixed concentration of preosteoclasts


class binding_constant:
    """ This class defines the binding constants of RANK RANKL and OPG.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | RANKL_OPG        |:math:`K_{a, [RANKL-OPG]}` |1/pM     |
    +------------------+---------------------------+---------+
    | RANKL_RANK       |:math:`K_{a, [RANKL-RANK]}`|1/pM     |
    +------------------+---------------------------+---------+

    :param RANKL_OPG: association binding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: association binding constant for RANKL-RANK
    :type RANKL_RANK: float"""
    def __init__(self):
        # Association binding constant for RANKL-OPG [(pM)^{-1}]
        # -> K_{a, [RANKL-OPG]}
        self.RANKL_OPG = 1.00e-3
        # Association binding constant for RANKL-RANK [(pM)^{-1}]
        # -> K_{a, [RANKL-RANK]}
        self.RANKL_RANK = 3.411764705882353e-002
        # dissociation binding coefficient of TGFb with its receptor
        # [pM] value of OC to get half differentiation flux
        # -> C^s
        # self.TGFb_OC = 5.00e-3
        # [(pM day)^{-1}] rate of PTH binding with its receptor on OB
        # -> k_5
        # self.PTH_OB = 2.00e-2


class production_rate:
    r""" This class defines the intrinsic/ endogenous production rates of the different factors.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------+---------------------------+---------+
    | Parameter Name       | Symbol                    | Units   |
    +======================+===========================+=========+
    | intrinsic_PTH        |:math:`\beta_{PTH}`        |  pM/day |
    +----------------------+---------------------------+---------+
    | intrinsic_RANKL      |:math:`\beta_{RANKL}`      |  pM/day |
    +----------------------+---------------------------+---------+
    | min_OPG_per_cell     |:math:`p_{OB}^{OPG}`       |  pM     |
    +----------------------+---------------------------+---------+
    | bool_OBp_produce_OPG | -                         |  -      |
    +----------------------+---------------------------+---------+
    | bool_OBa_produce_OPG | -                         |  -      |
    +----------------------+---------------------------+---------+
    | max_RANKL_per_cell   |:math:`N_{RANKL}^{OB}`     |  pM     |
    +----------------------+---------------------------+---------+
    | max_RANK_per_cell    |:math:`C_{RANK}`           |  pM     |
    +----------------------+---------------------------+---------+
    |bool_OBp_produce_RANKL| -                         |  -      |
    +----------------------+---------------------------+---------+
    |bool_OBa_produce_RANKL| -                         |  -      |
    +----------------------+---------------------------+---------+

    :param intrinsic_PTH: intrinsic production rate of PTH
    :type intrinsic_PTH: float
    :param intrinsic_RANKL: intrinsic production rate of RANKL
    :type intrinsic_RANKL: float
    :param min_OPG_per_cell: minimal rate of OPG production per cell
    :type min_OPG_per_cell: float
    :param bool_OBp_produce_OPG: boolean variable determining which cells produce OPG
    :type bool_OBp_produce_OPG: int
    :param bool_OBa_produce_OPG: boolean variable determining which cells produce OPG
    :type bool_OBa_produce_OPG: int
    :param max_RANKL_per_cell: production rate of RANKL per cell
    :type max_RANKL_per_cell: float
    :param max_RANK_per_cell: production rate of RANK per cell
    :type max_RANK_per_cell: float
    :param bool_OBp_produce_RANKL: boolean variable determining which cells produce RANKL
    :type bool_OBp_produce_RANKL: int
    :param bool_OBa_produce_RANKL: boolean variable determining which cells produce RANKL
    :type bool_OBa_produce_RANKL: int
    """
    def __init__(self):
        # Intrinsic production rate of PTH [pM/day] (assumed to be constant)
        # -> beta_PTH
        self.intrinsic_PTH = 2.907
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
        # C_RANK
        self.max_RANK_per_cell = 1.000e+004
        # Boolean variables determining which cells produce RANKL
        self.bool_OBp_produce_RANKL = 1  # 1=yes
        self.bool_OBa_produce_RANKL = 0  # 0=no


class bone_volume:
    r""" This class defines the parameters relevant for bone volume of the bone model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------+-----------------------+---------+
    | Parameter Name       | Symbol                | Units   |
    +======================+=======================+=========+
    | formation_rate       |:math:`k_{form}`       |  1/day  |
    +----------------------+-----------------------+---------+
    | resorption_rate      |:math:`k_{res}`        |  1/day  |
    +----------------------+-----------------------+---------+
    |stored_TGFb_content   |:math:`\alpha`         |  -      |
    +----------------------+-----------------------+---------+
    |vascular_pore_fraction|:math:`f_{vas}`        |  %      |
    +----------------------+-----------------------+---------+
    | bone_fraction        |:math:`f_{bm}`         |  %      |
    +----------------------+-----------------------+---------+

    :param formation_rate: formation rate of bone volume
    :type formation_rate: float
    :param resorption_rate: resorption rate of bone volume
    :type resorption_rate: float
    :param stored_TGFb_content: proportionality constant expressing the TGF-β content stored in bone volume
    :type stored_TGFb_content: float
    :param vascular_pore_fraction: fraction of vascular pores in bone volume in percentage
    :type vascular_pore_fraction: float
    :param bone_fraction: fraction of bone matrix in bone volume in percentage
    :type bone_fraction: float"""
    def __init__(self):
        # -> k_form
        self.formation_rate = 9.011
        # -> k_res
        self.resorption_rate = 566.7
        # -> alpha
        self.stored_TGFb_content = 0.1  # proportionality constant expressing the TGF-β content stored in bone volume
        self.vascular_pore_fraction = 5  # fraction of vascular pores in bone volume in percentage
        self.bone_fraction = 95  # fraction of bone matrix in bone volume in percentage


class mechanics:
    r""" This class defines the parameters relevant for mechanics of the bone model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +-------------------------------------+--------------------------------------+---------+
    | Parameter Name                      | Symbol                               | Units   |
    +=====================================+======================================+=========+
    | strain_effect_on_OBp                |:math:`\check{\Pi}_{act, OB_p}^{mech}`|  -      |
    +-------------------------------------+--------------------------------------+---------+
    | strain_energy_density               |:math:`\check{\psi}_{bm}`             |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | update_OBp_proliferation_rate       | -                                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | OBu_differentiation_rate            |:math:`a_{P_{OB_p}}`                  |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | RANKL_production                    |:math:`P_{RANKL}`                     |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | bulk_modulus_water                  |:math:`k_{H_2O}`                      |  GPa    |
    +-------------------------------------+--------------------------------------+---------+
    | shear_modulus_water                 |:math:`\mu_{H_2O}`                    |  GPa    |
    +-------------------------------------+--------------------------------------+---------+
    | volumetric_part_of_unit_tensor      |:math:`\mathbb{J}`                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | unit_tensor_as_matrix               |:math:`\mathbb{I}`                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | deviatoric_part_of_unit_tensor      |:math:`\mathbb{K}`                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | stiffness_tensor_vascular_pores     |:math:`\mathbb{c}_{vas}`              |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | stiffness_tensor_bone_matrix        |:math:`\mathbb{c}_{bm}`               |  -      |
    +-------------------------------------+--------------------------------------+---------+
    |step_size_for_Hill_tensor_integration| -                                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | hill_tensor_cylindrical_inclusion   |:math:`\mathbb{P}_{r}^{bm}`           |        -|
    +-------------------------------------+--------------------------------------+---------+
    | stress_tensor_normal_loading        |:math:`\mathbb{\Sigma}_{cort}`        |  -      |
    +-------------------------------------+--------------------------------------+---------+

    :param strain_effect_on_OBp_steady_state: strain effect on OBp steady state
    :type strain_effect_on_OBp_steady_state: float
    :param strain_energy_density_steady_state: strain energy density steady state
    :type strain_energy_density_steady_state: float
    :param update_OBp_proliferation_rate: update OBp proliferation rate
    :type update_OBp_proliferation_rate: bool
    :param fraction_of_OBu_differentiation_rate: fraction of OBu differentiation rate
    :type fraction_of_OBu_differentiation_rate: float
    :param RANKL_production: RANKL production
    :type RANKL_production: float
    :param bulk_modulus_water: bulk modulus of water
    :type bulk_modulus_water: float
    :param shear_modulus_water: shear modulus of water
    :type shear_modulus_water: float
    :param volumetric_part_of_unit_tensor: volumetric part of unit tensor
    :type volumetric_part_of_unit_tensor: numpy.ndarray
    :param unit_tensor_as_matrix: unit tensor as matrix
    :type unit_tensor_as_matrix: numpy.ndarray
    :param deviatoric_part_of_unit_tensor: deviatoric part of unit tensor
    :type deviatoric_part_of_unit_tensor: numpy.ndarray
    :param stiffness_tensor_vascular_pores: stiffness tensor of vascular pores
    :type stiffness_tensor_vascular_pores: numpy.ndarray
    :param stiffness_tensor_bone_matrix: stiffness tensor of bone matrix
    :type stiffness_tensor_bone_matrix: numpy.ndarray
    :param step_size_for_Hill_tensor_integration: step size for Hill tensor integration
    :type step_size_for_Hill_tensor_integration: float
    :param hill_tensor_cylindrical_inclusion: Hill tensor of cylindrical inclusion
    :type hill_tensor_cylindrical_inclusion: numpy.ndarray
    :param stress_tensor_normal_loading: stress tensor of normal/ habitual loading
    :type stress_tensor_normal_loading: numpy.ndarray"""
    def __init__(self):
        # \breve{\Pi}_{act, OB_p}^{mech}
        self.strain_effect_on_OBp_steady_state = 0.5
        # \breve{\psi}_{bm}
        self.strain_energy_density_steady_state = None
        self.update_OBp_proliferation_rate = True
        # a_{P_{OB_p}}
        self.fraction_of_OBu_differentiation_rate = 0.1
        # P_{RANKL}
        self.RANKL_production = 0
        # k_{H_2O}
        self.bulk_modulus_water = 2.3  # [GPa]
        # \mu_{H_2O}
        self.shear_modulus_water = 0  # [GPa]
        # \mathbb{J}
        self.volumetric_part_of_unit_tensor = np.array([[1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0]]) * (1 / 3)
        # \mathbb{I}
        self.unit_tensor_as_matrix = np.array([[1, 0, 0, 0, 0, 0],
                                               [0, 1, 0, 0, 0, 0],
                                               [0, 0, 1, 0, 0, 0],
                                               [0, 0, 0, 1, 0, 0],
                                               [0, 0, 0, 0, 1, 0],
                                               [0, 0, 0, 0, 0, 1]])
        # \mathbb{K}
        self.deviatoric_part_of_unit_tensor = self.unit_tensor_as_matrix - self.volumetric_part_of_unit_tensor
        # \mathbb{c}_vas
        self.stiffness_tensor_vascular_pores = 3 * self.bulk_modulus_water * self.volumetric_part_of_unit_tensor + 2 * self.shear_modulus_water  * self.deviatoric_part_of_unit_tensor
        # \mathbb{c}_{bm}
        self.stiffness_tensor_bone_matrix = np.array([[18.5, 10.3, 10.4, 0, 0, 0],
                                                    [10.3, 20.8, 11.0, 0, 0, 0],
                                                    [10.4, 11.0, 28.4, 0, 0, 0],
                                                    [0, 0, 0, 12.9, 0, 0],
                                                    [0, 0, 0, 0, 11.5, 0],
                                                    [0, 0, 0, 0, 0, 9.3]])
        self.step_size_for_Hill_tensor_integration = 2 * np.pi / 50
        # \mathbb{P}_{r}^{bm}
        self.hill_tensor_cylindrical_inclusion = None
        # \Sigma_{cort}
        self.stress_tensor_normal_loading = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -30]]) * (10 ** -3)  # [GPa]


class calibration:
    def __init__(self):
        # delta
        self.turnover = 0.255
        # beta
        self.OCa = 0.09
        # gamma
        self.OBa = 1.132


class Lerebours_Parameters:
    """ This class defines the parameters of the bone model.

    :param differentiation_rate: differentiation rates of the different cell types
    :type differentiation_rate: differentiation_rate
    :param apoptosis_rate: apoptosis rates of the different cell types
    :type apoptosis_rate: apoptosis_rate
    :param activation_coefficient: activation coefficients of respective receptor-ligand binding
    :type activation_coefficient: activation_coefficient
    :param repression_coefficient: repression coefficients of respective receptor-ligand binding
    :type repression_coefficient: repression_coefficient
    :param degradation_rate: degradation rates of the different factors
    :type degradation_rate: degradation_rate
    :param concentration: fixed concentrations
    :type concentration: concentration
    :param binding_constant: binding constants of RANK RANKL and OPG
    :type binding_constant: binding_constant
    :param production_rate: intrinsic/ endogenous production rates of the different factors
    :type production_rate: production_rate
    :param mechanics: parameters relevant for mechanics of the bone model
    :type mechanics: mechanics
    :param bone_volume: parameters relevant for bone volume of the bone model
    :type bone_volume: bone_volume"""
    def __init__(self):
        self.differentiation_rate = differentiation_rate()
        self.apoptosis_rate = apoptosis_rate()
        self.activation_coefficient = activation_coefficient()
        self.repression_coefficient = repression_coefficient()
        self.degradation_rate = degradation_rate()
        self.concentration = concentration()
        self.binding_constant = binding_constant()
        self.production_rate = production_rate()
        self.mechanics = mechanics()
        self.proliferation_rate = proliferation_rate()
        self.bone_volume = bone_volume()
        self.calibration = calibration()
