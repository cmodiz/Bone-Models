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
    | OCu          | :math:`D_{OC_u}`         |  1/day   |
    +--------------+--------------------------+----------+

    :param OBu: differentiation rate of uncommitted osteoblasts
    :type OBu: float
    :param OBp: differentiation rate of precursor osteoblasts
    :type OBp: float
    :param OCp: differentiation rate of precursor osteoclasts
    :type OCp: float
    :param OCu: differentiation rate of uncommitted osteoclasts
    :type OCu: float """

    def __init__(self):
        self.OBu = 0.7
        self.OBp = 0.165696312976030
        self.OCp = 2.1
        self.OCu = 0.42


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
        self.OBa = 0.211072625806496
        self.OCa = 5.64874468409633


class proliferation_rate:
    """ This class defines the proliferation rates. The proliferation rate of OBp depends the mechanics effect and is
    thus computed in the model (Eq. (22) in the paper).

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
        self.OBp = 3.5e-3

class activation_coefficient:
    """ This class defines the activation coefficients of respective receptor-ligand binding.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | TGFb_OBu         |:math:`k^{TGFb}_{OB_u}`    |  pM     |
    +------------------+---------------------------+---------+
    | TGFb_OCa         |:math:`k^{TGFb}_{OC_a}`    |  pM     |
    +------------------+---------------------------+---------+
    | PTH_OB           | :math:`k^{PTH}_{OB}`      |  pM     |
    +------------------+---------------------------+---------+
    | RANKL_RANK       |:math:`k_{RANK}^{RANKL}`   |  pM     |
    +------------------+---------------------------+---------+
    | MCSF_OCu         |:math:`k_{OC_u}^{MCSF}`    |  pM     |
    +------------------+---------------------------+---------+

    :param TGFb_OBu: parameter for TGF-beta binding on OBp
    :type TGFb_OBu: float
    :param TGFb_OCa: parameter for TGF-beta binding on OCa
    :type TGFb_OCa: float
    :param PTH_OB: parameter for PTH binding to osteoblasts (activation)
    :type PTH_OB: float
    :param RANKL_RANK: parameter for RANKL binding on RANK
    :type RANKL_RANK: float
    :param MCSF_OCu: parameter for MCSF binding on OCu
    :type MCSF_OCu: float"""
    def __init__(self):
        self.TGFb_OBu = 0.000563278809675429
        self.TGFb_OCa = 0.000563278809675429
        self.PTH_OB = 150
        self.RANKL_RANK = 16.65
        self.MCSF_OCu = 0.001


class repression_coefficient:
    """ This class defines the repression coefficients of respective receptor-ligand binding.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | TGFb_OBp         |:math:`k^{TGFb}_{OBp}`     |  pM     |
    +------------------+---------------------------+---------+
    | PTH_OB           | :math:`k^{PTH}_{OB}`      |  pM     |
    +------------------+---------------------------+---------+

    :param TGFb_OBp: parameter for TGF-beta binding on OBp
    :type TGFb_OBp: float
    :param PTH_OB: parameter for PTH binding on osteoblasts (repressor)
    :type PTH_OB: float"""
    def __init__(self):
        self.TGFb_OBp = 0.00189
        self.PTH_OB = 0.222581427709954


class degradation_rate:
    r""" This class defines the degradation rates of the different factors.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+------------------------------+---------+
    | Parameter Name   | Symbol                       | Units   |
    +==================+==============================+=========+
    | PTH              |:math:`{D}_{PTH}`             |  1/day  |
    +------------------+------------------------------+---------+
    | OPG              |:math:`{D}_{OPG}`             |  1/day  |
    +------------------+------------------------------+---------+
    | RANKL            |:math:`{D}_{RANKL}`           |  1/day  |
    +------------------+------------------------------+---------+
    | TGFb             |:math:`{D}_{TGF-\beta}`       |  1/day  |
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
        self.PTH = 86
        self.OPG = 3.50e-1
        self.RANKL = 1.0132471014805027e+1
        self.TGFb = 2


class concentration:
    """ This class defines fixed concentrations.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | OPG_max          |:math:`{OPG}_{sat}`        |  pM     |
    +------------------+---------------------------+---------+
    | MCSF             |-                          |  pM     |
    +------------------+---------------------------+---------+
    | RANK             |-                          |  pM     |
    +------------------+---------------------------+---------+

    :param OPG_max: maximum concentration of OPG
    :type OPG_max: float
    :param MCSF: concentration of MCSF
    :type MCSF: float
    :param RANK: concentration of RANK
    :type RANK: float
    """
    def __init__(self):
        self.OPG_max = 2.00e+8
        self.MCSF = 0.001
        self.RANK = 1.00e+1

class binding_constant:
    """ This class defines the binding constants of RANK RANKL and OPG.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+---------------------------+---------+
    | Parameter Name   | Symbol                    | Units   |
    +==================+===========================+=========+
    | RANKL_OPG        |:math:`k^{RANKL}_{OPG}`    |1/pM     |
    +------------------+---------------------------+---------+
    | RANKL_RANK       |:math:`k^{RANKL}_{RANK}`   |1/pM     |
    +------------------+---------------------------+---------+

    :param RANKL_OPG: binding constant for RANKL-OPG
    :type RANKL_OPG: float
    :param RANKL_RANK: binding constant for RANKL-RANK
    :type RANKL_RANK: float"""
    def __init__(self):
        self.RANKL_OPG = 1.00e-3
        self.RANKL_RANK = 3.411764705882353e-002

class production_rate:
    r""" This class defines the intrinsic/ endogenous production rates of the different factors.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------+-----------------------------+---------+
    | Parameter Name       | Symbol                      | Units   |
    +======================+=============================+=========+
    | intrinsic_PTH        |:math:`P_{PTH}`              |  pM/day |
    +----------------------+-----------------------------+---------+
    | intrinsic_RANKL      |:math:`\beta^{RANKL}_{OB_p}` |  pM/day |
    +----------------------+-----------------------------+---------+
    | min_OPG_per_cell     |:math:`\beta^{OPG}_{OB_a}`   |  pM     |
    +----------------------+-----------------------------+---------+
    | bool_OBp_produce_OPG | code-specific - set to 0    |  -      |
    +----------------------+-----------------------------+---------+
    | bool_OBa_produce_OPG | code-specific - set to 1    |  -      |
    +----------------------+-----------------------------+---------+
    | max_RANKL_per_cell   |:math:`N_{RANKL}^{OB}`       |  pM     |
    +----------------------+-----------------------------+---------+
    | max_RANK_per_cell    |:math:`N^{RANK}_{OC_p}`      |  pM     |
    +----------------------+-----------------------------+---------+
    |bool_OBp_produce_RANKL| code-specific - set to 1    |  -      |
    +----------------------+-----------------------------+---------+
    |bool_OBa_produce_RANKL| code-specific - set to 0    |  -      |
    +----------------------+-----------------------------+---------+

    :param intrinsic_PTH: systemic concentration of PTH
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
        self.intrinsic_PTH = 2.907
        self.intrinsic_RANKL = 1.684195714712206e+5
        self.min_OPG_per_cell = 1.624900337835679e+008
        self.bool_OBp_produce_OPG = 0
        self.bool_OBa_produce_OPG = 1
        self.max_RANKL_per_cell = 27e+5
        self.max_RANK_per_cell = 1.000e+004
        self.bool_OBp_produce_RANKL = 1  # 1=yes
        self.bool_OBa_produce_RANKL = 0  # 0=no


class bone_volume:
    r""" This class defines the parameters relevant for bone volume of the bone model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------+----------------------------+---------+
    | Parameter Name       | Symbol                     | Units   |
    +======================+============================+=========+
    | formation_rate       |:math:`k_{form}`            |  1/day  |
    +----------------------+----------------------------+---------+
    | resorption_rate      |:math:`k_{res}`             |  1/day  |
    +----------------------+----------------------------+---------+
    |stored_TGFb_content   |:math:`n_{TGF\beta}^{bone}` |  pM     |
    +----------------------+----------------------------+---------+

    :param formation_rate: formation rate of bone volume
    :type formation_rate: float
    :param resorption_rate: resorption rate of bone volume
    :type resorption_rate: float
    :param stored_TGFb_content: proportionality constant expressing the TGF-β content stored in bone volume
    :type stored_TGFb_content: float """
    def __init__(self):
        self.formation_rate = 40.0
        self.resorption_rate = 200.0
        self.stored_TGFb_content = 0.01


class mechanics:
    r""" This class defines the parameters relevant for mechanics of the bone model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +-------------------------------------+--------------------------------------+---------+
    | Parameter Name                      | Symbol                               | Units   |
    +=====================================+======================================+=========+
    | strain_effect_on_OBp_steady_state   |:math:`\mu(r,t)`                      |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | strain_effect_on_OBp                |:math:`\mu(r,t)`                      |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | strain_energy_density               |:math:`\check{\psi}(r,t))`            |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | update_OBp_proliferation_rate       | -                                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | RANKL_production                    |:math:`\beta_{RANKL}^{mech}`          |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | unit_tensor_as_matrix               |:math:`\mathbb{I}`                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | deviatoric_part_of_unit_tensor      |:math:`\mathbb{K}`                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | stiffness_tensor_vascular_pores     |:math:`\mathbb{c}^{micro}_{vas}`      |  GPa    |
    +-------------------------------------+--------------------------------------+---------+
    | stiffness_tensor_bone_matrix        |:math:`\mathbb{c}^{micro}_{bm}`       |  GPa    |
    +-------------------------------------+--------------------------------------+---------+
    |step_size_for_Hill_tensor_integration| -                                    |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | hill_tensor_cylindrical_inclusion   |:math:`\mathbb{P}_{r}^{bm}`           |        -|
    +-------------------------------------+--------------------------------------+---------+
    | stress_tensor_normal_loading        |:math:`\sigma^{tissue}`               |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | biomech_transduction_strength       |:math:`\lambda`                       |  -      |
    +-------------------------------------+--------------------------------------+---------+
    | biomech_transduction_strength_RANKL |:math:`\kappa`                        | pM/day  |
    +-------------------------------------+--------------------------------------+---------+
    | correction_factor                   |:math:`K`                             |  GPa    |
    +-------------------------------------+--------------------------------------+---------+
    """
    def __init__(self):
        self.strain_effect_on_OBp_steady_state = None
        self.strain_energy_density_steady_state = None
        self.update_OBp_proliferation_rate = True
        self.RANKL_production = 0
        self.unit_tensor_as_matrix = np.array([[1, 0, 0, 0, 0, 0],
                                               [0, 1, 0, 0, 0, 0],
                                               [0, 0, 1, 0, 0, 0],
                                               [0, 0, 0, 1, 0, 0],
                                               [0, 0, 0, 0, 1, 0],
                                               [0, 0, 0, 0, 0, 1]])
        self.stiffness_tensor_vascular_pores = np.array([[1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [1, 1, 1, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0],
                                                       [0, 0, 0, 0, 0, 0]]) * 2.3
        self.stiffness_tensor_bone_matrix = np.array([[18.5, 10.3, 10.4, 0, 0, 0],
                                                    [10.3, 20.8, 11.0, 0, 0, 0],
                                                    [10.4, 11.0, 28.4, 0, 0, 0],
                                                    [0, 0, 0, 12.9, 0, 0],
                                                    [0, 0, 0, 0, 11.5, 0],
                                                    [0, 0, 0, 0, 0, 9.3]])
        self.step_size_for_Hill_tensor_integration = 2 * np.pi / 50
        self.hill_tensor_cylindrical_inclusion = None
        self.stress_tensor_normal_loading = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -30]]) * (10 ** -3)
        self.biomech_transduction_strength = 0.5
        self.biomech_transduction_strength_RANKL = 18
        self.correction_factor = 1.0e-6


class calibration:
    r""" This class defines the parameters relevant for calibration of the bone model.
    The following table provides a mapping between the model parameters and the original names from the publication:

    +----------------------+-----------------------+---------+
    | Parameter Name       | Symbol                | Units   |
    +======================+=======================+=========+
    | turnover             |:math:`\tau_{turnover}`|  -      |
    +----------------------+-----------------------+---------+
    | steady_state_turnover|:math:`\delta`         |  -      |
    +----------------------+-----------------------+---------+
    | OCa                  |:math:`\beta`          |  -      |
    +----------------------+-----------------------+---------+
    | OBa                  |:math:`\gamma`         |  -      |
    +----------------------+-----------------------+---------+

    :param turnover: calibration for turnover to specific surface
    :type turnover: float
    :param steady_state_turnover: calibration for steady state turnover rate of the bone model
    :type steady_state_turnover: float
    :param OCa: calibration coefficient for OCa of the bone model
    :type OCa: float
    :param OBa: calibration coefficient for OBa of the bone model
    :type OBa: float
    """
    def __init__(self):
        # not stated
        # self.turnover = 0.00395
        self.turnover = 5.961e-03
        self.steady_state_turnover = 0.395
        self.OCa = 0.09
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
    :param proliferation_rate: proliferation rates of the different cell types
    :type proliferation_rate: proliferation_rate
    :param bone_volume: parameters relevant for bone volume of the bone model
    :type bone_volume: bone_volume
    :param calibration: parameters relevant for calibration of the bone model
    :type calibration: calibration
    """
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
