class Calcium:
    r"""
    Class to hold parameters for the calcium grid.

    +-----------------------+-----------------------------------------------+-----------------------------+
    | **Parameter (Type)**  | **Description**                               | **Symbol in Paper**         |
    +=======================+===============================================+=============================+
    | ``minimum_content``   | Minimum calcium content in % wt.              | :math:`c_{0}`               |
    | *(float)*             |                                               |                             |
    +-----------------------+-----------------------------------------------+-----------------------------+
    | ``maximum_content``   | Maximum calcium content in % wt.              | :math:`c_{\max}`            |
    | *(float)*             |                                               |                             |
    +-----------------------+-----------------------------------------------+-----------------------------+
    """
    def __init__(self):
        self.minimum_content = 0
        self.maximum_content = 31


class Reference_BMDD:
    r""" Class to hold reference BMDD parameters for skewed Gaussian distribution if model is initialized with BMDD.
    The parameters were taken from the original paper and adjusted to fit a realistic BMDD.
    In real life, this should be initialized with a BMDD/ characteristics measured from an image.

    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | **Parameter (Type)**                         | **Description**                               | **Symbol in Paper**         |
    +==============================================+===============================================+=============================+
    | ``standard_deviation_lower_percentile``      | Standard deviation of the lower percentiles   | :math:`\sigma_{1}`          |
    | *(float)*                                    | of the BMDD in % wt.                          |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | ``standard_deviation_upper_percentile``      | Standard deviation of the upper percentiles   | :math:`\sigma_{2}`          |
    | *(float)*                                    | of the BMDD in % wt.                          |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | ``peak``                                     | Peak position of the BMDD in % wt.            | :math:`Ca_{peak}`           |
    | *(float)*                                    |                                               |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    """
    def __init__(self):
        self.standard_deviation_lower_percentile = 4.93 / 1.5
        self.standard_deviation_upper_percentile = 5.55 / 1.5
        self.peak = 22.94


class Rate:
    r""" Class to hold reference BMDD parameters for formation and resorption rates. The parameters are initialized here,
    but can be changed later depending on the initialized bone volume.

    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | **Parameter (Type)**                         | **Description**                               | **Symbol in Paper**         |
    +==============================================+===============================================+=============================+
    | ``initial_formation``                        | Initial formation rate in :math: `mm^3/ year` | :math:`\\j_{OB}`            |
    | *(float)*                                    |                                               |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | ``final_formation``                          | Final formation rate in :math: `mm^3/ year`   | :math:`\\j_{OB}(t)`         |
    | *(float)*                                    |                                               |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | ``initial_resorption``                       | Initial resorption rate in 1/year             | :math:`\omega_{OC}(c)`      |
    | *(float)*                                    |                                               |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    | ``final_resorption``                         | Final resorption rate in 1/year               | :math:`\omega_{OC}(c,t)`    |
    | *(float)*                                    |                                               |                             |
    +----------------------------------------------+-----------------------------------------------+-----------------------------+
    """
    def __init__(self):
        self.initial_formation = 0.1
        self.final_formation = 0.1
        self.initial_resorption = 0.1
        self.final_resorption = 0.1


class Mineralization_Law:
    r"""
    Class representing the mineralization law parameters used in the mineralization model.

    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    | **Parameter (Type)**                    | **Description**                               | **Symbol in Paper**         |
    +=========================================+===============================================+=============================+
    | ``turnover_time`` *(float)*             | Turnover time of bone packets in years        | :math:`t_{TO}`              |
    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    | ``primary_mineral_content`` *(float)*   | Primary mineral content in % wt.              | :math:`c_{1}`               |
    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    | ``maximum_mineral_content`` *(float)*   | Secondary mineral content contribution in % wt| :math:`c_{2}`               |
    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    | ``primary_apposition_rate`` *(float)*   | Time constant of the primary mineralization   | :math:`\tau_{1}`            |
    |                                         | process in years                              |                             |
    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    | ``secondary_apposition_rate`` *(float)* | Time constant of the secondary mineralization | :math:`\tau_{2}`            |
    |                                         | process years                                 |                             |
    +-----------------------------------------+-----------------------------------------------+-----------------------------+
    """
    def __init__(self):
        self.turnover_time = 5  # [years]
        self.primary_mineral_content = 11.9
        self.maximum_mineral_content = 31.1 - self.primary_mineral_content
        self.primary_apposition_rate = 150
        self.secondary_apposition_rate = 0.2


class Bone_Volume:
    r"""
    Class storing the bone volume parameters required for initializing
    the mineralization model.

    +-----------------------------------+-------------------------------------------+-------------------------+
    | **Parameter (Type)**              | **Description**                           | **Symbol in Paper**     |
    +===================================+===========================================+=========================+
    | ``initial_value`` *(float)*       | Initial bone volume used when initializing| BV                      |
    |                                   | the model from a mineralization law       |                         |
    |                                   | in :math:`mm^3`                           |                         |
    +-----------------------------------+-------------------------------------------+-------------------------+
    """
    def __init__(self):
        self.initial_value = 2.4  # [mm^3]


class Ruffoni_Parameters:
    r"""
    Class to hold all parameters for the Ruffoni mineralization model.
    """
    def __init__(self):
        self.calcium = Calcium()
        self.rate = Rate()
        self.reference_bmdd = Reference_BMDD()
        self.mineralization_law = Mineralization_Law()
        self.bone_volume = Bone_Volume()
