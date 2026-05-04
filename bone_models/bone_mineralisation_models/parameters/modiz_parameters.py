class bone_volume:
    r""" This class defines the parameters relevant for bone volume.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------+--------------------------+---------+
    | Parameter Name       | Symbol                   | Units   |
    +======================+==========================+=========+
    | formation_rate       |:math:`\mathrm{BFR}`      |  1/day  |
    +----------------------+--------------------------+---------+
    | resorption_rate      |:math:`\mathrm{BRs.R}`    |  1/day  |
    +----------------------+--------------------------+---------+

    :param formation_rate: formation rate of bone volume
    :type formation_rate: float
    :param resorption_rate: resorption rate of bone volume
    :type resorption_rate: float """
    def __init__(self):
        self.formation_rate = 9.011
        self.resorption_rate = 566.7


class mineralisation:
    r"""
    This class defines the parameters relevant for mineralisation of the bone model.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +----------------------------+---------------------------+---------+
    | Parameter Name             | Symbol                    | Units   |
    +============================+===========================+=========+
    | density_organic            | :math:`\rho_{O}`          | g/cm³   |
    +----------------------------+---------------------------+---------+
    | density_mineral            | :math:`\rho_{M}`          | g/cm³   |
    +----------------------------+---------------------------+---------+
    | volume_fraction_organic    | :math:`f_O`               | -       |
    +----------------------------+---------------------------+---------+
    | lag_time                   | :math:`Mlt`               | days    |
    +----------------------------+---------------------------+---------+
    | primary_phase_duration     | :math:`Mpt`               | days    |
    +----------------------------+---------------------------+---------+
    | primary_mineral_content    | :math:`M_{prim}`          | -       |
    +----------------------------+---------------------------+---------+
    | maximum_mineral_content    | :math:`M_{max}`           | -       |
    +----------------------------+---------------------------+---------+
    | length_of_queue            | :math:`n`                 | -       |
    +----------------------------+---------------------------+---------+
    | reference_apposition_rate  | :math:`\mathrm{MAR}_{ref}`| -       |
    +----------------------------+---------------------------+---------+
    | resorption_bias_age        | :math:`t_{pref}`          | -       |
    +----------------------------+---------------------------+---------+

    :param density_organic: Density of organic material.
    :type density_organic: float
    :param density_mineral: Density of mineral material.
    :type density_mineral: float
    :param volume_fraction_organic: Volume fraction of organic material.
    :type volume_fraction_organic: float
    :param lag_time: Lag time before mineralisation begins.
    :type lag_time: float
    :param primary_phase_duration: Duration of the primary mineralisation phase.
    :type primary_phase_duration: float
    :param primary_mineral_content: Mineral content after the primary phase.
    :type primary_mineral_content: float
    :param maximum_mineral_content: Maximum achievable mineral content.
    :type maximum_mineral_content: float
    :param length_of_queue: Length of the ageing queue.
    :type length_of_queue: int
    :param reference_apposition_rate: Reference mineral apposition rate.
    :type reference_apposition_rate: float
    :param resorption_bias_age: Preferred age for resorption bias.
    :type resorption_bias_age: float
    """
    def __init__(self):
        self.density_organic = 1.1
        self.density_mineral = 3.2
        self.volume_fraction_organic = 0.34
        self.lag_time = 12
        self.primary_phase_duration = 10
        self.primary_mineral_content = 0.121
        self.maximum_mineral_content = 0.442
        self.length_of_queue = 30000
        self.reference_apposition_rate = 0.00039772586198511
        self.resorption_bias_age = 1000

class Modiz_Parameters:
    def __init__(self):
        self.mineralisation = mineralisation()
        self.bone_volume = bone_volume()
