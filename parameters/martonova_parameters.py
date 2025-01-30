class kinematics:
    """  This class contains the kinematic parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name        | Symbol    | Units
    +==================+===========+==================+
    | receptor_desensitized | k_1       |    [1/min]
    | receptor_resensitized | k_{-1}    |    [1/min]
    | complex_desensitized  | k_2       |    [1/min]
    | complex_resensitized  | k_{-2}    |      [1/min]
    | active_complex_binding | k_r       |  [-]
    | active_complex_unbinding | k_{-r}    |    [-]
    | inactive_complex_binding | k_d       |    [-]
    | receptor              | K_1           | [-]
    | complex               | K_2           | [-]
    | active_binding_unbinding | K_r        | [nM]
    | inactive_complex_unbinding | k_{-d}    | [1/h]
    | inactive_binding_unbinding | K_d       | [nM]

    :param receptor_desensitized: rate of receptor desensitization
    :type receptor_desensitized: float
    :param receptor_resensitized: rate of receptor resensitization
    :type receptor_resensitized: float
    :param complex_desensitized: rate of complex desensitization
    :type complex_desensitized: float
    :param complex_resensitized: rate of complex resensitization
    :type complex_resensitized: float
    :param active_complex_binding: rate of active complex binding
    :type active_complex_binding: float
    :param active_complex_unbinding: rate of active complex unbinding
    :type active_complex_unbinding: float
    :param inactive_complex_binding: rate of inactive complex binding
    :type inactive_complex_binding: float
    :param receptor: combined kinetic receptor constant
    :type receptor: float
    :param complex: combined kinetic complex constant
    :type complex: float
    :param active_binding_unbinding: active binding unbinding constant
    :type active_binding_unbinding: float
    :param inactive_complex_unbinding: rate of inactive complex unbinding
    :type inactive_complex_unbinding: float
    :param inactive_binding_unbinding: inactive binding unbinding constant
    :type inactive_binding_unbinding: float """
    def __init__(self):
        """ Constructor method. """
        # -> k_1
        self.receptor_desensitized = 0.012
        # -> k_{-1}
        self.receptor_resensitized = 0.104
        # k_2
        self.complex_desensitized = 0.222
        # k_{-2}
        self.complex_resensitized = 0.055
        # k_r
        self.active_complex_binding = 1
        # k_{-r}
        self.active_complex_unbinding = 1000
        # k_d
        self.inactive_complex_binding = 1
        # K_1
        self.receptor = self.receptor_resensitized / self.receptor_desensitized
        # K_2
        self.complex = self.complex_resensitized / self.complex_desensitized
        # K_r
        self.active_binding_unbinding = self.active_complex_unbinding / self.active_complex_binding
        # k_{-d}
        self.inactive_complex_unbinding = self.active_binding_unbinding / (self.receptor / self.complex) * self.active_complex_binding
        # K_d
        self.inactive_binding_unbinding = self.inactive_complex_unbinding / self.inactive_complex_binding


class activity:
    """ This class contains the activity parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name        | Symbol    | Units
    +==================+===========+==================+
    | active_receptor        | a_1       |    [nM]
    | active_complex         | a_2       |    [nM]
    | inactive_complex       | a_3       |    [nM]
    | inactive_receptor      | a_4       |    [nM]

    :param active_receptor: activity constant of active receptor
    :type active_receptor: float
    :param active_complex: activity constant of active complex
    :type active_complex: float
    :param inactive_complex: activity constant of inactive complex
    :type inactive_complex: float
    :param inactive_receptor: activity constant of inactive receptor
    :type inactive_receptor: float """
    def __init__(self):
        """ Constructor method. """
        # -> a_1
        self.active_receptor = 100 * 20
        # -> a_2
        self.active_complex = None  # computed in init method
        # -> a_3
        self.inactive_complex = 100 * 10
        # -> a_4
        self.inactive_receptor = 100 * 1


class basal_PTH_pulse:
    """ This class contains the basal PTH pulse parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name        | Symbol    | Units
    +==================+===========+==================+
    | min                   | gamma_off |    [nM]
    | max                   | gamma_on  |    [nM]
    | off_duration          | tau_off   |    [min]
    | on_duration           | tau_on    |    [min]
    | period                | T         |    [min]

    :param min: concentration of non-pulsatile share of PTH pulse
    :type min: float
    :param max: concentration of pulsatile share of PTH pulse
    :type max: float
    :param off_duration: duration of non-pulsatile share of PTH pulse
    :type off_duration: float
    :param on_duration: duration of pulsatile share of PTH pulse
    :type on_duration: float
    :param period: duration of the PTH pulse period (pulse + non-pulse)
    :type period: float """
    def __init__(self):
        """ Constructor method. """
        # -> gamma_off
        self.min = None
        # -> gamme_on
        self.max = None
        # -> tau_off
        self.off_duration = None
        # -> tau_on
        self.on_duration = None
        # -> T
        self.period = None


class injected_PTH_pulse:
    """ This class contains the injected PTH pulse parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name        | Symbol    | Units
    +==================+===========+==================+
    | min                   | gamma_off |    [nM]
    | max                   | gamma_on  |    [nM]
    | on_duration           | tau_on    |    [min]

    :param min: concentration of non-pulsatile share of PTH pulse
    :type min: float
    :param max: concentration of pulsatile share of PTH pulse
    :type max: float
    :param on_duration: duration of pulsatile share of PTH pulse
    :type on_duration: float """
    def __init__(self):
        """ Constructor method. """
        # -> gamma_off
        self.min = None
        # -> gamma_on
        self.max = None
        # -> tau_on
        self.on_duration = None


class pharmacokinetics:
    """ This class contains the pharmacokinetic parameters of the two-state receptor model by Martonova et al. for injected PTH.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name        | Symbol    | Units
    +==================+===========+==================+
    | absorption_rate       | k_a        |    [1/h]
    | elimination_rate      | k_e        |    [1/h]
    | volume_of_distribution| V         |    [L]
    | bioavailability       | F         |    [-]

    :param absorption_rate: rate of absorption of the injected PTH
    :type absorption_rate: float
    :param elimination_rate: rate of elimination of the injected PTH
    :type elimination_rate: float
    :param volume_of_distribution: volume of distribution of the injected PTH
    :type volume_of_distribution: float
    :param bioavailability: bioavailability of the injected PTH
    :type bioavailability: float """
    def __init__(self):
        """ Constructor method. """
        # -> ka
        self.absorption_rate = 4.38  # [1/h]
        # -> ke
        self.elimination_rate = 0.693  # [1/h]
        # -> V
        self.volume_of_distribution = 110  # [L]
        # F
        self.bioavailability = 0.95  # [-]


class Parameters:
    """ This class contains all parameters of the two-state receptor model by Martonova et al.

    :param kinematics: kinematic parameters of the model, see :class:`kinematics` for details
    :type kinematics: kinematics
    :param activity: activity parameters of the model, see :class:`activity` for details
    :type activity: activity
    :param basal_PTH_pulse: basal PTH pulse parameters of the model, see :class:`basal_PTH_pulse` for details
    :type basal_PTH_pulse: basal_PTH_pulse
    :param injected_PTH_pulse: injected PTH pulse parameters of the model, see :class:`injected_PTH_pulse` for details
    :type injected_PTH_pulse: injected_PTH_pulse
    :param pharmacokinetics: pharmacokinetic parameters of the model, see :class:`pharmacokinetics` for details
    :type pharmacokinetics: pharmacokinetics """

    def __init__(self):
        self.kinematics = kinematics()
        self.activity = activity()
        self.basal_PTH_pulse = basal_PTH_pulse()
        self.injected_PTH_pulse = injected_PTH_pulse()
        self.pharmacokinetics = pharmacokinetics()
