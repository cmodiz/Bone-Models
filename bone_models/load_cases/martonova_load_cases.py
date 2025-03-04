class injected_PTH_pulse:
    """ This class contains the injected PTH pulse parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +------------------+-----------+------------------+
    | Parameter Name   | Symbol    | Units            |
    +==================+===========+==================+
    | max              | gamma_on  |    [nM]          |
    +------------------+-----------+------------------+
    | off_duration     | tau_off   |    [min]         |
    +------------------+-----------+------------------+
    | on_duration      | tau_on    |    [min]         |
    +------------------+-----------+------------------+
    | period           |   T       |    [min]         |
    +------------------+-----------+------------------+

    :param max: concentration of pulsatile share of injected PTH pulse
    :type max: float
    :param off_duration: duration of non-pulsatile share of injected PTH pulse
    :type off_duration: float
    :param on_duration: duration of pulsatile share of injected PTH pulse
    :type on_duration: float
    :param period: duration of the injected PTH pulse period (pulse + non-pulse)
    :type period: float """
    def __init__(self):
        self.max = None
        self.on_duration = None
        self.off_duration = None
        self.period = None


class Basal_PTH_pulse:
    """ This class contains the basal PTH pulse parameters of the two-state receptor model by Martonova et al.

    The following table provides a mapping between the model parameters
    and their original names from the publication:

    +-----------------------+-----------+-------------+
    | Parameter Name        | Symbol    | Units       |
    +=======================+===========+=============+
    | min                   | gamma_off |    [nM]     |
    +-----------------------+-----------+-------------+
    | max                   | gamma_on  |    [nM]     |
    +-----------------------+-----------+-------------+
    | off_duration          | tau_off   |    [min]    |
    +-----------------------+-----------+-------------+
    | on_duration           | tau_on    |    [min]    |
    +-----------------------+-----------+-------------+
    | period                | T         |    [min]    |
    +-----------------------+-----------+-------------+

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
        self.min = None
        self.max = None
        self.off_duration = None
        self.on_duration = None
        self.period = None


class Martonova_Healthy:
    """ This class contains the parameters for the healthy state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 4.2
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 6.4
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Hyperparathyroidism:
    """ This class contains the parameters for the hyperparathyroidism state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """

    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.01381
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00977
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 3.5
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 7.6
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Osteoporosis:
    """ This class contains the parameters for the osteoporosis state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.003321
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.0016967
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 24.6
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 5.2
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Postmenopausal_Osteoporosis:
    """ This class contains the parameters for the postmenopausal osteoporosis state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected_PTH_pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.9
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 0.9
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 4.2
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 6.4
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Hypercalcemia:
    """ This class contains the parameters for the hypercalcemia state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.25
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 0.12
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 4.2 / 0.68
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 6.4 / 0.68
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Hypocalcemia:
    """ This class contains the parameters for the hypocalcemia state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 2.57
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 13.1
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 4.2 / 1.96
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 6.4 / 1.96
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Glucocorticoid_Induced_Osteoporosis:
    """ This class contains the parameters for the glucocorticoid-induced osteoporosis state without injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.48
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 1.75
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 4.2 * 0.95
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 6.4 * 0.95
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Martonova_Hyperparathyroidism_With_Drug:
    """ This class contains the parameters for the hyperparathyroidism state with injection of the two-state receptor model by Martonova et al.

    :param basal_PTH_pulse: basal PTH pulse parameters
    :type basal_PTH_pulse: Basal_PTH_pulse
    :param drug_dose: dose of the drug in micrograms
    :type drug_dose: float
    :param injection_frequency: frequency of the drug injection (e.g. 24h)
    :type injection_frequency: float
    :param injected_PTH_pulse: injected PTH pulse parameters
    :type injected_PTH_pulse: injected_PTH_pulse """
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.01381
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00977
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 3.5
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 7.6
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = 20  # micrograms
        self.injection_frequency = 24  # injections per day
        self.injected_PTH_pulse = injected_PTH_pulse()
