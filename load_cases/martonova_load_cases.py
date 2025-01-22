class injected_PTH_pulse:
    def __init__(self):
        self.max = None
        self.on_duration = None
        self.off_duration = None
        self.period = None


class Basal_PTH_pulse:
    def __init__(self):
        self.min = None
        self.max = None
        self.off_duration = None
        self.on_duration = None
        self.period = None


class Healthy:
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


class Hyperparathyroidism:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.01381
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00977
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 7.6
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 3.5
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Osteoporosis:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.003321
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.0016967
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 5.2
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 24.6
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Postmenopausal_Osteoporosis:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.9
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 0.9
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 6.4
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 4.2
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Hypercalcemia:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.25
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 0.12
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 6.4 / 0.68
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 4.2 / 0.68
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Hypocalcemia:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 2.57
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 13.1
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 6.4 / 1.96
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 4.2 / 1.96
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()


class Glucocorticoid_Induced_Osteoporosis:
    def __init__(self):
        self.basal_PTH_pulse = Basal_PTH_pulse()
        # -> gamma_off
        self.basal_PTH_pulse.min = 0.00332 * 0.48
        # -> gamma_on
        self.basal_PTH_pulse.max = 0.00276 * 1.75
        # -> tau_off
        self.basal_PTH_pulse.off_duration = 6.4 * 0.95
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 4.2 * 0.95
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()

