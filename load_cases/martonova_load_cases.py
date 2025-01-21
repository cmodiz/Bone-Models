class injected_PTH_pulse:
    def __init__(self):
        self.on_duration = None
        self.off_duration = None
        self.max = None


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
        self.basal_PTH_pulse.off_duration = 6.4
        # -> tau_on
        self.basal_PTH_pulse.on_duration = 4.2
        # -> T
        self.basal_PTH_pulse.period = self.basal_PTH_pulse.off_duration + self.basal_PTH_pulse.on_duration
        # drug
        self.drug_dose = None
        self.injection_frequency = None
        self.injected_PTH_pulse = injected_PTH_pulse()

