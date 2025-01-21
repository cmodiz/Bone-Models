class Healthy:
    def __init__(self):
        # -> gamma_off
        self.basal_PTH_pulse_min = 0.00332
        # -> gamme_on
        self.basal_PTH_pulse_max = 0.00276
        # -> tau_off
        self.basal_PTH_pulse_off_duration = 6.4
        # -> tau_on
        self.basal_PTH_pulse_on_duration = 4.2
        # -> T
        self.basal_PTH_pulse_period = None
        # drug
        self.drug_dose = None
        self.injection_frequency = None

