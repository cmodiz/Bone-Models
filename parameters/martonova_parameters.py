class kinematics:
    def __init__(self):
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
        self.receptor = self.receptor_desensitized / self.receptor_resensitized
        # K_2
        self.complex = self.complex_desensitized / self.complex_resensitized
        # K_r
        self.active_binding_unbinding = self.active_complex_unbinding / self.active_complex_binding
        # k_{-d}
        self.inactive_complex_unbinding = self.active_binding_unbinding / (self.receptor / self.complex) * self.active_complex_binding
        # K_d
        self.inactive_binding_unbinding = self.inactive_complex_unbinding / self.inactive_complex_binding


class activity:
    def __init__(self):
        # -> a_1
        self.active_receptor = 100 * 20
        # -> a_2
        self.active_complex = None  # computed in init method
        # -> a_3
        self.inactive_complex = 100 * 10
        # -> a_4
        self.inactive_receptor = 100 * 1


class basal_PTH_pulse:
    def __init__(self):
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
    def __init__(self):
        # -> gamma_off
        self.min = None
        # -> gamme_on
        self.max = None
        # -> tau_on
        self.on_duration = None


class Parameters:
    def __init__(self):
        self.kinematics = kinematics()
        self.activity = activity()
        self.basal_PTH_pulse = basal_PTH_pulse()
        self.injected_PTH_pulse = injected_PTH_pulse()
