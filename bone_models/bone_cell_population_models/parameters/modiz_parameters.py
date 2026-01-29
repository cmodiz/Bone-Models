from .lemaire_parameters import Lemaire_Parameters


class Calibration_Parameters:
    def __init__(self):
        self.cellular_responsiveness = 0.030259870370592704
        self.integrated_activity = 0.0007172096391750288


class Calibration_Parameters_Only_For_Healthy_State:
    def __init__(self):
        self.cellular_responsiveness = 0.02066869136769719
        self.integrated_activity = 0.0005288141739266334


class Modiz_Parameters(Lemaire_Parameters):
    def __init__(self):
        super().__init__()
        self.calibration = Calibration_Parameters()
        self.calibration_only_for_healthy_state = Calibration_Parameters_Only_For_Healthy_State()
