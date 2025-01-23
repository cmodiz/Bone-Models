from parameters.lemaire_parameters import Parameters as Lemaire_Parameters


class Calibration_Parameters:
    def __init__(self):
        self.cellular_responsiveness = 0.030259870370592704
        self.integrated_activity = 0.0007172096391750288


class Parameters(Lemaire_Parameters):
    def __init__(self):
        super().__init__()
        self.calibration_parameters = Calibration_Parameters()
