import numpy as np


class Lerebours_Load_Case:
    def __init__(self):
        self.start_time = 50
        self.end_time = 2050
        self.stress_tensor = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -25]]) * 10 ** -3   # [GPa]
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.TGFb_injection = 0
