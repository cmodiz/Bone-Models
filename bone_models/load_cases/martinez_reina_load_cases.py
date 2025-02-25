import numpy as np


class Martinez_Reina_Load_Case:
    """ Load cases of the Martinez-Reina 2019 model """

    def __init__(self):
        self.start_time = 1
        self.end_time = 1000 + 2*365
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

        self.start_postmenopausal_osteoporosis = 0
        self.end_postmenopausal_osteoporosis = 1000 + 2*365

        self.start_denosumab_treatment = self.start_postmenopausal_osteoporosis + 365  # leave 1 year untreated
        self.end_denosumab_treatment = self.start_postmenopausal_osteoporosis + 365 * 2
        self.treatment_period = 183   # every 6 months
        self.denosumab_dose = 60 * (10**6) / 60  # [ng]/kg body weight (60 kg as reference body weight) every 6months