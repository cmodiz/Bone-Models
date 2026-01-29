import numpy as np


class Martinez_Reina_Load_Case:
    """ Load cases of the Martinez-Reina 2019 model. Mechanical loading is applied to the bone (compression with 5 MPa, same as habitual loading),
    external injections like the Lemaire et al., 2004 model can be given (but are zero in this case), and postmenopausal
    osteoporosis (from beginning to end) and denosumab treatment (from one year onwards with a treatment period of 183 days) are simulated.

    :param start_time: start time of the mechanical loading
    :type start_time: float
    :param end_time: end time of the mechanical loading
    :type end_time: float
    :param stress_tensor: stress tensor applied to the bone
    :type stress_tensor: numpy.ndarray
    :param OBp_injection: external injection of precursor osteoblasts
    :type OBp_injection: float
    :param OBa_injection: external injection of active osteoblasts
    :type OBa_injection: float
    :param OCa_injection: external injection of active osteoclasts
    :type OCa_injection: float
    :param PTH_injection: external injection of parathyroid hormone
    :type PTH_injection: float
    :param OPG_injection: external injection of osteoprotegerin
    :type OPG_injection: float
    :param RANKL_injection: external injection of receptor activator of nuclear factor kappa-B ligand
    :type RANKL_injection: float
    :param TGFb_injection: external injection of transforming growth factor beta
    :type TGFb_injection: float
    :param start_postmenopausal_osteoporosis: start time of postmenopausal osteoporosis
    :type start_postmenopausal_osteoporosis: float
    :param end_postmenopausal_osteoporosis: end time of postmenopausal osteoporosis
    :type end_postmenopausal_osteoporosis: float
    :param start_denosumab_treatment: start time of denosumab treatment
    :type start_denosumab_treatment: float
    :param end_denosumab_treatment: end time of denosumab treatment
    :type end_denosumab_treatment: float
    :param treatment_period: period of denosumab treatment
    :type treatment_period: float
    :param denosumab_dose: dose of denosumab in ng/kg body weight
    :type denosumab_dose: float
    """

    def __init__(self):
        self.start_time = 0
        self.end_time = 4000
        self.stress_tensor = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -5]]) * 10 ** -3   # [GPa]
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
        self.end_postmenopausal_osteoporosis = 4000

        # self.start_denosumab_treatment = self.start_postmenopausal_osteoporosis + 365  # leave 1 year untreated
        # self.end_denosumab_treatment = self.start_postmenopausal_osteoporosis + 365 * 2
        self.treatment_period = 183   # every 6 months
        self.start_denosumab_treatment = 365
        self.end_denosumab_treatment = 4000
        self.denosumab_dose = 60 * (10**6)  # [ng]/kg body weight (60 kg as reference body weight) every 6months