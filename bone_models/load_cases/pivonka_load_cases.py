import numpy as np


class Pivonka_Load_Case_1:
    """ Load case for Pivonka et al. 2008 model: injection of PTH during a specified time period. The same load cases as
    in Lemaire et al. 2004 are implemented.

    :param OBp_injection: Osteoblast precursor injection rate
    :type OBp_injection: float
    :param OBa_injection: Osteoblast injection rate
    :type OBa_injection: float
    :param OCa_injection: Osteoclast injection rate
    :type OCa_injection: float
    :param PTH_injection: Parathyroid hormone injection rate
    :type PTH_injection: float
    :param OPG_injection: Osteoprotegerin injection rate
    :type OPG_injection: float
    :param RANKL_injection: Receptor activator of nuclear factor kappa-B ligand injection rate
    :type RANKL_injection: float
    :param TGFb_injection: Transforming growth factor beta injection rate
    :type TGFb_injection: float
    :param differentiation_rate_OCp_multiplier: Osteoclast precursor differentiation rate multiplier
    :type differentiation_rate_OCp_multiplier: float
    :param start_time: Start time of the load case
    :type start_time: float
    :param end_time: End time of the load case
    :type end_time: float
    """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 1.0e+3
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.TGFb_injection = 0
        self.differentiation_rate_OCp_multiplier = np.exp(1)
        self.start_time = 0
        self.end_time = 100