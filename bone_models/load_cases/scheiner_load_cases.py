import numpy as np


class Scheiner_Load_Case:
    """ Load cases of the Scheiner et al., 2013 model: disuse simulation in the time interval [50, 2050] days using a
    stress tensor with reduced unixial compressive stress (-25 vs. -30 MPa in habitual loading according to paper).

    :param start_time: start time of the load case
    :type start_time: float
    :param end_time: end time of the load case
    :type end_time: float
    :param stress_tensor: applied stress tensor
    :type stress_tensor: numpy.ndarray
    :param OBp_injection: injection of precursor osteoblasts
    :type OBp_injection: float
    :param OBa_injection: injection of active osteoblasts
    :type OBa_injection: float
    :param OCa_injection: injection of active osteoclasts
    :type OCa_injection: float
    :param PTH_injection: injection of parathyroid hormone
    :type PTH_injection: float
    :param OPG_injection: injection of osteoprotegerin
    :type OPG_injection: float
    :param RANKL_injection: injection of receptor activator of nuclear factor kappa-B ligand
    :type RANKL_injection: float
    :param TGFb_injection: injection of transforming growth factor beta
    :type TGFb_injection: float"""
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
