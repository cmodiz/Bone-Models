class Lemaire_Load_Case_1:
    """ Load case 1 for the Lemaire model: injection of active osteoblasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 1.0e-4
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_2:
    """ Load case 2 for the Lemaire model: retraction of active osteoblasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = -8.3e-5
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_3:
    """ Load case 3 for the Lemaire model: injection of parathyroid hormone in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 1.0e+3
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_4:
    """ Load case 4 for the Lemaire model: injection of active osteoclasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 1.0e-4
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_5:
    """ Load case 5 for the Lemaire model: retraction of active osteoclasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = -2.9e-4
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_6:
    """ Load case 6 for the Lemaire model: injection of osteoprotegerin in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 2.0e+5
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_7:
    """ Load case 7 for the Lemaire model: injection of precursor osteoblasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 1.0e-4
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_8:
    """ Load case 8 for the Lemaire model: retraction of precursor osteoblasts in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = -1.2e-3
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 80


class Lemaire_Load_Case_9:
    """ Load case 9 for the Lemaire model: injection of receptor activator of nuclear factor kappa-B ligand in a specific time interval.

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
    :param start_time: start time of the external injections
    :type start_time: float
    :param end_time: end time of the external injections
    :type end_time: float
    """
    def __init__(self):
        """ Constructor method """
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 10
        self.start_time = 20
        self.end_time = 80
