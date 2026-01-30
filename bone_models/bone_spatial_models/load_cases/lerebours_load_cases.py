class Lerebours_Load_Case_Osteoporosis:
    """ Load case representing osteoporosis conditions (altered PTH levels) based on the publication.
    Relevant to this specific model are the parameters PTH_injection, force_reduction, and moment_reduction.
    The latter two are zero in this case, as no mechanical unloading is considered. All other parameters are remainders
    from the base class.

    :param start_time: Start time of the load case [days]
    :type start_time: int
    :param end_time: End time of the load case [days]
    :type end_time: int
    :param PTH_injection: Injection rate of Parathyroid Hormone (PTH)
    :type PTH_injection: float
    :param force_reduction: Reduction factor multiplied with applied forces (1 - reduction_fraction)
    :type force_reduction: float
    :param moment_reduction: Reduction factor multiplied with applied moments (1 - reduction_fraction)
    :type moment_reduction: float
    """
    def __init__(self):
        self.start_time = 0
        self.end_time = 500000000000000
        self.stress_tensor = None   # [GPa]
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        self.PTH_injection = 0.047
        self.OPG_injection = 0
        self.RANKL_injection = 0
        self.TGFb_injection = 0
        self.force_reduction = 1 - 0  # 0% reduction in force
        self.moment_reduction = 1 - 0  # 0% reduction in moment


class Lerebours_Load_Case_Spaceflight:
    """ Load case representing microgravity conditions (reduced loading) based on the publication.
    Relevant to this specific model are the parameters PTH_injection, force_reduction, and moment_reduction.
    The first one is zero in this case. All other parameters are remainders
    from the base class.

    :param start_time: Start time of the load case [days]
    :type start_time: int
    :param end_time: End time of the load case [days]
    :type end_time: int
    :param PTH_injection: Injection rate of Parathyroid Hormone (PTH)
    :type PTH_injection: float
    :param force_reduction: Reduction factor multiplied with applied forces (1 - reduction_fraction)
    :type force_reduction: float
    :param moment_reduction: Reduction factor multiplied with applied moments (1 - reduction_fraction)
    :type moment_reduction: float
    """
    def __init__(self):
        self.start_time = 0
        self.end_time = 500000000000000
        self.stress_tensor = None   # [GPa]
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
        self.force_reduction = 1 - 0.8  # 80% reduction in force
        self.moment_reduction = 1 - 0.8  # 80% reduction in force
