from .martonova_load_cases import Martonova_Healthy, Martonova_Hyperparathyroidism, Martonova_Osteoporosis, Martonova_Postmenopausal_Osteoporosis, \
    Martonova_Hypercalcemia, Martonova_Hypocalcemia, Martonova_Glucocorticoid_Induced_Osteoporosis


class Lemaire_Load_Case:
    """ Load case for the reference class Lemaire model with external injections. """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.start_time = 20
        self.end_time = 385


class Modiz_Healthy_to_Hyperparathyroidism:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Hyperparathyroidism.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Hyperparathyroidism load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Hyperparathyroidism"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Hyperparathyroidism()


class Modiz_Healthy_to_Osteoporosis:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Osteoporosis.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Osteoporosis load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Osteoporosis"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Osteoporosis()


class Modiz_Healthy_to_Postmenopausal_Osteoporosis:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Postmenopausal Osteoporosis.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Postmenopausal Osteoporosis load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Postmenopausal_Osteoporosis"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Postmenopausal_Osteoporosis()


class Modiz_Healthy_to_Hypercalcemia:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Hypercalcemia.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Hypercalcemia load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Hypercalcemia"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Hypercalcemia()


class Modiz_Healthy_to_Hypocalcemia:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Hypocalcemia.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Hypocalcemia load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Hypocalcemia"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Hypocalcemia()


class Modiz_Healthy_to_Glucocorticoid_Induced_Osteoporosis:
    """ Load case for the subclass Modiz_Model for the transition from Healthy to Glucocorticoid-Induced Osteoporosis.
    It contains the load case class for the Lemaire model and the Martonova model.
    The latter is directly imported from the Martonova load cases.

    :param lemaire: Lemaire load case for reference model
    :type lemaire: Lemaire_Load_Case
    :param martonova: Glucocorticoid-Induced Osteoporosis load case for Martonova model
    :type martonova: bone_models.load_cases.martonova_load_cases.Martonova_Glucocorticoid_Induced_Osteoporosis"""
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Martonova_Glucocorticoid_Induced_Osteoporosis()


class Modiz_Reference_Healthy_to_Hyperparathyroidism:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Hyperparathyroidism.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 3.8782894736842097
        self.start_time = 20
        self.end_time = 385


class Modiz_Reference_Healthy_to_Osteoporosis:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Osteoporosis.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 0.8252796052631578
        self.start_time = 20
        self.end_time = 385


class Modiz_Reference_Healthy_to_Postmenopausal_Osteoporosis:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Postmenopausal Osteoporosis.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 0.9
        self.start_time = 20
        self.end_time = 385


class Modiz_Reference_Healthy_to_Hypercalcemia:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Hypercalcemia.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 0.19098684210526315
        self.start_time = 20
        self.end_time = 385


class Modiz_Reference_Healthy_to_Hypocalcemia:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Hypocalcemia.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 7.3500657894736845
        self.start_time = 20
        self.end_time = 385


class Modiz_Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis:
    """ Load case for the subclass Reference_Lemaire_Model for the transition from Healthy to Glucocorticoid-Induced Osteoporosis.
    It contains a PTH_elevation factor, which is used to align disease states with the Modiz_Model to make them comparable.
    The PTH_elevation factor is calculated by the ratio of the maximum PTH values in :func:`bone_models.models.modiz_model.calculate_elevation_parameter`.

    :param PTH_elevation: PTH_elevation factor to align the disease states
    :type PTH_elevation: float """
    def __init__(self):
        self.OBp_injection = 0
        self.OBa_injection = 0
        self.OCa_injection = 0
        # -> I_P
        self.PTH_injection = 0
        # -> I_O
        self.OPG_injection = 0
        # -> I_L
        self.RANKL_injection = 0
        self.PTH_elevation = 1.0565131578947367
        self.start_time = 20
        self.end_time = 385


