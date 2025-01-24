from load_cases.martonova_load_cases import Healthy, Hyperparathyroidism, Osteoporosis, Postmenopausal_Osteoporosis, \
    Hypercalcemia, Hypocalcemia, Glucocorticoid_Induced_Osteoporosis


class Lemaire_Load_Case:
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
        self.end_time = 80


class Healthy_to_Hyperparathyroidism:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Hyperparathyroidism()


class Healthy_to_Osteoporosis:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Osteoporosis()


class Healthy_to_Postmenopausal_Osteoporosis:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Postmenopausal_Osteoporosis()


class Healthy_to_Hypercalcemia:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Hypercalcemia()


class Healthy_to_Hypocalcemia:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Hypocalcemia()


class Healthy_to_Glucocorticoid_Induced_Osteoporosis:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Glucocorticoid_Induced_Osteoporosis()


# Load cases for reference Lemaire model
class Reference_Healthy_to_Hyperparathyroidism:
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
        # 1.0, 3.8782894736842097, 0.8252796052631578, 0.9, 0.19098684210526315, 7.3500657894736845, 1.0565131578947367
        self.start_time = 20
        self.end_time = 80


class Reference_Healthy_to_Osteoporosis:
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
        self.end_time = 80


class Reference_Healthy_to_Postmenopausal_Osteoporosis:
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
        self.end_time = 80


class Reference_Healthy_to_Hypercalcemia:
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
        self.end_time = 80


class Reference_Healthy_to_Hypocalcemia:
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
        self.end_time = 80


class Reference_Healthy_to_Glucocorticoid_Induced_Osteoporosis:
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
        self.end_time = 80


