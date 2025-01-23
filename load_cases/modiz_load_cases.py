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


class Healthy_to_Glucocorticoid_Induced_Osteoporosis:
    def __init__(self):
        self.lemaire = Lemaire_Load_Case()
        self.martonova = Glucocorticoid_Induced_Osteoporosis()



