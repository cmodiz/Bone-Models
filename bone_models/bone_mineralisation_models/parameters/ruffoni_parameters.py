class Calcium:
    def __init__(self):
        self.minimum_content = 0  # [%wt]
        self.maximum_content = 31  # %


class Reference_BMDD:
    def __init__(self):
        self.standard_deviation_lower_percentile = 4.93 / 1.5
        self.standard_deviation_upper_percentile = 5.55 / 1.5
        self.peak = 22.94


class Rate:
    def __init__(self):
        self.initial_formation = 0.1  # [mm^3/ year]
        self.final_formation = 0.1  # [mm^3/ year]
        self.initial_resorption = 0.1  # [1/ year]
        self.final_resorption = 0.1  # [1/ year]


class Mineralization_Law:
    def __init__(self):
        self.turnover_time = 5  # [years]
        self.primary_mineral_content = 11.9  # [%wt]
        self.maximum_mineral_content = 31.1 - self.primary_mineral_content  # [%wt]
        self.primary_apposition_rate = 150  # [years]
        self.secondary_apposition_rate = 0.1  # [years]


class Bone_Volume:
    def __init__(self):
        # needed when initializing the model with a mineralization law
        self.initial_value = 2.4  # [mm^3]


class Ruffoni_Parameters:
    def __init__(self):
        self.calcium = Calcium()
        self.rate = Rate()
        self.reference_bmdd = Reference_BMDD()
        self.mineralization_law = Mineralization_Law()
        self.bone_volume = Bone_Volume()
