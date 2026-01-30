class mechanics:
    """ This class defines the mechanical loading parameters of the bone model.

    The following table summarizes the parameters and matches them to the variables used in the publication.

    +-------------------+----------------------+
    | Parameter         |Variable name         |
    +===================+======================+
    | axial_force       | :math:`N`            |
    +-------------------+----------------------+
    | bending_moment_y  | :math:`M_y`          |
    +-------------------+----------------------+
    | bending_moment_z  | :math:`M_z`          |
    +-------------------+----------------------+

    :param axial_force: axial force applied to the bone model (N)
    :type axial_force: float
    :param bending_moment_y: bending moment applied to the bone model around the y-axis
    :type bending_moment_y: float
    :param bending_moment_z: bending moment applied to the bone model around the z-axis
    :type bending_moment_z: float """
    def __init__(self):
        self.axial_force = -700  # N
        self.bending_moment_y = 50  # Nm
        # # self.bending_moment_z = 50  # Nm
        self.bending_moment_z = 0  # Nm
        # body_weight = 70  # kg
        # force_body_weight = body_weight * 9.81  # N
        # self.axial_force = - 1.12 * force_body_weight  # N
        # self.bending_moment_y = - 0.05 * force_body_weight  # Nm
        # # self.bending_moment_z = 50  # Nm
        # self.bending_moment_z = 0.02 * force_body_weight  # Nm

class cross_section:
    """ This class defines the cross section geometry parameters of the bone model.

    The following table summarizes the parameters and matches them to the variables used in the publication.

    +-------------------+----------------------+
    | Parameter         | Variable name        |
    +===================+======================+
    | delta_y           | :math:`\Delta_y`     |
    +-------------------+----------------------+
    | delta_z           | :math:`\Delta_z`     |
    +-------------------+----------------------+
    | outer_radius      | -                    |
    +-------------------+----------------------+
    | inner_radius      | -                    |
    +-------------------+----------------------+

    :param delta_y: height of the cross section in y-direction (m)
    :type delta_y: float
    :param delta_z: width of the cross section in z-direction (m)
    :type delta_z: float
    :param outer_radius: outer radius of the bone cross section (cortical) (mm)
    :type outer_radius: float
    :param inner_radius: inner radius of the bone cross section (trabecular) (mm)
    :type inner_radius: float
    """
    def __init__(self):
        self.delta_y = 0.8 * 1e-3  # mm in m
        self.delta_z = 0.8 * 1e-3  # mm in m
        self.outer_radius = 17  # mm
        self.inner_radius = 7

class Lerebours_Parameters:
    """ This class defines the parameters of the bone model.

    :param mechanics: parameters relevant for mechanics of the bone model
    :type mechanics: mechanics
    :param cross_section: parameters relevant for the cross section geometry of the bone model
    :type cross_section: cross_section"""

    def __init__(self):
        self.mechanics = mechanics()
        self.cross_section = cross_section()
