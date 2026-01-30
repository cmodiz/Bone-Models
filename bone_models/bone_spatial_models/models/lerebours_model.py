import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from bone_models.bone_spatial_models.parameters.lerebours_parameters import Lerebours_Parameters
from bone_models.bone_cell_population_models.models.lerebours_model import Lerebours_Model as Lerebours_Bone_Cell_Model


class Lerebours_Model:
    """
    Orchestrates the multiscale spatial simulation of bone remodeling.

    This class manages a 2D cross-section of bone (e.g., a femur midshaft)
    represented as a grid of Representative Volume Elements (RVEs). It couples
    macroscopic mechanical loading with local Bone Cell Population Models (BCPM).

    .. note::
       **Source Publication**:
       Lerebours, C., Buenzli, P. R., Scheiner, S., & Pivonka, P. (2016).
       *A multiscale mechanobiological model of bone remodeling predicts
       site-specific bone loss in the femur during osteoporosis and mechanical disuse.*
       Biomechanics and Modeling in Mechanobiology, 15(1), 43-67.
       :doi:`10.1007/s10237-015-0705-x`
    """
    def __init__(self, load_case, duration_of_simulation=3):
        """
        Initializes the spatial model parameters and simulation timing.

        :param load_case: Object defining the mechanical loading scenario.
        :type load_case: Load_Case
        :param duration_of_simulation: Simulation length in years.
        :type duration_of_simulation: int
        """
        self.parameters = Lerebours_Parameters()
        self.load_case = load_case
        self.duration_of_simulation = duration_of_simulation * 365
        self.time_for_mechanics_update = 1 * 365

    def initialize_circular_cross_section(self):
        """
        Generates an idealized circular bone cross-section with radial density zones.

        Creates a grid-based geometry with an outer cortical ring and a
        transitional inner zone using a simple radial mask.

        :return: DataFrame containing coordinates (y, z) and initial BV/TV values.
        :rtype: pandas.DataFrame
        """
        # parameters
        grid_size_mm = 0.8
        num_elements = 40
        length_mm = num_elements * grid_size_mm
        y = np.linspace(-length_mm / 2, length_mm / 2, num_elements)
        z = np.linspace(-length_mm / 2, length_mm / 2, num_elements)
        Y, Z = np.meshgrid(y, z)
        radius = np.sqrt(Y ** 2 + Z ** 2)

        # create circular mask (in mm)
        outer_radius = self.parameters.cross_section.outer_radius
        inner_radius = self.parameters.cross_section.inner_radius
        mask = (radius >= inner_radius) & (radius <= outer_radius)

        # fill in BV/TV values based on radial zones (cortical, transitional, trabecular)
        bv_tv_matrix = np.zeros_like(Y)
        rand_vals_outer = 0.8 + 0.2 * np.random.rand(*Y.shape)
        rand_vals_outer = np.clip(rand_vals_outer, 0.01, 0.99)
        bv_tv_matrix[(radius > 10) & (radius <= 17)] = rand_vals_outer[(radius > 10) & (radius <= 17)]
        rand_vals_mid = 0.3 + 0.1 * np.random.rand(*Y.shape)
        rand_vals_mid = np.clip(rand_vals_mid, 0.01, 0.99)
        bv_tv_matrix[(radius > 7) & (radius <= 10)] = rand_vals_mid[(radius > 7) & (radius <= 10)]
        bv_tv_matrix[~mask] = np.nan

        # plotting
        plt.figure(figsize=(6, 5))
        c = plt.pcolormesh(Z, Y, bv_tv_matrix, cmap='seismic', shading='auto', vmin=0, vmax=1)
        plt.colorbar(c, label='BV/TV')
        plt.title('Femur-like Circular BV/TV Field')
        plt.xlabel('z [mm]')
        plt.ylabel('y [mm]')
        plt.axis('equal')
        plt.tight_layout()
        plt.show()
        # return dataframe
        flat_y = Y.flatten()
        flat_z = Z.flatten()
        flat_bvtv = bv_tv_matrix.flatten()
        valid = ~np.isnan(flat_bvtv)
        df = pd.DataFrame({
            'y': flat_y[valid],
            'z': flat_z[valid],
            'BV/TV': flat_bvtv[valid]
        })
        return df

    def initialize_elliptical_cross_section(self):
        """
        Generates an idealized elliptical femur midshaft cross-section.

        Uses elliptical boundaries to define three distinct zones: dense cortical
        bone, transitional cortex, and the medullary cavity (marrow).

        :return: DataFrame containing coordinates (y, z) and initial BV/TV values.
        :rtype: pandas.DataFrame
        """
        # parameters
        grid_size_mm = 0.8
        num_elements = 40
        length_mm = num_elements * grid_size_mm
        y = np.linspace(-length_mm / 2, length_mm / 2, num_elements)
        z = np.linspace(-length_mm / 2, length_mm / 2, num_elements)
        Y, Z = np.meshgrid(y, z)

        # periosteal parameters (outer boundary)
        peri_y_radius = 17.0  # Mediolateral (wider)
        peri_z_radius = 14.0  # Anteroposterior (narrower)
        # mid-cortical parameters (middle boundary)
        mid_y_radius = 10.0
        mid_z_radius = 7.0
        # endosteal parameters (inner boundary)
        endo_y_radius = 7.0
        endo_z_radius = 4.0

        # create elliptical masks and cortical and transitional zones
        mask_peri = ((Y / peri_y_radius) ** 2 + (Z / peri_z_radius) ** 2) <= 1
        mask_mid = ((Y / mid_y_radius) ** 2 + (Z / mid_z_radius) ** 2) <= 1
        mask_endo = ((Y / endo_y_radius) ** 2 + (Z / endo_z_radius) ** 2) <= 1

        zone_cortical = mask_peri & ~mask_mid
        zone_transitional = mask_mid & ~mask_endo

        # fill in BV/TV values based on zones
        bv_tv_matrix = np.full_like(Y, np.nan)
        rand_vals_cortical = 0.8 + 0.2 * np.random.rand(*Y.shape)
        rand_vals_cortical = np.clip(rand_vals_cortical, 0.01, 0.99)
        rand_vals_trans = 0.3 + 0.1 * np.random.rand(*Y.shape)
        rand_vals_trans = np.clip(rand_vals_trans, 0.01, 0.99)
        bv_tv_matrix[zone_cortical] = rand_vals_cortical[zone_cortical]
        bv_tv_matrix[zone_transitional] = rand_vals_trans[zone_transitional]

        # plotting
        plt.figure(figsize=(6, 5))
        c = plt.pcolormesh(Z, Y, bv_tv_matrix, cmap='seismic', shading='auto', vmin=0, vmax=1)
        plt.colorbar(c, label='BV/TV')
        plt.title('Realistic Femur Midshaft BV/TV Field')
        plt.xlabel('z (A-P) [mm]')
        plt.ylabel('y (M-L) [mm]')
        plt.axis('equal')
        plt.tight_layout()
        plt.show()
        # return dataframe
        flat_y = Y.flatten()
        flat_z = Z.flatten()
        flat_bvtv = bv_tv_matrix.flatten()
        valid = ~np.isnan(flat_bvtv)
        df = pd.DataFrame({
            'y': flat_y[valid],
            'z': flat_z[valid],
            'BV/TV': flat_bvtv[valid]
        })
        print(f"Generated {len(df)} valid RVEs for the midshaft.")
        return df

    def solve_spatial_model(self, only_initialize=False):
        """
        Executes the multiscale simulation loop over the specified duration.

        The solver iterates through time intervals, updating local cell
        populations via BCPMs and periodically recalculating the macroscopic
        stress distribution based on changed bone density/stiffness.

        For initialization, each RVE gets one BCPM instance with mechanical environment depending on the local bone
        volume fraction and the steady-state is calculated. All BCPMs in the cross-section are solved over time intervals,
        after each the macroscopic mechanics ("global" stress tensor) are updated.
        The results are saved at the end of each interval and represent the initial condition for the next interval with
        the "new" macro stress tensor.

        :param only_initialize: If True, returns after setup without running the ODE solver.
        :type only_initialize: bool
        :return: Final cross-section state and a dictionary of temporal results.
        :rtype: tuple(pandas.DataFrame, dict)
        """
        cross_section = self.initialize_elliptical_cross_section()
        models = [Lerebours_Bone_Cell_Model(self.load_case, porosity=1 - bone_volume_fraction) for bone_volume_fraction
                  in cross_section['BV/TV']]
        cross_section['y'] = cross_section['y'] * 1e-3  # Convert to meters
        cross_section['z'] = cross_section['z'] * 1e-3  # Convert to meters
        cross_section['models'] = models
        cross_section = self.update_stress_tensor(cross_section, t=None)

        number_of_intervals = self.duration_of_simulation // self.time_for_mechanics_update
        t_start, t_end = 0, self.time_for_mechanics_update
        results = {t_start: {RVE_index: {} for RVE_index in cross_section.index}}
        # Initialise BCPMs
        for RVE_index, row in cross_section.iterrows():
            bone_volume_fraction, bone_cell_model, stress_xx = (row['BV/TV'], row['models'], row['stress_xx'])
            bone_cell_model.set_macroscopic_stress_tensor(stress_xx * 1e-9, 0, 0, steady_state=True)
            # Calculate steady state
            bone_cell_model.calculate_steady_state(1 - bone_volume_fraction)
            results[t_start][RVE_index]['OBp'] = bone_cell_model.steady_state.OBp
            results[t_start][RVE_index]['OBa'] = bone_cell_model.steady_state.OBa
            results[t_start][RVE_index]['OCp'] = bone_cell_model.steady_state.OCp
            results[t_start][RVE_index]['OCa'] = bone_cell_model.steady_state.OCa
            results[t_start][RVE_index]['porosity'] = 1 - bone_volume_fraction
            results[t_start][RVE_index]['BV/TV'] = bone_volume_fraction
            results[t_start][RVE_index]['SED_bm'] = bone_cell_model.parameters.mechanics.strain_energy_density_steady_state
            results[t_start][RVE_index]['strain_effect_on_OBp'] = bone_cell_model.strain_effect_on_OBp
            results[t_start][RVE_index]['stress_xx'] = stress_xx

        if only_initialize:
            return cross_section, results

        for interval in range(number_of_intervals):
            print(f"--- Interval {interval + 1}: {t_start} to {t_end} days ---")
            results[t_end] = {}
            results[t_end] = {RVE_index: {} for RVE_index in cross_section.index}
            for RVE_index, row in cross_section.iterrows():
                y, z, bone_volume_fraction, element_stiffness, stress_xx, bone_cell_model = (
                    row['y'], row['z'],
                    row['BV/TV'], row['Stiffness'],
                    row['stress_xx'], row['models'])
                # set new stress tensor for the bone cell model
                bone_cell_model.set_macroscopic_stress_tensor(stress_xx * 1e-9, 0, 0)
                # Update the last solution for this RVE (for initial condition for the next interval)
                OBp, OBa, OCp, OCa, porosity, bone_volume_fraction = (results[t_start][RVE_index]['OBp'],
                                                                      results[t_start][RVE_index]['OBa'],
                                                                      results[t_start][RVE_index]['OCp'],
                                                                      results[t_start][RVE_index]['OCa'],
                                                                      results[t_start][RVE_index]['porosity'],
                                                                      results[t_start][RVE_index]['BV/TV'])
                initial_conditions = np.array([OBp, OBa, OCp, OCa, porosity, bone_volume_fraction])
                # bone_cell_model.update_mechanical_effects = True
                bone_cell_model.apply_mechanical_effects(OBp, OBa, OCa, porosity, bone_volume_fraction, t_start)
                # bone_cell_model.update_mechanical_effects = False
                solution = bone_cell_model.solve_bone_cell_population_model(tspan=[t_start, t_end],
                                                                            porosity=1 - bone_volume_fraction,
                                                                            initial_conditions=initial_conditions)
                [OBp, OBa, OCp, OCa, porosity, bone_volume_fraction] = solution.y[:, -1]
                # 2. Save solution at end time point in the results dictionary
                results[t_end][RVE_index]['OBp'] = OBp
                results[t_end][RVE_index]['OBa'] = OBa
                results[t_end][RVE_index]['OCp'] = OCp
                results[t_end][RVE_index]['OCa'] = OCa
                results[t_end][RVE_index]['porosity'] = porosity
                results[t_end][RVE_index]['BV/TV'] = bone_volume_fraction
                results[t_end][RVE_index]['SED_bm'] = bone_cell_model.strain_energy_density
                results[t_end][RVE_index]['strain_effect_on_OBp'] = bone_cell_model.strain_effect_on_OBp

                cross_section.at[RVE_index, 'BV/TV'] = solution.y[5][-1]
            cross_section = self.update_stress_tensor(cross_section, t_end)
            t_start = t_end
            t_end += self.time_for_mechanics_update
        return cross_section, results

    def update_stress_tensor(self, cross_section, t):
        """
        Recalculates the local axial stress for every RVE based on current stiffness.

        Uses beam theory to decompose global forces and moments into local
        strains, then applies the constitutive law to find local stress.

        .. note::
           Mapped to internal z-direction (stress_xx index [2,2]) to align
           the longitudinal axis across different coordinate conventions (Scheiner 2013 v. Lerebours 2016).

        :param cross_section: Current spatial data of the bone.
        :type cross_section: pandas.DataFrame
        :param t: Current simulation time.
        :type t: float or None
        :return: Updated cross-section with 'stress_xx' values.
        :rtype: pandas.DataFrame
        """
        cross_section = self.calculate_stiffness_for_all_RVEs(cross_section)
        axial_strain, curvature_y, curvature_z, y_normal_force_center, z_normal_force_center = self.calculate_strain_decomposition(
            cross_section, t)
        cross_section['stress_xx'] = np.nan
        for index, row in cross_section.iterrows():
            y, z, bone_volume_fraction, element_stiffness = row['y'], row['z'], row['BV/TV'], row['Stiffness']
            strain_xx = axial_strain - curvature_y * (y - y_normal_force_center) + curvature_z * (
                    z - z_normal_force_center)
            stress_xx = element_stiffness * strain_xx
            cross_section.at[index, 'stress_xx'] = stress_xx
        return cross_section

    def calculate_strain_decomposition(self, cross_section, t):
        """
        Decomposes global forces into axial strain and curvatures.

        Solves the system of equations linking axial force and bending moments
        to the geometric stiffness matrix of the cross-section: determine normal force center, axial stiffness and
        second moments; insert in moment of area matrix and invert to calculate strain decomposition.

        :param cross_section: Current spatial data.
        :type cross_section: pandas.DataFrame
        :param t: Current simulation time.
        :type t: float
        :return: Tuple of (axial_strain, curvature_y, curvature_z, y_nf_center, z_nf_center).
        :rtype: tuple
        """
        y_normal_force_center, z_normal_force_center, axial_stiffness = self.calculate_normal_force_center(
            cross_section)
        second_moment_y, second_moment_z, second_moment_yz = self.calculate_second_moments_of_area(cross_section,
                                                                                                   y_normal_force_center,
                                                                                                   z_normal_force_center)
        if t is None or t <= self.load_case.start_time or t >= self.load_case.end_time:
            axial_force = self.parameters.mechanics.axial_force
            bending_moment_y = self.parameters.mechanics.bending_moment_y
            bending_moment_z = self.parameters.mechanics.bending_moment_z
        else:
            axial_force = self.load_case.force_reduction * self.parameters.mechanics.axial_force
            bending_moment_y = self.load_case.moment_reduction * self.parameters.mechanics.bending_moment_y
            bending_moment_z = self.load_case.moment_reduction * self.parameters.mechanics.bending_moment_z
        force_and_moments = np.array([axial_force, bending_moment_y, bending_moment_z])
        moments_of_area_matrix = np.array([[axial_stiffness, 0, 0],
                                           [0, second_moment_y, -second_moment_yz],
                                           [0, -second_moment_yz, second_moment_z]])
        strain_decomposition = np.linalg.inv(moments_of_area_matrix) @ force_and_moments
        axial_strain, curvature_z, curvature_y = strain_decomposition[0], strain_decomposition[1], strain_decomposition[
            2]
        return axial_strain, curvature_y, curvature_z, y_normal_force_center, z_normal_force_center

    def calculate_normal_force_center(self, cross_section):
        """
        Calculates the stiffness-weighted centroid (neutral axis) of the bone.

        :param cross_section: Current spatial data (coordinates and stiffness for each RVE).
        :type cross_section: pandas.DataFrame
        :return: Tuple of (y_nf_center, z_nf_center, total_axial_stiffness).
        :rtype: tuple
        """
        y_normal_force_center = np.sum(cross_section['y'] * cross_section['Stiffness']) / np.sum(
            cross_section['Stiffness'])
        z_normal_force_center = np.sum(cross_section['z'] * cross_section['Stiffness']) / np.sum(
            cross_section['Stiffness'])
        axial_stiffness = np.sum(
            cross_section['Stiffness']) * self.parameters.cross_section.delta_y * self.parameters.cross_section.delta_z
        return y_normal_force_center, z_normal_force_center, axial_stiffness

    def calculate_second_moments_of_area(self, cross_section, y_normal_force_center, z_normal_force_center):
        """
        Computes the stiffness-weighted second moments and product of area.

        :param cross_section: Current spatial data (coordinates and stiffness for each RVE).
        :type cross_section: pandas.DataFrame
        :param y_normal_force_center: Y-coordinate of the neutral axis.
        :type y_normal_force_center: float
        :param z_normal_force_center: Z-coordinate of the neutral axis.
        :type z_normal_force_center: float
        :return: Tuple of (Iyy, Izz, Iyz).
        :rtype: tuple
        """
        second_moment_y = np.sum(
            cross_section['Stiffness'] * ((cross_section['z'] - z_normal_force_center) ** 2) * self.parameters.cross_section.delta_y * self.parameters.cross_section.delta_z)
        second_moment_z = np.sum(
            cross_section['Stiffness'] * ((cross_section['y'] - y_normal_force_center) ** 2) * self.parameters.cross_section.delta_y * self.parameters.cross_section.delta_z)
        second_moment_yz = np.sum(
            cross_section['Stiffness'] * (cross_section['y'] - y_normal_force_center) * (
                    cross_section['z'] - z_normal_force_center) *
            self.parameters.cross_section.delta_y * self.parameters.cross_section.delta_z)
        return second_moment_y, second_moment_z, second_moment_yz

    def calculate_stiffness_for_all_RVEs(self, cross_section):
        """
        Performs homogenization to find the effective longitudinal stiffness of each RVE.

        Calls the BCPM homogenization methods to determine the macroscopic
        stiffness tensor based on local bone volume fraction.

        :param cross_section: Current spatial data.
        :type cross_section: pandas.DataFrame
        :return: Updated cross-section with 'Stiffness' column (in Pa).
        :rtype: pandas.DataFrame
        """
        cross_section['Stiffness'] = np.nan
        for index, row in cross_section.iterrows():
            bone_volume_fraction, bone_cell_model = row['BV/TV'], row['models']
            strain_concentration_tensor_bone_matrix, strain_concentration_tensor_vascular_pores = bone_cell_model.calculate_strain_concentration_tensors(
                bone_volume_fraction*100)
            macroscopic_stiffness_tensor = bone_cell_model.calculate_macroscopic_stiffness_tensor(
                strain_concentration_tensor_bone_matrix,
                strain_concentration_tensor_vascular_pores,
                (1 - bone_volume_fraction)*100,
                bone_volume_fraction*100)
            cross_section.at[index, 'Stiffness'] = macroscopic_stiffness_tensor[2, 2] * 1e9  # Convert to Pa
        return cross_section
