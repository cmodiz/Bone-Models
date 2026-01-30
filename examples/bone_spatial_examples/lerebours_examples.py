from bone_models.bone_spatial_models.models.lerebours_model import Lerebours_Model
from bone_models.bone_spatial_models.load_cases.lerebours_load_cases import Lerebours_Load_Case_Spaceflight
import matplotlib.pyplot as plt
import numpy as np


def plot_initial_sed_distribution(cross_section, results):
    """
    Plots the initial Strain Energy Density (SED_bm) distribution.

    :param cross_section: DataFrame containing spatial coordinates.
    :param results: Dictionary containing the solution for each interval and RVE.
    """
    print("--- Generating Initial SED Distribution Plot ---")
    y_coords = cross_section['y'].values * 1e3  # Convert to mm
    z_coords = cross_section['z'].values * 1e3  # Convert to mm

    # Get results for the initial interval (t=0)
    try:
        initial_results = results[0]
    except KeyError:
        print("Error: results[0] not found. Cannot plot initial SED.")
        return

    # Create a meshgrid for plotting
    unique_y = np.unique(y_coords)
    unique_z = np.unique(z_coords)
    Y, Z = np.meshgrid(unique_y, unique_z)

    # Reshape SED values to match the grid
    SED_grid = np.full_like(Y, np.nan, dtype=float)

    # Check if 'SED_bm' is available (it should be, but good to check)
    if not any('SED_bm' in rve_data for rve_data in initial_results.values()):
        print("\n*** Error: 'SED_bm' not found in initial results (results[0]). ***")
        print("Skipping initial SED distribution plot.")
        return

    for RVE_index, solution in initial_results.items():
        sed = solution['SED_bm']  # Extract the SED value for this RVE
        y, z = cross_section.loc[RVE_index, ['y', 'z']].values * 1e3
        y_idx = np.where(unique_y == y)[0][0]
        z_idx = np.where(unique_z == z)[0][0]
        SED_grid[z_idx, y_idx] = sed

    # Plot the SED distribution for this interval
    plt.figure(figsize=(6, 5))
    # Use 'viridis' or 'hot' colormap for energy (usually positive)
    c = plt.pcolormesh(Z, Y, SED_grid, cmap='viridis', shading='auto')
    plt.colorbar(c, label='SED_bm')
    plt.title(f'Initial SED (bm) Distribution')
    plt.xlabel('z [mm]')
    plt.ylabel('y [mm]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Plots/straight_beam_uniaxial/initial_sed_distribution.pdf', dpi=500)


def plot_bvtv_for_intervals(cross_section, results):
    """
    Plots the BV/TV distribution for each interval from the results.

    :param cross_section: DataFrame containing spatial coordinates and initial BV/TV values.
    :param results: Dictionary containing the solution for each interval and RVE.
    """
    y_coords = cross_section['y'].values * 1e3  # Convert to mm
    z_coords = cross_section['z'].values * 1e3  # Convert to mm

    for interval, interval_results in results.items():
        # Create a meshgrid for plotting
        unique_y = np.unique(y_coords)
        unique_z = np.unique(z_coords)
        Y, Z = np.meshgrid(unique_y, unique_z)

        # Reshape BV/TV values to match the grid
        BV_TV_grid = np.full_like(Y, np.nan, dtype=float)
        for RVE_index, solution in interval_results.items():
            bvtv = solution['BV/TV']  # Extract the BV/TV value for this RVE
            y, z = cross_section.loc[RVE_index, ['y', 'z']].values * 1e3
            y_idx = np.where(unique_y == y)[0][0]
            z_idx = np.where(unique_z == z)[0][0]
            BV_TV_grid[z_idx, y_idx] = bvtv

        # Plot the BV/TV distribution for this interval
        plt.figure(figsize=(6, 5))
        c = plt.pcolormesh(Z, Y, BV_TV_grid, cmap='seismic', shading='auto', vmin=0, vmax=1)
        plt.colorbar(c, label='BV/TV')
        plt.title(f'BV/TV Distribution - Interval {interval + 1}')
        plt.xlabel('z [mm]')
        plt.ylabel('y [mm]')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(f'bvtv_distribution_interval_{interval + 1}.pdf')


def plot_initial_stress_distribution(cross_section, results):
    """
    Plots the initial stress (zz) distribution.

    :param cross_section: DataFrame containing spatial coordinates.
    :param results: Dictionary containing the solution for each interval and RVE.
    """
    print("--- Generating Initial Stress Distribution Plot ---")
    y_coords = cross_section['y'].values * 1e3  # Convert to mm
    z_coords = cross_section['z'].values * 1e3  # Convert to mm

    # Get results for the initial interval (t=0)
    try:
        initial_results = results[0]
    except KeyError:
        print("Error: results[0] not found. Cannot plot initial stress.")
        return

    # Create a meshgrid for plotting
    unique_y = np.unique(y_coords)
    unique_z = np.unique(z_coords)
    Y, Z = np.meshgrid(unique_y, unique_z)

    # Reshape Stress values to match the grid
    Stress_grid = np.full_like(Y, np.nan, dtype=float)

    # Check if 'stress_xx' is available
    if not any('stress_xx' in rve_data for rve_data in initial_results.values()):
        print("\n*** Error: 'stress_xx' not found in initial results (results[0]). ***")
        print("Skipping initial stress distribution plot.")
        return

    for RVE_index, solution in initial_results.items():
        stress = solution['stress_xx']  # Extract the stress value for this RVE
        y, z = cross_section.loc[RVE_index, ['y', 'z']].values * 1e3
        y_idx = np.where(unique_y == y)[0][0]
        z_idx = np.where(unique_z == z)[0][0]
        Stress_grid[z_idx, y_idx] = stress

    # Plot the Stress distribution for this interval
    plt.figure(figsize=(6, 5))
    # Use 'seismic' colormap, remove vmin/vmax
    c = plt.pcolormesh(Z, Y, Stress_grid, cmap='seismic', shading='auto')
    plt.colorbar(c, label='Stress xx(Pa)')
    plt.title(f'Initial Stress (xx) Distribution')
    plt.xlabel('z [mm]')
    plt.ylabel('y [mm]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Plots/straight_beam_uniaxial/initial_stress_xx_distribution.pdf', dpi=500)


def plot_initial_strain_effect_distribution(cross_section, results):
    """
    Plots the initial Strain Effect on OBp (strain_effect_on_OBp) distribution.

    :param cross_section: DataFrame containing spatial coordinates.
    :param results: Dictionary containing the solution for each interval and RVE.
    """
    print("--- Generating Initial Strain Effect on OBp Distribution Plot ---")
    y_coords = cross_section['y'].values * 1e3  # Convert to mm
    z_coords = cross_section['z'].values * 1e3  # Convert to mm

    # Get results for the initial interval (t=0)
    try:
        initial_results = results[0]
    except KeyError:
        print("Error: results[0] not found. Cannot plot initial Strain Effect on OBp.")
        return

    # Create a meshgrid for plotting
    unique_y = np.unique(y_coords)
    unique_z = np.unique(z_coords)
    Y, Z = np.meshgrid(unique_y, unique_z)

    # Reshape Strain Effect values to match the grid
    StrainEffect_grid = np.full_like(Y, np.nan, dtype=float)

    # Check if 'strain_effect_on_OBp' is available
    if not any('strain_effect_on_OBp' in rve_data for rve_data in initial_results.values()):
        print("\n*** Error: 'strain_effect_on_OBp' not found in initial results (results[0]). ***")
        print("Skipping initial Strain Effect on OBp distribution plot.")
        return

    for RVE_index, solution in initial_results.items():
        strain_effect = solution['strain_effect_on_OBp']  # Extract the value for this RVE
        y, z = cross_section.loc[RVE_index, ['y', 'z']].values * 1e3
        y_idx = np.where(unique_y == y)[0][0]
        z_idx = np.where(unique_z == z)[0][0]
        StrainEffect_grid[z_idx, y_idx] = strain_effect

    # Plot the Strain Effect distribution
    plt.figure(figsize=(6, 5))
    # Using 'plasma' colormap which works well for positive, high-to-low values
    c = plt.pcolormesh(Z, Y, StrainEffect_grid, cmap='plasma', shading='auto')
    plt.colorbar(c, label='Strain Effect on $OB_p$')
    plt.title(f'Initial Strain Effect on Osteoblasts ($OB_p$) Distribution')
    plt.xlabel('z [mm]')
    plt.ylabel('y [mm]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Plots/straight_beam_uniaxial/initial_strain_effect_on_OBp_distribution.pdf', dpi=500)


def plot_profiles(cross_section, results, model):
    # Helper to extract profiles along y=0 and z=0
    def extract_profile(df, value_col, axis='z'):
        # if axis == 'z':
        #     profile = df[np.isclose(df['y'], 0, rtol=1e-02, atol=1e-03)]
        #     x = profile['z'].values * 1e3
        # else:
        #     profile = df[np.isclose(df['z'], 0, rtol=1e-02, atol=1e-03)]
        #     x = profile['y'].values * 1e3
        if axis == 'z':
            # get 40 rows where y is closest to 0, keep original order
            idx = (df['y'] - 0).abs().argsort()[:40]
            profile = df.iloc[idx]
            # sort by z (the plotting axis)
            profile = profile.sort_values(by='z')
            x = profile['z'].values * 1e3
        else:
            # get 40 rows where z is closest to 0, keep original order
            idx = (df['z'] - 0).abs().argsort()[:40]
            profile = df.iloc[idx]
            # sort by y (the plotting axis)
            profile = profile.sort_values(by='y')
            x = profile['y'].values * 1e3
        y = profile[value_col].values
        return x, y

    # Prepare DataFrames for start and end
    cross_section['BV/TV_start'] = [results[0][i]['BV/TV'] for i in range(len(cross_section))]
    cross_section['BV/TV_end'] = [results[max(results.keys())][i]['BV/TV'] for i in range(len(cross_section))]
    cross_section['BV/TV_diff'] = cross_section['BV/TV_start'] - cross_section['BV/TV_end']

    # Prepare DataFrames for start and end
    cross_section['SED_bm_start'] = [results[0][i]['SED_bm'] for i in range(len(cross_section))]
    cross_section['SED_bm_end'] = [results[max(results.keys())][i]['SED_bm'] for i in range(len(cross_section))]
    cross_section['SED_bm_diff'] = cross_section['SED_bm_start'] - cross_section['SED_bm_end']

    # Prepare DataFrames for start and end
    cross_section['strain_effect_on_OBp_start'] = [results[0][i]['strain_effect_on_OBp'] for i in
                                                   range(len(cross_section))]
    cross_section['strain_effect_on_OBp_end'] = [results[max(results.keys())][i]['strain_effect_on_OBp'] for i in
                                                 range(len(cross_section))]
    cross_section['strain_effect_on_OBp_diff'] = cross_section['strain_effect_on_OBp_start'] - cross_section[
        'strain_effect_on_OBp_end']

    # Repeat for stress_xx, strain energy density, and strain effect on OBp
    for key, label in [
        ('BV/TV', 'Bone Volume Fraction'),
        ('SED_bm', 'Strain Energy Density'),
        ('strain_effect_on_OBp', 'Strain Effect on OBp')
    ]:
        if key == 'BV/TV':
            for axis in ['z', 'y']:
                plt.figure()

                # loop over start/end/diff
                for suffix, title, style in zip(
                        ['start', 'end', 'diff'],
                        ['Start', 'End', 'Difference'],
                        ['-o', '-s', '--^']  # different line/marker styles
                ):
                    x, y = extract_profile(cross_section, f'{key}_{suffix}', axis)
                    plt.plot(x, y, style, label=f'{title}')

                plt.title(f'{label} Profile (axis={axis})')
                plt.xlabel('z [mm]' if axis == 'z' else 'y [mm]')
                plt.ylabel(label)
                plt.legend()
                plt.savefig(f'Plots/straight_beam_uniaxial/BV_TV_profile_axis_{axis}.pdf', dpi=500)
        else:
            for axis in ['z', 'y']:
                plt.figure()

                # loop over start/end/diff
                for suffix, title, style in zip(
                        ['start', 'end', 'diff'],
                        ['Start', 'End', 'Difference'],
                        ['-o', '-s', '--^']  # different line/marker styles
                ):
                    x, y = extract_profile(cross_section, f'{key}_{suffix}', axis)
                    plt.plot(x, y, style, label=f'{title}')

                plt.title(f'{label} Profile (axis={axis})')
                plt.xlabel('z [mm]' if axis == 'z' else 'y [mm]')
                plt.ylabel(label)
                plt.legend()
                plt.savefig(f'Plots/straight_beam_uniaxial/{key}_profile_axis_{axis}.pdf', dpi=500)
pass


def plot_initial_stress_sed_vs_bvtv(results):
    """
    Plots the initial stress and SED vs. bone volume fraction as scatter plots.

    :param results: Dictionary containing the solution for each interval and RVE.
    """
    print("--- Generating Initial State Scatter Plots ---")

    # Extract initial data (t=0)
    initial_results = results[0]

    bvtv_values = []
    sed_values = []
    stress_values = []

    # Check if 'stress_xx' is available before proceeding
    if not any('stress_xx' in rve_data for rve_data in initial_results.values()):
        print("\n*** Error: 'stress_xx' not found in initial results (results[0]). ***")
        print("Skipping initial stress/SED scatter plots.")
        return

    for rve_data in initial_results.values():
        bvtv_values.append(rve_data['BV/TV'])
        sed_values.append(rve_data['SED_bm'])
        stress_values.append(rve_data['stress_xx'])

    # --- Plot 1: SED vs BV/TV ---
    plt.figure(figsize=(8, 6))
    plt.scatter(bvtv_values, sed_values, alpha=0.7, s=15, c='blue', edgecolors='k', linewidth=0.5)
    plt.title('Initial Strain Energy Density (SED) vs. Bone Volume Fraction', fontsize=16)
    plt.xlabel('Bone Volume Fraction (BV/TV)', fontsize=12)
    plt.ylabel('Strain Energy Density (SED_bm)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('Plots/straight_beam_uniaxial/initial_sed_vs_bvtv.pdf', dpi=500)

    # --- Plot 2: Stress vs BV/TV ---
    plt.figure(figsize=(8, 6))
    plt.scatter(bvtv_values, stress_values, alpha=0.7, s=15, c='red', edgecolors='k', linewidth=0.5)
    plt.title('Initial Stress (xx) vs. Bone Volume Fraction', fontsize=16)
    plt.xlabel('Bone Volume Fraction (BV/TV)', fontsize=12)
    plt.ylabel('Stress xx (Pa)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('Plots/straight_beam_uniaxial/initial_stress_xx_vs_bvtv.pdf', dpi=500)


def plot_rve_time_series(cross_section, results):
    """
    Plots the time series of strain_effect_on_OBp, stress_xx, and SED_bm
    for 6 selected RVE elements based on their spatial location in the cross-section.

    :param cross_section: DataFrame containing spatial coordinates.
    :param results: Dictionary containing the solution for each interval and RVE.
    """
    print("--- Generating RVE Time Series Plots ---")

    # 1. Select 6 RVE indices based on distinct spatial locations
    cross_section_df = cross_section.copy()
    # Calculate radius from origin (assuming cross-section is centered at 0,0)
    cross_section_df['r'] = np.sqrt(cross_section_df['y']**2 + cross_section_df['z']**2)

    # Indices based on physical location extremes
    try:
        # 1. Inner Ring (Smallest Radius)
        idx_inner = cross_section_df['r'].idxmin()
        # 2. Outer Ring (Largest Radius)
        idx_outer = cross_section_df['r'].idxmax()
        # 3. Maximum Z (Top side of bone)
        idx_max_z = cross_section_df['z'].idxmax()
        # 4. Minimum Z (Bottom side of bone)
        idx_min_z = cross_section_df['z'].idxmin()
        # 5. Maximum Y (Right side of bone)
        idx_max_y = cross_section_df['y'].idxmax()
        # 6. Minimum Y (Left side of bone)
        idx_min_y = cross_section_df['y'].idxmin()

        # Collect unique indices (in case max/min/inner/outer overlap)
        rve_indices_to_track = list(set([idx_inner, idx_outer, idx_max_z, idx_min_z, idx_max_y, idx_min_y]))
        rve_indices_to_track.sort() # Sort for deterministic plot order

    except Exception as e:
        print(f"Error selecting RVE indices based on spatial location: {e}")
        print("Falling back to sequential indices.")
        rve_indices_to_track = [10, 50, 90, 130, 170, 210]

    # Ensure selected indices are valid
    max_index = len(cross_section) - 1
    valid_indices = [i for i in rve_indices_to_track if i <= max_index]

    if len(valid_indices) < 3:
        print("Error: Too few valid RVE indices selected. Skipping time series plot.")
        return

    # Extract time points (assuming interval keys represent years)
    time_points = sorted(results.keys())
    time_years = np.array(time_points)

    # Define the parameters to plot
    parameters = [
        ('strain_effect_on_OBp', 'Strain Effect on $OB_p$', 'Activity'),
        ('stress_xx', 'Stress $z z$', 'Stress (Pa)'),
        ('SED_bm', 'Strain Energy Density (SED)', 'SED (Pa or J/m$^3$)')
    ]

    # Iterate through parameters and create plots
    for param_key, title_label, y_label in parameters:
        plt.figure(figsize=(10, 6))
        plt.title(f'{title_label} Over Time for Selected RVEs', fontsize=16)
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)

        for RVE_index in valid_indices:
            # Extract the trajectory for this RVE
            rve_time_series = []
            for t_key in time_points:
                try:
                    rve_time_series.append(results[t_key][RVE_index][param_key])
                except KeyError:
                    # Use NaN to handle missing data point if a key isn't present
                    rve_time_series.append(np.nan)

            # Get RVE coordinates for the legend
            y_mm = cross_section.loc[RVE_index, 'y'] * 1e3
            z_mm = cross_section.loc[RVE_index, 'z'] * 1e3
            label = f'RVE {RVE_index} (y={y_mm:.1f}mm, z={z_mm:.1f}mm)'

            # Plot
            plt.plot(time_years, rve_time_series, label=label, marker='o', markersize=4, linestyle='-')

        plt.legend(loc='best', fontsize=9)
        plt.tight_layout()
        plt.savefig(f'Plots/straight_beam_uniaxial/rve_time_series_{param_key}.pdf', dpi=500)


load_case = Lerebours_Load_Case_Spaceflight()
model = Lerebours_Model(load_case, duration_of_simulation=4)

cross_section, results = model.solve_spatial_model()
# plot_rve_time_series(cross_section, results)
# plot_initial_strain_effect_distribution(cross_section, results)
# plot_initial_stress_distribution(cross_section, results)
# plot_initial_sed_distribution(cross_section, results)
# plot_initial_stress_sed_vs_bvtv(results)
plot_bvtv_for_intervals(cross_section, results)
# plot_profiles(cross_section, results, model)
