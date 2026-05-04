from bone_models.bone_mineralisation_models.models.modiz_model import Modiz_Model
from bone_models.bone_mineralisation_models.load_cases.modiz_load_cases import Modiz_Load_Case
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import minimize
from matplotlib.colors import LinearSegmentedColormap

# --------------- Set up the plot style --------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,   # embed TrueType fonts
    "ps.fonttype": 42
})
figsize_single = (90/25.4, 90/25.4 * 0.75)  # width, height
figsize_double = (190/25.4, 190/25.4 * 0.55)

plt.rcParams.update({
    "font.size": 7,          # base font
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 7
})

plt.rcParams.update({
    "lines.linewidth": 1.0,
    "lines.markersize": 4,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3
})

# -------------- Calculate the fitting parameter for specific surface to turnover (Figure 3) --------------
calculate_fitting_parameter_specific_surface_to_turnover = False
if calculate_fitting_parameter_specific_surface_to_turnover:
    load_case = Modiz_Load_Case()
    model = Modiz_Model(load_case, 0)
    fitting_parameter = model.fit_specific_surface_to_turnover_data(show_plot=True)
    print(f"Fitting parameter: {fitting_parameter}")

# -------------- Plot the specific surface curve to include dispersion --------------
plot_specific_surface_with_dispersion = False
if plot_specific_surface_with_dispersion:
    load_case = Modiz_Load_Case()
    model = Modiz_Model(load_case, 0)
    bone_volume_fractions = np.linspace(0, 1, 200)
    specific_surface_values = model.lerebours_model.specific_surface(bone_volume_fractions)

    plt.figure(figsize=figsize_single)
    plt.plot(bone_volume_fractions, 1.4826 * specific_surface_values, color='#4a4a4aff',
             label='1.4826 x Specific Surface', linestyle='dashed')
    plt.plot(bone_volume_fractions, specific_surface_values, color='#4a4a4aff', label='Specific Surface')
    plt.plot(bone_volume_fractions, 0.723217 * specific_surface_values, color='#4a4a4aff',
             label='0.723217 x Specific Surface', linestyle='dashed')
    plt.xlabel('Bone Volume Fraction [-]')
    plt.ylabel(r'Specific Surface $\left[\frac{1}{mm}\right]$')
    plt.title('Specific Surface vs Bone Volume Fraction')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('specific_surface_plot_with_dispersion.pdf')
    plt.show()

# -------------- Plotting the experimental data (Figure 1) --------------
plot_experimental_data_color_coded = False
plot_experimental_data_annotated = False
base_dir = os.path.dirname(__file__)
relative_path = '../utils/data/Zioupos_Experimental_Data.xlsx'
file_path = os.path.join(base_dir, relative_path)
experimental_data = pd.read_excel(file_path, sheet_name=1)
porosity = experimental_data.iloc[:, 3]
bone_volume_fraction = 1 - porosity / 100
if plot_experimental_data_annotated:
    plt.figure(figsize=(7, 4.5))
    plt.scatter(experimental_data.iloc[:, 2], experimental_data.iloc[:, 1], label='Experimental Data',
                color='#3d92ceff')
    for i, bvf in enumerate(bone_volume_fraction):
        if i % 4 == 0:
            plt.annotate(f'{bvf:.2f}', (experimental_data.iloc[i, 2], experimental_data.iloc[i, 1]))
    plt.ylabel(r'Apparent Density g/cm$^3$')
    plt.xlabel(r'Material Density g/cm$^3$')
    plt.grid(True)
    plt.savefig('material_vs_apparent_density_experimental_data.pdf')
    plt.show()
if plot_experimental_data_color_coded:
    # Define a colormap
    blues_truncated = plt.cm.Blues(np.linspace(0.3, 1.0, 256))
    colormap = LinearSegmentedColormap.from_list('blues_dark', blues_truncated)
    plt.figure(figsize=figsize_single)
    scatter = plt.scatter(experimental_data.iloc[:, 2], experimental_data.iloc[:, 1],
                           c=bone_volume_fraction, cmap=colormap, label='Experimental Data')
    # Define target values to annotate
    targets = [0.94, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    for target in targets:
        closest_idx = np.argmin(np.abs(bone_volume_fraction - target))
        closest_value = bone_volume_fraction[closest_idx]
        plt.annotate(f'{closest_value:.2f}',
            (experimental_data.iloc[closest_idx, 2],
             experimental_data.iloc[closest_idx, 1]),
                     fontsize=6,
            xytext=(5, 5),  # Offset the text
            textcoords='offset points',
            bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9, edgecolor='black',
                               linewidth=0.4),
            arrowprops=dict(
                         arrowstyle='-',
                         linewidth=0.4,
                         connectionstyle='arc3,rad=0'
                     ))
    max_idx = np.argmax(bone_volume_fraction)
    max_value = bone_volume_fraction[max_idx]
    plt.annotate(f'{max_value:.2f}',
                (experimental_data.iloc[max_idx, 2],
                 experimental_data.iloc[max_idx, 1]), xytext=(5, 5),fontsize=6,  # Offset the text
                textcoords='offset points',
                bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9, edgecolor='black',linewidth=0.4),
                 arrowprops=dict(
                     arrowstyle='-',
                     linewidth=0.4,
                     connectionstyle='arc3,rad=0'
                 ))
    cbar = plt.colorbar(scatter)
    cbar.set_label('Bone Volume Fraction [-]')
    plt.ylabel(r'Apparent Density [g/cm$^3$]')
    plt.xlabel(r'Material Density [g/cm$^3$]')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.savefig('material_vs_apparent_density_experimental_data_colorcoded.pdf', bbox_inches='tight')
    plt.show()

# -------------- Plotting the experimental data from Martin Specific Surface (Figure 6) --------------
plot_specific_surface_dispersion = False
if plot_specific_surface_dispersion:
    base_dir = os.path.dirname(__file__)
    relative_path = '../utils/data/Martin_Data.xlsx'
    file_path = os.path.join(base_dir, relative_path)
    Martin_data = pd.ExcelFile(file_path)
    sheet_names = Martin_data.sheet_names
    marker_styles = {
        'Full Circles': 'o',  # Filled circle
        'Empty Circles': 'o',  # Empty circle
        'Triangles': '^',  # Triangle
        'Full Square': 's',  # Square
        'Empty Square': 's',  # Empty square
        'Diamond': 'D'  # Diamond
    }
    fill_styles = {
        'Full Circles': 'full',
        'Empty Circles': 'none',
        'Triangles': 'full',
        'Full Square': 'full',
        'Empty Square': 'none',
        'Diamond': 'full'
    }
    plt.figure(figsize=figsize_single)
    for sheet_name in sheet_names:
        df = pd.read_excel(Martin_data, sheet_name=sheet_name)
        x_values = df.iloc[:, 0].values
        y_values = df.iloc[:, 1].values
        marker = marker_styles.get(sheet_name, 'o')  # Default to circle if not specified
        fillstyle = fill_styles.get(sheet_name, 'full')  # Default to filled if not specified
        if sheet_name == 'Sv':
            continue
        if fillstyle == 'none':
            plt.scatter(1-x_values, y_values,
                        marker=marker,
                        facecolors='none',
                        edgecolors='black',
                        label=sheet_name,
                        linewidths=0.6,
                        s=8,
                        alpha=0.7)  # Transparency
        else:
            plt.scatter(1-x_values, y_values,
                        marker=marker,
                        color='black',
                        label=sheet_name,
                        s=8,
                        alpha=0.7)  # Transparency

    bone_volume_fraction = np.linspace(0, 1, 200)
    load_case = Modiz_Load_Case()
    model_1 = Modiz_Model(load_case, 0, specific_surface_multiplier=1.459167495)
    model_2 = Modiz_Model(load_case, 0, specific_surface_multiplier=1)
    model_3 = Modiz_Model(load_case, 0, specific_surface_multiplier=0.734321304)
    specific_surface_1 = model_1.lerebours_model.specific_surface(bone_volume_fraction)
    specific_surface_2 = model_2.lerebours_model.specific_surface(bone_volume_fraction)
    specific_surface_3 = model_3.lerebours_model.specific_surface(bone_volume_fraction)

    plt.plot(bone_volume_fraction, specific_surface_1, label='Multiplier = 1.459167495', linestyle='dashed', color='#4a4a4aff')
    plt.plot(bone_volume_fraction, specific_surface_2, label='Multiplier = 1', linestyle='solid', color='#4a4a4aff')
    plt.plot(bone_volume_fraction, specific_surface_3, label='Multiplier = 0.734321304', linestyle='dashed', color='#4a4a4aff')

    plt.ylabel('Specific Surface [mm$^{-1}$]')
    plt.xlabel('Bone Volume Fraction [-]')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig('Sv_curves_with_dispersion.pdf', bbox_inches='tight')
    plt.show()

# -------------- Fitting the mineral apposition rate --------------
def objective_function(apposition_rate, model, data_points):
    total_error = 0
    for bone_volume_fraction, target_material_density, target_apparent_density in data_points:
        target_material_density_range = data_points[:, 1].max() - data_points[:, 1].min()
        target_apparent_density_range = data_points[:, 2].max() - data_points[:, 2].min()
        # Set the apposition rate in the model
        model.parameters.mineralisation.reference_apposition_rate = apposition_rate
        # Calculate the model's material and apparent densities
        predicted_material_density, predicted_apparent_density = model.calculate_mineral_densities(bone_volume_fraction)
        # Calculate the error for this data point
        error = ((predicted_material_density - target_material_density)/target_material_density_range) ** 2 + \
                ((predicted_apparent_density - target_apparent_density)/target_apparent_density_range) ** 2
        total_error += error
    return total_error


fit_mineral_apposition_rate = False
if fit_mineral_apposition_rate:
    data_points = experimental_data.iloc[:, [3, 2, 1]].values  # Porosity, Material Density, Apparent Density
    data_points = pd.DataFrame(data_points).dropna().values
    data_points[:, 0] = 1 - data_points[:, 0] / 100  # Convert porosity to bone volume fraction
    load_case = Modiz_Load_Case()
    model = Modiz_Model(load_case, 0)  # Bone volume fraction will be set dynamically
    initial_apposition_rate = 0.000363855075413035
    result = minimize(
        objective_function,
        x0=initial_apposition_rate,
        args=(model, data_points),
        bounds=[(1e-5, 1e-3)],
        options={'disp': True}
    )
    # Optimized apposition rate
    optimized_apposition_rate = result.x[0]
    print(f"Optimized reference apposition rate: {optimized_apposition_rate}")

# -------------- Initialise model and save results in lists --------------
load_case = Modiz_Load_Case()
material_density_list = []
apparent_density_list = []
ageing_queues = []
mineral_queues = []
OCu_list, OBu_list, OCa_list, OBa_list = [], [], [], []
all_ageing_queues, all_mineralisation_values, all_mineral_queues = {}, {}, {}
bone_volume_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98]
# bone_volume_fractions = [0.3, 0.8]
fig, axs = plt.subplots(2, 2, figsize=(7, 5.5))
line_styles = {0.3: 'dashed', 0.8: 'solid'}
colors = {0.3: '#b4196fff', 0.8: '#3d92ceff'}

time_vector = np.arange(0, 40000)
for bone_volume_fraction in bone_volume_fractions:
    model = Modiz_Model(load_case, bone_volume_fraction)
    if bone_volume_fraction >= 0.8:
        model.parameters.mineralisation.length_of_queue = 45000
    material_density, apparent_density = model.calculate_mineral_densities(bone_volume_fraction)
    material_density_list.append(material_density)
    apparent_density_list.append(apparent_density)
    OCu_list.append(model.lerebours_model.steady_state.OCu)
    OBu_list.append(model.lerebours_model.steady_state.OBu)
    OCa_list.append(model.lerebours_model.steady_state.OCa)
    OBa_list.append(model.lerebours_model.steady_state.OBa)
    mineral_queue = model.calculate_mineral_queue(model.ageing_queue, bone_volume_fraction)
    mineralisation_values = [model.calculate_mineralisation_law(t, bone_volume_fraction) for t in time_vector]

    non_zero_indices = [i for i, value in enumerate(model.ageing_queue) if value != 0]
    non_zero_ageing_queue = [model.ageing_queue[i] for i in non_zero_indices]
    non_zero_mineral_queue = [mineral_queue[i] for i in non_zero_indices]

    all_ageing_queues[bone_volume_fraction] = (non_zero_ageing_queue, non_zero_mineral_queue)
    all_mineralisation_values[bone_volume_fraction] = mineralisation_values
    all_mineral_queues[bone_volume_fraction] = (non_zero_mineral_queue, non_zero_ageing_queue)

# ----------- Plot Queuing Details -----------
plot_queuing_details = False
# --------- Plot Ageing Queue -----------
if plot_queuing_details:
    plt.figure(figsize=figsize_single)
    for bone_volume_fraction in reversed(bone_volume_fractions):
        non_zero_ageing_queue, _ = all_ageing_queues[bone_volume_fraction]
        plt.bar(range(len(non_zero_ageing_queue)), non_zero_ageing_queue, width=0.6, alpha=0.5,
                color=colors[bone_volume_fraction], label=f'$f_{{\\mathrm{{Tt.B}}}}={bone_volume_fraction}$')
    plt.xlabel(r'Element Index $j$ [-]')
    plt.ylabel(r'$B_j^{t_R}$ [-]')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.xlim(-10, 30000)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join('ageing_queue.pdf'), bbox_inches='tight')
    plt.close()

    # --------- Plot Mineralisation Law -----------
    plt.figure(figsize=figsize_single)
    for bone_volume_fraction in reversed(bone_volume_fractions):
        mineralisation_values = all_mineralisation_values[bone_volume_fraction]
        plt.plot(time_vector, mineralisation_values, linestyle=line_styles[bone_volume_fraction],
                 color=colors[bone_volume_fraction], label=f'$f_{{\\mathrm{{Tt.B}}}}={bone_volume_fraction}$')
    plt.xlabel(r'Time $t$ [days]')
    plt.ylabel(r'$M(t)$')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.xlim(-10, 30000)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join('mineralisation_law.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()

    # --------- Plot Mineral Queue -----------
    plt.figure(figsize=figsize_single)
    for bone_volume_fraction in reversed(bone_volume_fractions):
        non_zero_mineral_queue, _ = all_mineral_queues[bone_volume_fraction]
        plt.bar(range(len(non_zero_mineral_queue)), non_zero_mineral_queue, width=0.6, alpha=0.5,
                color=colors[bone_volume_fraction], label=f'$f_{{\\mathrm{{Tt.B}}}}={bone_volume_fraction}$')
    plt.xlabel(r'Element Index $j$ [-]')
    plt.ylabel(r'$B_j^{t_R} \cdot M(j)$')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.xlim(-10, 30000)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join('mineral_queue.pdf'), bbox_inches='tight')
    plt.close()

    # --------- Plot Boomerang-shaped plot -----------
    plt.figure(figsize=figsize_single)
    for bone_volume_fraction in reversed(bone_volume_fractions):
        non_zero_mineral_queue, non_zero_ageing_queue = all_mineral_queues[bone_volume_fraction]
        plt.plot(non_zero_mineral_queue, non_zero_ageing_queue, linestyle=line_styles[bone_volume_fraction],
                 color=colors[bone_volume_fraction], label=f'$f_{{\\mathrm{{Tt.B}}}}={bone_volume_fraction}$')
    plt.xlabel(r'$B_j^{t_R} \cdot M(j)$')
    plt.ylabel(r'$B_j^{t_R}$ [-]')
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join('bone_volume_vs_mineral_content.pdf'), bbox_inches='tight')
    plt.close()

# -------------- Plotting the hypothesis 2.1 and 2.2 --------------
plot_results_for_hypothesis_2_1_and_2_2 = False
if plot_results_for_hypothesis_2_1_and_2_2:
    load_case = Modiz_Load_Case()
    material_density_list = []
    apparent_density_list = []
    bone_volume_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98]

    time_vector = np.arange(0, 40000)
    hypotheses = [2.1, 2.2]
    for hypothesis in hypotheses:
        for bone_volume_fraction in bone_volume_fractions:
            model = Modiz_Model(load_case, bone_volume_fraction, hypothesis=hypothesis)
            if bone_volume_fraction >= 0.8:
                model.parameters.mineralisation.length_of_queue = 45000
            material_density, apparent_density = model.calculate_mineral_densities(bone_volume_fraction)
            material_density_list.append(material_density)
            apparent_density_list.append(apparent_density)

    # Plot material and apparent density for both hypotheses
    plt.figure(figsize=figsize_single)
    plt.plot(material_density_list[:len(bone_volume_fractions)],
                apparent_density_list[:len(bone_volume_fractions)],
                label='Age Bias', color='#b4196fff', marker='d', markersize=3, linestyle='--',)
    plt.plot(material_density_list[len(bone_volume_fractions):],
                apparent_density_list[len(bone_volume_fractions):],
                label='Restricted Resorption', color='#620c54ff', marker='*', markersize=3, linestyle='--',)
    plt.scatter(experimental_data.iloc[:, 2], experimental_data.iloc[:, 1],
                label='Experimental Data', alpha=0.6, s=8,
                edgecolors='#3d92ceff', facecolors='#3d92ceff', linewidth=0.6)
    for i, bvf in enumerate(bone_volume_fractions):
        plt.annotate(f'{bvf:.2f}',
                     (material_density_list[i], apparent_density_list[i]),
                     color='black')
    for i, bvf in enumerate(bone_volume_fractions):
        plt.annotate(f'{bvf:.2f}',
                     (material_density_list[len(bone_volume_fractions)+i], apparent_density_list[len(bone_volume_fractions)+i]),
                     color='black')
    plt.xlabel(r'Material Density [g/cm$^3$]')
    plt.ylabel(r'Apparent Density [g/cm$^3$]')
    plt.legend()
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig('material_vs_apparent_density_hypotheses_21_and_22.pdf', bbox_inches='tight')
    plt.show()

#  -------------- Plotting the steady-state cell concentrations (Figure 2) --------------
plot_steady_state_cell_concentrations = False
if plot_steady_state_cell_concentrations:
    fig, ax1 = plt.subplots(figsize=figsize_single)
    line1, = ax1.plot(bone_volume_fractions, OBa_list, color='#3699c7ff', linestyle='solid', label=r'$C_\mathrm{OB_a}$')
    line2, = ax1.plot(bone_volume_fractions, OCa_list, color='#e574c1ff', linestyle='solid', label=r'$C_\mathrm{OC_a}$')
    ax1.set_xlabel('Bone Volume Fraction [-]')
    ax1.set_ylabel('Active Cell Concentrations [pM]', color='black')
    ax1.set_xlim([0.001, 0.999])
    ax1.tick_params(axis='y', labelcolor='black')
    ax2 = ax1.twinx()
    line3, = ax2.plot(bone_volume_fractions, OBu_list, color='#bcd9e5ff', linestyle='dashed',
                      label=r'$C_\mathrm{OB_u}$')
    line4, = ax2.plot(bone_volume_fractions, OCu_list, color='#e6c1c8ff', linestyle='dashed',
                      label=r'$C_\mathrm{OC_u}$')
    ax2.set_ylabel('Uncommitted Cell Concentrations [pM]', color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    lines = [line1, line2, line3, line4]
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels)
    ax1.grid(True, linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.savefig('steady_state_cell_concentrations.pdf', bbox_inches='tight')
    plt.show()

# -------------- Plotting the apparent density vs material density --------------
plot_apparent_density_vs_material_density = False
if plot_apparent_density_vs_material_density:
    plt.figure(figsize=figsize_single)
    plt.plot(material_density_list, apparent_density_list,
             color='#b4196fff', marker='d', markersize=3, linestyle='--', label='Model Results')
    for i, bvf in enumerate(bone_volume_fractions):
        plt.annotate(f'{bvf:.2f}', (material_density_list[i], apparent_density_list[i]), color='black')
    plt.xlabel(r'Material Density [g/cm$^3$]')
    plt.ylabel(r'Apparent Density [g/cm$^3$]')
    plt.scatter(experimental_data.iloc[:, 2], experimental_data.iloc[:, 1], label='Experimental Data', alpha=0.6, s=8,
                edgecolors='#3d92ceff', facecolors='#3d92ceff', linewidth=0.6)
    plt.legend()
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.savefig('material_vs_apparent_density_hypo1.pdf', bbox_inches='tight')
    plt.show()

# -------------- Plotting the apparent density vs material density for specific surface boundaries --------------
plot_apparent_density_vs_material_density_for_boundaries = False
if plot_apparent_density_vs_material_density_for_boundaries:
    bone_volume_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98]
    material_density_high, apparent_density_high = [], []
    material_density_low, apparent_density_low = [], []
    high_multiplier = 1.459167495
    low_multiplier = 0.734321304
    for bone_volume_fraction in bone_volume_fractions:
        model_high = Modiz_Model(load_case, bone_volume_fraction, specific_surface_multiplier=high_multiplier)
        mat_density_high, app_density_high = model_high.calculate_mineral_densities(bone_volume_fraction)
        material_density_high.append(mat_density_high)
        apparent_density_high.append(app_density_high)
        model_low = Modiz_Model(load_case, bone_volume_fraction, specific_surface_multiplier=low_multiplier)
        mat_density_low, app_density_low = model_low.calculate_mineral_densities(bone_volume_fraction)
        material_density_low.append(mat_density_low)
        apparent_density_low.append(app_density_low)

    plt.figure(figsize=figsize_single)
    plt.plot(material_density_high, apparent_density_high, linestyle='dashed', color='#b4196fff', marker='d',
             markersize=3, label='Model Results - Upper Bound')
    for i, bvf in enumerate(bone_volume_fractions):
        plt.annotate(f'{bvf:.2f}', (material_density_high[i], apparent_density_high[i]), color='black')
    plt.plot(material_density_low, apparent_density_low, label='Model Results - Lower Bound', color='#b4196fff', linestyle='dashed', marker='s',
             markersize=3)
    for i, bvf in enumerate(bone_volume_fractions):
        plt.annotate(f'{bvf:.2f}', (material_density_low[i], apparent_density_low[i]), color='black')

    plt.scatter(experimental_data.iloc[:, 2], experimental_data.iloc[:, 1], label='Experimental Data', alpha=0.6, s=8,
                edgecolors='#3d92ceff', facecolors='#3d92ceff', linewidth=0.6)
    plt.xlabel(r'Material Density [g/cm$^3$]')
    plt.ylabel(r'Apparent Density [g/cm$^3$]')
    plt.legend()
    plt.grid(True, linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig('apparent_vs_material_density_boundaries.pdf', bbox_inches='tight')
    plt.show()

# -------------- Plot mineralisation law for exemplary BVF --------------
plot_mineralisation_law = False
if plot_mineralisation_law:
    time_range = np.linspace(0, 5000, 300)
    bone_volume_fraction = 0.2
    mineralisation_values = [model.calculate_mineralisation_law(t, bone_volume_fraction) for t in time_range]
    plt.figure(figsize=(7, 4.5))
    plt.plot(time_range, mineralisation_values, color='#7c4d76ff', label='Mineralisation Law')
    plt.xlabel('Time [days]')
    plt.ylabel('Mineral content')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('mineralisation_law_bvf_0.2.pdf')
    plt.show()


