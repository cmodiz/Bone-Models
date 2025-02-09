import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# plt.style.use('custom_whitegrid.mplstyle')


# def plot_bone_cell_concentrations(time, OBp, OBa, OCa, title):
#     plt.figure()
#     plt.plot(time, OBp, label=r'$OB_p$', color='#2066a8', linestyle='dotted', linewidth=2)
#     plt.plot(time, OBa, label=r'$OB_a$', color='#2066a8', linewidth=2)
#     plt.plot(time, OCa, label=r'$OC_a$', color='#ae282c', linestyle='dashed', linewidth=2)
#     plt.ylabel('Concentration [pM]')
#     plt.xlabel('Time [days]')
#     plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
#     plt.title(title)
#     plt.show()


def plot_bone_volume_fraction(time, bone_volume_fractions, title):
    plt.figure()
    # time = np.insert(time, 0, 0)
    plt.plot(time, bone_volume_fractions, label='Bone Volume Fraction', color='#2066a8', linewidth=2)
    plt.ylabel('Bone Volume Fraction')
    plt.xlabel('Time [days]')
    plt.title(title)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.show()


def plot_bone_cell_concentrations(solutions, Disease_Case, Reference_Case, calibration_type):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12), sharey=True)
    model_types = ['cellular_responsiveness', 'integrated_activity']
    calibration_types = ['calibration_type_only_for_healthy_state', 'calibration_type_all']
    titles = [[r'$\tilde{k}_R \cdot \alpha_R$', r'$k_R \cdot \alpha_R$'],
              [r'$\tilde{k}_T \cdot \alpha_T$', r'$k_T \cdot \alpha_T$']]

    for i, model_type in enumerate(model_types):
        for j, calibration_type in enumerate(calibration_types):
            time = solutions[Disease_Case][model_type][calibration_type]['t']
            OBp = solutions[Disease_Case][model_type][calibration_type]['y'][0]
            OBa = solutions[Disease_Case][model_type][calibration_type]['y'][1]
            OCa = solutions[Disease_Case][model_type][calibration_type]['y'][2]
            axs[i][j].plot(time, OBp, label=r'$OB_p$', linestyle='dotted', linewidth=2, color='#2066a8')
            axs[i][j].plot(time, OBa, label=r'$OB_a$', linewidth=2, color='#2066a8')
            axs[i][j].plot(time, OCa, label=r'$OC_a$', linestyle='dashed', linewidth=2, color='#ae282c')
            time = solutions[Reference_Case]['old_activation']['calibration_type_all']['t']
            OBp = solutions[Reference_Case]['old_activation']['calibration_type_all']['y'][0]
            OBa = solutions[Reference_Case]['old_activation']['calibration_type_all']['y'][1]
            OCa = solutions[Reference_Case]['old_activation']['calibration_type_all']['y'][2]
            axs[i][j].plot(time, OBp, label=r'$OB_p^{ref}$', linestyle='dotted', linewidth=2, color='#90B8DD')
            axs[i][j].plot(time, OBa, label=r'$OB_a^{ref}$', linewidth=2, color='#90B8DD')
            axs[i][j].plot(time, OCa, label=r'$OC_a^{ref}$', linestyle='dashed', linewidth=2, color='#E39A9C')
            axs[i][j].set_title(titles[i][j])
            axs[i][j].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            axs[i][j].grid(True)
            if j == 0:
                axs[i][j].set_ylabel('Concentration [pM]')
            if i == 1:
                axs[i][j].set_xlabel('Time [days]')

    # Create a shared legend
    handles, labels = axs[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=6, fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def plot_bone_volume_fractions(solutions):
    mpl.rcParams['font.size'] = 18
    fig, axs = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
    model_types = ['cellular_responsiveness', 'integrated_activity', 'old_activation']
    titles = [r'$k_R \cdot \alpha_R$', r'$k_T \cdot \alpha_T$', r'$\pi_{PTH}$']
    colors = ['#59a89c', '#0b81a2', '#7E4794', '#e25759', 'grey', '#9d2c00']
    labels = ['HPT', 'OP', 'PMO', 'HyperC', 'HypoC', 'GIO']
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted', (0, (3, 1, 1, 1)), (0, (5, 10))]

    for i, model_type in enumerate(model_types):
        for j, (disease_name, data) in enumerate(solutions.items()):
            if model_type in data:
                if model_type == 'old_activation':
                    time = data[model_type]['calibration_type_all']['t']
                    bone_volume_fraction = data[model_type]['calibration_type_all']['bone_volume_fraction']
                    axs[i].plot(time, bone_volume_fraction, label=labels[j % len(labels)], linewidth=3, color=colors[j % len(colors)], linestyle=linestyles[j % len(linestyles)])
                else:
                    time = data[model_type]['calibration_type_all']['t']
                    bone_volume_fraction = data[model_type]['calibration_type_all']['bone_volume_fraction']
                    axs[i].plot(time, bone_volume_fraction, label=labels[j % len(labels)], linewidth=3, color=colors[j % len(colors)], linestyle=linestyles[j % len(linestyles)])
        axs[i].set_title(titles[i])
        axs[i].set_xlabel('Time [days]')
        axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axs[i].grid(True)
        if i == 0:
            axs[i].set_ylabel('Bone Volume Fraction [-]')

    # Create a shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=6, fontsize=18)
    plt.tight_layout()
    fig.savefig(r'C:\Users\n12050199\OneDrive - Queensland University of Technology\Bone Model Code\bone_volume_fractions_calibration_all.pdf', bbox_inches='tight')
    plt.show()

def plot_all_model_options(solutions, Disease_Case, Reference_Case):
    mpl.rcParams['font.size'] = 18
    fig, axs = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
    model_types = ['cellular_responsiveness', 'integrated_activity']
    calibration_types = ['calibration_type_only_for_healthy_state', 'calibration_type_all']
    titles = [r'$OB_p$', r'$OB_a$', r'$OC_a$']
    colors = ['#59a89c', '#0b81a2', '#7E4794', '#e25759']
    linestyles = ['solid', 'dashed', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]
    labels = [r'$\tilde{k}_R \cdot \alpha_R$', r'$k_R \cdot \alpha_R$', r'$\tilde{k}_T \cdot \alpha_T$', r'$k_T \cdot \alpha_T$']
    for i, title in enumerate(titles):
        for model_idx, model_type in enumerate(model_types):
            for calib_idx, calibration_type in enumerate(calibration_types):
                time = solutions[Disease_Case][model_type][calibration_type]['t']
                data = solutions[Disease_Case][model_type][calibration_type]['y'][i]
                color = colors[model_idx * len(calibration_types) + calib_idx]
                linestyle = linestyles[model_idx * len(calibration_types) + calib_idx]
                label = labels[model_idx * len(calibration_types) + calib_idx]
                axs[i].plot(time, data, label=label, linewidth=3, color=color, linestyle=linestyle)
        # Add reference data
        time = solutions[Reference_Case]['old_activation']['calibration_type_all']['t']
        data = solutions[Reference_Case]['old_activation']['calibration_type_all']['y'][i]
        axs[i].plot(time, data, label=r'$\pi_{PTH}$', linestyle='dotted', linewidth=3, color='#c8c8c8')
        axs[i].set_title(title)
        axs[i].set_xlabel('Time [days]')
        axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axs[i].grid(True)
        if i == 0:
            axs[i].set_ylabel('Concentration [pM]')
    # Create a shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.085), ncol=5, fontsize=18)
    plt.tight_layout()
    fig.savefig(r'C:\Users\n12050199\OneDrive - Queensland University of Technology\Bone Model Code\bone_cells_plot_HPT.pdf', bbox_inches='tight')
    plt.show()


def plot_PTH_activation_for_all_disease_states(cellular_responsiveness_calibration_all, integrated_activity_calibration_all, cellular_responsiveness_calibration_healthy, integrated_activity_calibration_healthy, old_activation):
    PTH_activations = {}
    PTH_activations[r'$\tilde{k}_R \cdot \alpha_R$'] = cellular_responsiveness_calibration_healthy
    PTH_activations[r'$k_R \cdot \alpha_R$'] = cellular_responsiveness_calibration_all
    PTH_activations[r'$\tilde{k}_T \cdot \alpha_T$'] = integrated_activity_calibration_healthy
    PTH_activations[r'$k_T \cdot \alpha_T$'] = integrated_activity_calibration_all
    PTH_activations[r'$\pi_{PTH}$'] = old_activation

    mpl.rcParams['font.size'] = 18
    colors = ['#59a89c', '#0b81a2', '#7E4794', '#e25759', '#c8c8c8']
    patterns = [ "/" , "x", "-" , "+" , ""]
    disease_states = ['Healthy', 'HPT', 'OP', 'PMO', 'HyperC', 'HypoC', 'GIO']
    x = np.arange(len(disease_states))  # the label locations
    width = 0.15  # the width of the bars
    spacing = 0.05  # space between groups
    multiplier = 0

    fig, ax = plt.subplots(figsize=(18, 7))
    for attribute, measurement in PTH_activations.items():
        offset = width * multiplier
        ax.bar(x + offset, measurement, width, label=attribute, color=colors[multiplier], hatch=patterns[multiplier])
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Values [-]')
    ax.set_xticks(x + (width + spacing) * (len(PTH_activations) - 1) / 2)
    ax.set_xticklabels(disease_states)
    ax.legend(loc='upper left', ncols=3)
    fig.savefig(r'C:\Users\n12050199\OneDrive - Queensland University of Technology\Bone Model Code\calibration_results.pdf', bbox_inches='tight')
    plt.tight_layout()
    plt.show()