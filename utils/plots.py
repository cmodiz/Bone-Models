import matplotlib.pyplot as plt
import numpy as np

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
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    for i, model_type in enumerate(['cellular_responsiveness', 'integrated_activity']):
        if model_type in solutions[Disease_Case]:
            time = solutions[Disease_Case][model_type][calibration_type]['t']
            OBp = solutions[Disease_Case][model_type][calibration_type]['y'][0]
            OBa = solutions[Disease_Case][model_type][calibration_type]['y'][1]
            OCa = solutions[Disease_Case][model_type][calibration_type]['y'][2]
            axs[i].plot(time, OBp, label=r'$OB_p$', linestyle='dotted', linewidth=2, color='#2066a8')
            axs[i].plot(time, OBa, label=r'$OB_a$', linewidth=2, color='#2066a8')
            axs[i].plot(time, OCa, label=r'$OC_a$', linestyle='dashed', linewidth=2, color='#ae282c')
        axs[i].set_title(model_type)
        axs[i].set_xlabel('Time [days]')
        axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axs[i].grid(True)
        if i == 0:
            axs[i].set_ylabel('Concentration [pM]')
        axs[i].legend()

    time = solutions[Reference_Case]['old_activation']['t']
    OBp = solutions[Reference_Case]['old_activation']['y'][0]
    OBa = solutions[Reference_Case]['old_activation']['y'][1]
    OCa = solutions[Reference_Case]['old_activation']['y'][2]
    axs[2].plot(time, OBp, label=r'$OB_p$', linestyle='dotted', linewidth=2, color='#2066a8')
    axs[2].plot(time, OBa, label=r'$OB_a$', linewidth=2, color='#2066a8')
    axs[2].plot(time, OCa, label=r'$OC_a$', linestyle='dashed', linewidth=2, color='#ae282c')
    axs[2].set_title('old_activation')
    axs[2].set_xlabel('Time [days]')
    axs[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.show()


def plot_bone_volume_fractions(solutions, Disease_Cases_List, model_type='cellular_responsiveness', calibration_type='all'):
    plt.figure()
    for disease in Disease_Cases_List:
        disease_name = disease.__name__
        time = solutions[disease_name][model_type][calibration_type]['t']
        bone_volume_fractions = solutions[disease_name][model_type][calibration_type]['bone_volume_fraction']
        plt.plot(time, bone_volume_fractions, label=f'{disease_name}', linewidth=2)
    plt.ylabel('Bone Volume Fraction')
    plt.xlabel('Time [days]')
    plt.title(model_type)
    plt.legend()
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.show()
