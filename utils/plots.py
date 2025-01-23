import matplotlib.pyplot as plt
# plt.style.use('custom_whitegrid.mplstyle')


def plot_bone_cell_concentrations(time, OBp, OBa, OCa, title):
    plt.figure()
    plt.plot(time, OBp, label=r'$OB_p$', color='#2066a8', linestyle='dotted', linewidth=2)
    plt.plot(time, OBa, label=r'$OB_a$', color='#2066a8', linewidth=2)
    plt.plot(time, OCa, label=r'$OC_a$', color='#ae282c', linestyle='dashed', linewidth=2)
    plt.ylabel('Concentration [pM]')
    plt.xlabel('Time [days]')
    plt.title(title)
    plt.show()