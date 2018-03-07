import numpy as np
import pandas as pd
mac = True
if mac:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

import seaborn as sns
sns.set_context(rc={'lines.markeredgewidth': 0.1})
sns.set(style="whitegrid")


def add_margin(ax, x=0.05, y=0.05):
    # This will, by default, add 5% to the x and y margins. You
    # can customise this using the x and y arguments when you call it.

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xmargin = (xlim[1]-xlim[0])*x
    ymargin = (ylim[1]-ylim[0])*y

    ax.set_xlim(xlim[0]-xmargin, xlim[1]+xmargin)
    ax.set_ylim(ylim[0]-ymargin, ylim[1]+ymargin)


def hydro_plot(data):
    fig, ax1 = plt.subplots()
    time = data['Date']
    obs_vol = data['VolTot.L']
    sim_DPLF = data['LF_sim.L'] #+ data['Perc_sim.L']
    sim_vol = data['Q_sim.L'] #+ data['Perc_sim.L']

    legend = ["Q-obs (L/d)", "Q (L/d)", "DPLF (L/d)"]
    color_sequence = ['#d62728',
                      'darkviolet'
                      '#2ca02c',
                      '#1f77b4'
                      ]

    ax1.plot(time, obs_vol, label=legend[0])
    ax1.plot(time, sim_vol, linestyle='dashdot', label=legend[1])
    ax1.plot(time, sim_DPLF, linestyle='dashed', label=legend[2])
    # ax1.plot(time, obs_vol, color_sequence[0], linestyle='dashdot', label=legend[0])
    title = None
    if title:
        plt.title(title)
    # handles, labels = ax1.get_legend_handles_labels()
    # lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.3, 0.9), fancybox=True, framealpha=0.7,
    #                  ncol=2)
    # if save:
    #     plt.savefig('../Figures/' + str(title) + '.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

    # add_margin(ax1, x=0.01, y=0.01)
    plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Flux (L/day)')
    plt.show()

