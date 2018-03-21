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


def hydro_dt(data, m3=False):
    fig, ax1 = plt.subplots()
    # time = data['Date']
    time = data['Jdays']

    color_sequence = ['#d62728',
                      'darkviolet'
                      '#2ca02c',
                      '#1f77b4'
                      ]

    if m3:
        obs_vol = data['Q_obs.m3']
        sim_vol = data['Q_sim.m3']
        legend = ["Q-obs (m3/d)", "Q-sim (m3/d)"]

        ax1.plot(time, obs_vol, label=legend[0])
        ax1.plot(time, sim_vol, linestyle='dashdot', label=legend[1])
        ylabel = 'm3/day'

    else:
        obs_vol = data['VolTot.L']
        sim_DPLF = data['LF_sim.L'] #+ data['Perc_sim.L']
        sim_vol = data['Q_sim.L'] #+ data['Perc_sim.L']

        legend = ["Q-obs (L/d)", "Q (L/d)", "DPLF (L/d)"]

        ax1.plot(time, obs_vol, label=legend[0])
        ax1.plot(time, sim_vol, linestyle='dashdot', label=legend[1])
        ax1.plot(time, sim_DPLF, linestyle='dashed', label=legend[2])
        ylabel = 'L/day'

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
    ax1.set_ylabel(ylabel)
    plt.show()


def nash(data, n_tests):
    fig, ax1 = plt.subplots()
    # time = data.iloc['Jdays'][200:]
    time = data.iloc[:, 0]
    # n1 = data['Nash.c1'][200:]

    legend = []
    for i in range(1, n_tests+1):
        l = "N" + str(i)
        legend.append(l)

    for i in range(1, n_tests + 1):
        ax1.plot(time, data.iloc[:, i], label=legend[i-1])

    ylabel = 'Nash'
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel(ylabel)
    plt.show()
