# -*- coding: utf-8 -*-


mac = True
if mac:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns

sns.set_context(rc={'lines.markeredgewidth': 0.1})
sns.set(style="whitegrid")

# Get Observed data
PC = False
if PC:
    path = "Data/"
else:
    # path = "C:/Users/DayTimeChunks/Documents/PhD/HydrologicalMonitoring"
    path = "Data/"

path += "qmBlk_R.csv"
obs = pd.read_csv(path, sep=",")


def add_margin(ax, x=0.05, y=0.05):
    # This will, by default, add 5% to the x and y margins. You
    # can customise this using the x and y arguments when you call it.

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xmargin = (xlim[1]-xlim[0])*x
    ymargin = (ylim[1]-ylim[0])*y

    ax.set_xlim(xlim[0]-xmargin, xlim[1]+xmargin)
    ax.set_ylim(ylim[0]-ymargin, ylim[1]+ymargin)


def basic_plot(obs):
    fig, ax1 = plt.subplots()

    obs['DayMoYr'] = pd.to_datetime(obs['DayMoYr'])
    time = obs['DayMoYr']
    obs_vol = obs['VolTot.L']
    sim_vol = obs_vol*.10


    legend = ["m3 observed", "m3 simulated"]
    color_sequence = ['#d62728',
                      'darkviolet'
                      '#2ca02c',
                      '#1f77b4'
                      ]

    ax1.plot(time, obs_vol) #, linestyle='dashdot', label=legend[0])
    ax1.plot(time, sim_vol) #, linestyle='dashed', label=legend[1])

    # add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Volume (m3)')
    plt.show()

def plot_hydro(data, title=None, save=False):
    fig, ax1 = plt.subplots()
    data['DayMoYr'] = pd.to_datetime(data['DayMoYr'])
    time = data['DayMoYr']
    obs_vol = data['VolTot.L']

    legend = ["m3 - observed", "m3 - simulated"]
    color_sequence = ['#d62728',
                      'darkviolet'
                      '#2ca02c',
                      '#1f77b4'
                      ]

    ax1.plot(time, obs_vol, color_sequence[0], linestyle='dashdot', label=legend[0])
    ax1.plot(time, obs_vol*.10, color_sequence[1], linestyle='dashed', label=legend[1])

    # add_margin(ax1, x=0.01, y=0.01)

    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Volume (m3)')

    if title:
        plt.title(title)

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.3,  0.9), fancybox=True, framealpha=0.7, ncol=2)
    if save:
        plt.savefig('../Figures/' + str(title) + '.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

    plt.show()


basic_plot(obs)
