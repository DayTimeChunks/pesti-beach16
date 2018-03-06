# -*- coding: utf-8 -*-


mac = True
if mac:
    import matplotlib
    matplotlib.use('TkAgg')
    plt = matplotlib.pyplot
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
    path = "../Data/"
else:
    # path = "C:/Users/DayTimeChunks/Documents/PhD/HydrologicalMonitoring"
    path = "../Data/"

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
