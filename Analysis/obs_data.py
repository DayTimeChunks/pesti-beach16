# -*- coding: utf-8 -*-

import pandas as pd

# Get Observed data
PC = False
if PC:
    path = "Data/"
else:
    # path = "C:/Users/DayTimeChunks/Documents/PhD/HydrologicalMonitoring"
    path = "Data/"

use_tss = False
if use_tss:
    path += "q_obs_m3day.tss"
    obs = pd.read_table(path)
else:
    path += "qmBlk_R.csv"
    obs = pd.read_csv(path, sep=",")

test = 'test'


