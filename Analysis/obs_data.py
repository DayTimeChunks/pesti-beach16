# -*- coding: utf-8 -*-

import pandas as pd

# Get Observed data
PC = False
if PC:
    path = "Data/"
else:
    # path = "C:/Users/DayTimeChunks/Documents/PhD/HydrologicalMonitoring"
    path = "Data/"

path += "qmBlk_R.csv"
obs = pd.read_csv(path, sep=",")

test = 'test'


