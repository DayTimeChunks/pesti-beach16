# -*- coding: utf-8 -*-

import pandas as pd


def getWaterData(tss=True):

    path = "Data/"
    if tss:
        path += "q_obs_m3day.tss"
        obs = pd.read_table(path)
    else:
        path += "qmBlk_R.csv"
        obs = pd.read_csv(path, sep=",")

    return obs


def getSoilData(transect, tss=True):
    path = "Data/"

    if tss:
        if transect == "North":
            path += "north.tss"
            obs = pd.read_table(path)
        elif transect == "Valley":
            path += "valley.tss"
            obs = pd.read_table(path)
        else:
            path += "south.tss"
            obs = pd.read_table(path)
    else:
        path += "qmBlk_R.csv"
        obs = pd.read_csv(path, sep=",")

    return obs


def getDetailed(name):
    # Imports only tss files (read_table)
    path = "Data/BEACH_R/"
    path += name
    return pd.read_table(path)

