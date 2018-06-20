# -*- coding: utf-8 -*-

import pandas as pd


def getTime(sep=","):
    path = "Data/Time.csv"
    obs = pd.read_csv(path, sep=sep)
    return obs


def getTimeStamp(timestep, sep=','):
    path = "Data/Time.csv"
    obs = pd.read_csv(path, sep=sep)
    obs = obs[['Jdays', 'Date']]
    obs_dict = obs.to_dict(orient='split')
    return str(obs_dict['data'][timestep-1][1])

day = getTimeStamp(177, ";")
print(day.split("/")[2])


def getWaterData(tss=True):

    path = "Data/"
    if tss:
        path += "q_obs_m3day.tss"
        obs = pd.read_table(path, header=None)
        obs = obs.rename(index=str, columns={0: "Jdays", 1: "Qm3"})
    else:
        path += "qmBlk_R.csv"
        obs = pd.read_csv(path, sep=",")

    return obs


def getOutletData(name, var=None, tss=True):
    path = "Data/BEACH_R/"
    if var is None:
        var = name
    if tss:
        path += name + ".tss"
        obs = pd.read_table(path, header=None)
        obs = obs.rename(index=str, columns={0: "Jdays", 1: var})
    else:
        path += name + ".csv"
        obs = pd.read_csv(path, sep=",")

    return obs


def getETP():
    path = "Data/"
    path += "ET0.tss"
    obs = pd.read_table(path, header=None, skiprows=4)
    obs = obs.rename(index=str, columns={0: "Jdays", 1: "ETPmm"})
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


