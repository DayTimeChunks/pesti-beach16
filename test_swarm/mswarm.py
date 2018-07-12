import numpy as np
import pandas as pd
from pyswarm import pso
from model_v2 import *

def getOutletData(name, var=None, tss=True):
    path = '../Analysis/Data/BEACH_R/'
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


def get_outlet_merged():

    # Simulated
    filename = "resW_" + "tot_accVOl_m3" + ".tss"
    out_sim = pd.read_table(filename,
                            skiprows=4, delim_whitespace=True,
                            names=['Jdays', 'SIM'],
                            header=None)

    # Observed
    obs = getOutletData('q_obs_m3day', var='m3d')


    # Merge
    match = pd.merge(obs, out_sim, how='inner', on='Jdays')
    return match


def get_error(var):
    """

    :param var: 'm3d'
    :return:
    """
    # Considers only outlet volume discharged
    # Define observation dataframes to use

    df1 = get_outlet_merged()

    # Nash
    mean = df1[var].mean()
    # Diff sim vs. obs
    dfn = df1.assign(diff_sim=lambda row: (row['SIM'] - row[var]) ** 2)
    err_sim = dfn['diff_sim'].sum()
    # Variance
    dfn = dfn.assign(diff_obs=lambda row: (row[var] - mean) ** 2)
    err_obs = dfn['diff_obs'].sum()
    error = err_sim / err_obs  # nash = 1 - error

    return error


def objective(x):
    """
    1) run the model with x
    2) get an error result, return
    :param x:
    :return:
    """
    nTimeSteps = 200  # 360
    myAlteck16 = BeachModel("clone_nom.map", x)
    dynamicModel = DynamicFramework(myAlteck16,
                                    lastTimeStep=nTimeSteps,
                                    firstTimestep=firstTimeStep)
    dynamicModel.run()

    error_q = get_error('m3d')
    print(error_q)
    return error_q


def clf_constraint(x):
    cZ0Z1 = x[0]
    cZ = x[1]
    return cZ0Z1 - cZ


def gamma_constraint(x):
    gamma01 = x[2]
    gamma23 = x[3]
    return gamma01 - gamma23


def fc_constraints(x):
    fcZ2 = x[4]
    fcZ = x[4]
    return fcZ2 - fcZ


def sat_constraints(x):
    satZ2 = x[4]
    satZ = x[4]
    return satZ2 - satZ


def wp_constraints(x):
    WPz2 = x[1]
    WPz = x[1]
    return WPz2 - WPz


lb = [0.75,
      0.01, 0.01,
      0.0001,
      1100,
      0.001, 0.001,
      0.1,
      0.34, 0.34,
      0.54, 0.47,
      0.13, 0.13
      ]

ub = [0.99,
      1, 1,
      1,
      3650,
      1, 1,
      1,
      0.41, 0.41,
      0.61, 0.53,
      0.19, 0.19
      ]

cons = [clf_constraint, gamma_constraint, fc_constraints, sat_constraints, wp_constraints]
