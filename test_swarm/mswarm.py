import numpy as np
from pyswarm import pso
from model_v2 import *


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

    error_q = get_error('qm3')
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
