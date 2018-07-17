from SALib.sample import morris as mos
# from morris_analysis import *
import numpy as np


def saveMorris(Names, SiList, filename="MorrisResults.json"):
    import json
    mListSi = {}
    for i in range(len(Names)):
        var_name = Names[i]
        Si = SiList[i]
        mSi = dict()
        mSi["mu"] = Si["mu"].tolist()
        mSi["mu_star"] = Si["mu_star"].tolist()
        mSi['mu_star_conf'] = Si['mu_star_conf']
        mSi['sigma'] = Si['sigma'].tolist()
        mSi['names'] = Si['names']
        mListSi[var_name] = mSi

    with open(filename, 'w') as f:
        json.dump(mListSi, f)


def saveInputMatrix(param_values):
    np.savetxt("input_vectors.txt", param_values)


def scale_morris(bounds):
    scaled = []
    for e in bounds:
        se = []
        up = max(e[0], e[1])
        se.append(e[0] / up)
        se.append(e[1] / up)
        scaled.append(se)

    return scaled


def get_upper(bounds):
    upper = []
    for e in bounds:
        upper.append(e[1])
    return upper


Mini_TEST = False
if Mini_TEST:
    bounds = [[0.75, 0.99],
              [0.01, 1.0], [0.01, 1.0],  # c_lf
              [0.001, 1.0]]
    problem = {
        'num_vars': len(bounds),
        'names': ['z3_factor',
                  'cZ0Z1', 'cZ',
                  'c_adr'],
        'bounds': scale_morris(bounds)
    }
else:
    bounds = [[0.75, 0.99],  # z3_factor
              [0.01, 1.0], [0.01, 1.0],  # c_lf
              [0.001, 1.0],  # cadr
              [1100.0, 3650.0],  # k_g
              [0.0001, 1.0], [0, 1.0],  # gamma
              [0.0001, 1.0],  # f_transp
              [0.34, 0.41], [0.34, 0.41],  # FC
              [0.44, 0.61], [0.47, 0.53],  # SAT
              [0.13, 0.24], [0.13, 0.24],  # WP
              [0.39, 0.42], [0.34, 0.43], [0.35, 0.41],  # W100
              [0.01, 0.05], [110.0, 369.0],  # f_oc, k_oc
              [0.01, 1.0],  # beta_runoff
              [140.0, 7000.0],  # age_rate
              [130.0, 230.0],  # dt_50_ab
              [26.0, 37.0],  # dt_50_ref
              [1.0, 3.0],  # epsilon (in absolute, convert to negative!!)
              [0.01, 1.0]]  # beta_moisture

    names = ['z3_factor',
             'cZ0Z1', 'cZ',
             'c_adr',
             'k_g',
             'gamma01', 'gammaZ',
             'f_transp',
             'FCz2', 'FCZ',
             'SATz2', 'SATZ',
             'WPz2', 'WPZ',
             'W100_30mm', 'W100_1550mm', 'W100_1555mm',
             'f_oc', 'k_oc',
             'beta_runoff',
             'dt_50_aged',
             'dt_50_ab',
             'dt_50_ref',
             'epsilon_iso',
             'beta_moisture'
             ]

    problem = {
        'num_vars': len(bounds),
        'names': names,
        'bounds': scale_morris(bounds)
    }


def get_param_values(p=4.0, grid_jump=2, r=4):
    # r : Trajectories
    param_values = mos.sample(problem, r, num_levels=p, grid_jump=grid_jump)
    saveInputMatrix(param_values)

    return param_values


def get_runs(params):
    runs = int(params.shape[0])  # runs = rk+r = r(k+1)
    print("Runs: ", int(runs))
    return runs


def getInputVector(row, sample_matrix):
    """
    :param row: relevant sample row
    :param sample_matrix: numpy sample matrix
    :return: a numpy row with required input parameters
    """
    test_vector = sample_matrix[row]
    return test_vector
