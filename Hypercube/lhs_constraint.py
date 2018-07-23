
"""
Move chunks after import statements to phd-model-master
"""


def scale(bounds):
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


bounds = [[0.75, 0.99],  # z3_factor
          [0.01, 1.0], [0.01, 1.0],  # c_lf
          [0.001, 1.0],  # cadr
          [1100.0, 3650.0],  # k_g
          [0.0001, 1.0], [0, 1.0],  # gamma
          [0.0001, 1.0],  # f_transp
          [0.01, 0.05], [7.0, 16180],  # f_oc, k_oc
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
    'bounds': scale(bounds),
    'upper': get_upper(bounds)
}


# from SALib.sample import latin
from SALib.sample import latin
from copy import deepcopy
import numpy as np
from random import randint

# Define models to run
# problem = get_problem()
names = problem['names']
samples = 5
upper = problem['upper']


def get_constrained_matrix(smps):
    tot_len = 0
    while tot_len < 1:
        smps_values = latin.sample(problem, smps)
        new_values = deepcopy(smps_values)
        for row in range(len(smps_values))[::-1]:
            v = smps_values[row]
            if v[1] < v[2] or v[5] < v[6]:
                # print("Deleting: ", new_values[row])
                new_values = np.delete(new_values, row, axis=0)
        tot_len = len(new_values)
    return new_values

test_values = get_constrained_matrix(samples)

print("Start")
print(len(test_values))
print(test_values)

while len(test_values) < samples:
    resample = get_constrained_matrix(samples)
    test_values = np.vstack([test_values, resample[randint(0, len(resample)-1)]])

print("End")
print(len(test_values))
print(test_values)