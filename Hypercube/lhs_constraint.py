
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
import numpy as np

# Define models to run
# problem = get_problem()
names = problem['names']
samples = 2
upper = problem['upper']
param_values = latin.sample(problem, samples)
new_values = latin.sample(problem, samples)

updated_matrix = np.vstack([param_values, new_values[0]])


print(param_values.shape)
cprint(updated_matrix.shape)