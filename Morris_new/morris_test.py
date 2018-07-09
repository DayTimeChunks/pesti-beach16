# -*- coding: utf-8 -*-
from SALib.sample import morris as mos
from morris_analysis import *


def saveMorris(Names, SiList, filename="MorrisResults.json"):
    import json
    mListSi = {}
    for i in range(len(Names)):
        var_name = Names[i]
        Si = SiList[i]
        mSi = {}
        mSi["mu"] = Si["mu"].tolist()
        mSi["mu_star"] = Si["mu_star"].tolist()
        mSi['mu_star_conf'] = Si['mu_star_conf']
        mSi['sigma'] = Si['sigma'].tolist()
        mSi['names'] = Si['names']
        mListSi[var_name] = mSi

    with open(filename, 'w') as f:
        json.dump(mListSi, f)


problem = {
    'num_vars': 24,
    'names': ['z3_factor',
              'c30', 'c420',
              'c_adr',
              'k_g',
              'gamma0', 'gamma1', 'gammaZ',
              'f_transp',
              'FCz2', 'FCZ',
              'SATz2', 'SATZ',
              'W100_30mm', 'W100_1550mm', 'W100_1555mm',
              'f_oc', 'k_oc',
              'beta_runoff',
              'age_rate',
              'dt_50_ab',
              'dt_50_ref',
              'epsilon_iso',
              'beta_moisture'
              ],
    'bounds': [[0.75, 0.99],  # z3_factor
               [0.01, 1.0], [0.01, 1.0],  # c_lf
               [0.001, 1.0],  # cadr
               [1100, 3650],  # k_g
               [0, 1.0], [0, 1.0], [0, 1.0],  # gamma
               [0, 1.0],  # f_transp
               [0.34, 0.41], [0.34, 0.41],  # FC
               [0.44, 0.61], [0.47, 0.53],  # SAT
               [0.39, 0.42], [0.34, 0.43], [0.35, 0.41],  # W
               [0.01, 0.05], [110, 369],  # f_oc, k_oc
               [0.01, 1],  # beta_runoff
               [140, 7000],  # age_rate
               [130, 230],  # dt_50_ab
               [26, 37],  # dt_50_ref
               [-3.0, -1.0],  # epsilon
               [0, 1]]  # beta_moisture
}


global param_values
p = 4.0
grid_jump = 2
# r = 10  # Trajectories
r = 4

param_values = mos.sample(problem, r, num_levels=p, grid_jump=grid_jump)
saveInputMatrix(param_values)

runs = int(param_values.shape[0])  # runs = rk+r = r(k+1)
print("Runs: ", int(runs))