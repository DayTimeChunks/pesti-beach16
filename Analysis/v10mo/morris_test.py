# -*- coding: utf-8 -*-
from SALib.sample import morris as mos

from calib_v1.mVar_v1.morris.morris_analysis import *


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
    'num_vars': 13,
    'names': ['thsat0', 'thsat2',
              'thfc0', 'thfc2',
              'thwp',
              'gamma0', 'gamma1',
              's0', 's1',
              'clf0', 'clf1', 'clf2',
              'cadr1'
              ],
    'bounds': [[0.549, 0.61], [0.549, 0.61],  # sat
               [0.360, 0.40], [0.360, 0.40],  # fc
               [0.09,  0.15],  # wp
               [0.001, 1], [0.001, 1],  # gamma
               [0.004852, 1], [0.004852, 1],  # s
               [0.001, 1], [0.001, 1], [0.001, 1],  # clf
               [0.001, 1]  # cadr
               ]
}

# problem = {
#     'num_vars': 3,
#     'names': ['thsat0', 'thsat2',
#               'thfc0'
#               ],
#     'bounds': [[0.549, 0.61], [0.549, 0.61],  # sat
#                [0.360, 0.40]
#                ]
# }

global param_values
p = 4.0
grid_jump = 2
# r = 10  # Trajectories
r = 4

param_values = mos.sample(problem, r, num_levels=p, grid_jump=grid_jump)
saveInputMatrix(param_values)

runs = int(param_values.shape[0])  # runs = rk+r = r(k+2)
print("Runs: ", int(runs))