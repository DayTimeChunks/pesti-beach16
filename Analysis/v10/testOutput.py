# -*- coding: utf-8 -*-
import numpy as np
runs = 2
variable = "resM_global_mb_pest"

Y = np.zeros([runs])

for i in range(runs):
    folder = i+1
    path = str(folder) + "/" + variable + ".tss"
    # res = pd.read_table(path, header=None)
    res = np.loadtxt(path, float, skiprows=4, usecols=[1])
    Y[i] = res[-1]

print (Y)