
from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

# Correction of
# fields = readmap("fields_cover")
# out = readmap("outlet_true")
# res = readmap("resdt200_accu_runoff_m3.map")
# res = readmap("resdt200_cell_runoff_m3.map")
# aguila(res, fields, out)

res = readmap("ini_theta_z0")
res2 = readmap("6/z0theta0.279")

net = res - res2
aguila(res2)
