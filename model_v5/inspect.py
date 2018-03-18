
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
res1 = readmap("ini_theta_z1")
res2 = readmap("resV2274_theta_z2")

