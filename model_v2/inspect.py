
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

dem = readmap("dem_slope")  # 192 - 231 m a.s.l
dem_route = readmap("dem_ldd")  # To route surface run-off
zero_map = dem - dem
datum_depth = (dem - mapminimum(dem)) * scalar(10 ** 3)  # mm
z0 = zero_map + 10  # mm
z1 = zero_map + 140  # mm
z2 = datum_depth + 300 - z0 - z1  # mm (150mm at outlet)
tot_depth = z0 + z1 + z2

aguila(tot_depth, z2, res2)
# Landscape analysis
# z0_tss = TimeoutputTimeseries("res_z0_mm", self, "outlet.map", noHeader=False)
# z1_tss = TimeoutputTimeseries("res_z1_mm", self, "outlet.map", noHeader=False)
# z2_tss = TimeoutputTimeseries("res_z2_mm", self, "outlet.map", noHeader=False)
# ztot_tss = TimeoutputTimeseries("res_ztot_mm", self, "outlet.map", noHeader=False)
# z0_tss.sample(z0)
# z1_tss.sample(z1)
# z2_tss.sample(z2)
# ztot_tss.sample(tot_depth)