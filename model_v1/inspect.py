
from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

# Correction of
fields = readmap("fields_cover")
out = readmap("outlet_true")
res = readmap("resdt200_accu_runoff_m3.map")
res = readmap("resdt200_cell_runoff_m3.map")
aguila(res, fields, out)

# Landscape analysis
# self.z0_tss = TimeoutputTimeseries("res_z0_mm", self, "outlet.map", noHeader=False)
# self.z1_tss = TimeoutputTimeseries("res_z1_mm", self, "outlet.map", noHeader=False)
# self.z2_tss = TimeoutputTimeseries("res_z2_mm", self, "outlet.map", noHeader=False)
# self.ztot_tss = TimeoutputTimeseries("res_ztot_mm", self, "outlet.map", noHeader=False)
# self.z0_tss.sample(self.z0)
# self.z1_tss.sample(self.z1)
# self.z2_tss.sample(self.z2)
# self.ztot_tss.sample(self.tot_depth)