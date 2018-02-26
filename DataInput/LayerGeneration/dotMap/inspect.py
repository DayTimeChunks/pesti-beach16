from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

dem = readmap("dem_slope")
dem_ldd = readmap("dem_ldd")
slope_rad = sin(atan(max(slope(dem), 0.001)))

# Test for radians vs degrees
# aguila(slope_rad)
# angles = atan(max(slope(dem), 0.001))  # Results in a range of 0 - 360 (degrees)
# aguila(angles)
# slope_deg = slope_rad/0.0174533  # 1 deg = 0.0174533 rad
# aguila(slope_deg)

landuse = readmap("landuse")
clone = readmap("clone")

ini_theta = dem - dem + 0.20
tot_depth = (dem - mapminimum(dem))*scalar(10**3)

ldd1 = readmap("ldd")
ldd2 = lddcreate(dem, 10, 1e31, 1e31, 1e31)  # second param = "outflowdepth" ??

cell_area = cellarea()
up_area = accuflux(ldd1, cell_area)
up_area2 = accuflux(ldd2, cell_area)
wetness = ln(up_area / tan(slope_rad))
wetness2 = ln(up_area2 / tan(slope_rad))

wet = wetness - wetness2

aguila(wet)
