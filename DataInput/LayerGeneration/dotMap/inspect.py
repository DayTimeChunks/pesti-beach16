from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

dem = readmap("dem_slope")
slope_rad = sin(atan(max(slope(dem), 0.001)))

# Test for radians vs degrees
# aguila(slope_rad)
# angles = atan(max(slope(dem), 0.001))  # Results in a range of 0 - 360 (degrees)
# aguila(angles)
# slope_deg = slope_rad/0.0174533  # 1 deg = 0.0174533 rad
# aguila(slope_deg)

landuse = readmap("landuse")
clone = readmap("clone")
aguila(clone)
