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

# ldd1 = readmap("ldd")
# Surface
ldd1 = lddcreate(dem_ldd, 1e31, 1e31, 1e31, 1e31)  # second param = "outflowdepth" ??

# Subsurface
ldd2 = lddcreate(dem, 1e31, 1e31, 1e31, 1e31)

cell_area = cellarea()
up_area = accuflux(ldd1, cell_area)
up_area2 = accuflux(ldd2, cell_area)
wetness = ln(up_area / tan(slope_rad))
wetness2 = ln(up_area2 / tan(slope_rad))


wet = wetness - wetness2
mask = dem/dem
detail = readmap("detailed_smp")
# Applications Mass
# Product concentration (active ing.)
double = 2.0  # ~ Dosage for corn when growing beet
d_gold = 915 * 10 ** 6  # ug/L S-met
m_gold = 960 * 10 ** 6  # ug/L

# Dosages # L/Ha * 1Ha/1000m2 = L/m2
d_beet = None
d_corn = 2.1 * 1 / 10 ** 4  # L/Ha * 1 Ha / 10000 m2
m_beet = 0.6 * 1 / 10 ** 4 * double
m_corn = 2.0 * 1 / 10 ** 4
m_beet_Friess = 0.6 * 1 / 10 ** 4 * double  # (Likely larger dosage, early in the season)
m_beet_Mathis = 0.6 * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)

# Assign dosages based on Farmer-Crop combinations [g/m2]
fa_cr = readmap("farmer_crop")  # Contains codes to assign appropriate dosage
app_conc = ifthenelse(fa_cr == 1111, m_beet_Friess * m_gold * mask,  # 1111 (Friess, Beet)
                      # 1112 (Friess-Corn),
                      ifthenelse(fa_cr == 1122, m_corn * m_gold * mask,
                                 # 1212 (Speich-Corn),
                                 ifthenelse(fa_cr == 1212, m_corn * m_gold * mask,
                                            # 1312 (Mahler-Corn),
                                            ifthenelse(fa_cr == 1312, m_corn * m_gold * mask,
                                                       # 1412 (Schmitt-Corn)
                                                       ifthenelse(fa_cr == 1412, d_corn * d_gold * mask,
                                                                  # 1511 (Burger-Beet)
                                                                  ifthenelse(fa_cr == 1511,
                                                                             m_beet * m_gold * mask,
                                                                             # 1711 (Mathis-Beet),
                                                                             ifthenelse(fa_cr == 1711,
                                                                                        m_beet_Mathis * m_gold * mask,
                                                                                        # 1612 (Kopp-Beet)
                                                                                        ifthenelse(
                                                                                            fa_cr == 1611,
                                                                                            m_beet * m_gold * mask,
                                                                                            0 * mask))))))))

app2 = ifthenelse(fa_cr == 1511, 1*app_conc,  # 1511 (Burger-Beet)
                               ifthenelse(fa_cr == 1611, 1*app_conc,  # 1611 (Kopp-Beet),
                                          0*app_conc))
app3 = ifthenelse(fa_cr == 1112, 1 * app_conc,  # 1112 (Friess-Corn)
                               ifthenelse(fa_cr == 1212, 1 * app_conc,  # 1212 (Speich-Corn),
                                          ifthenelse(fa_cr == 1412, 1 * app_conc,  # 1412 (Schmitt-Corn),
                                                     ifthenelse(fa_cr == 1312, 1 * app_conc,  # 1312 (Mahler-Corn)
                                                                0 * app_conc))))
app2ug = app2 * cellarea()
app2delta = ifthenelse(app2 > 0, scalar(-32.3), scalar(-23.7))
# aguila(fa_cr, app2, app3)

# layon = readmap("landuseCopy")
aguila(landuse)
