
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
m_beet_Friess = 0.6 * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)
m_beet_Mathis = 0.6 * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)
mask = dem/dem
# Assign dosages based on Farmer-Crop combinations [ug/m2]
fa_cr = readmap("farmer_crop")  # Contains codes to assign appropriate dosage
app_conc = ( ifthenelse(fa_cr == 1111,  # 1111 (Friess, Beet)
               m_beet_Friess * m_gold * mask,
               ifthenelse(fa_cr == 1122,  # 1112 (Friess-Corn),
                          m_corn * m_gold * mask,
                          ifthenelse(fa_cr == 1212,  # 1212 (Speich-Corn),
                                     m_corn * m_gold * mask,
                                     ifthenelse(fa_cr == 1312,  # 1312 (Mahler-Corn),
                                                m_corn * m_gold * mask,
                                                ifthenelse(fa_cr == 1412,  # 1412 (Schmitt-Corn)
                                                           d_corn * d_gold * mask,
                                                           ifthenelse(fa_cr == 1511,  # 1511 (Burger-Beet)
                                                                      m_beet * m_gold * mask,
                                                                      # 1711 (Mathis-Beet),
                                                                      ifthenelse(fa_cr == 1711,
                                                                                 m_beet_Mathis * m_gold * mask,
                                                                                 # 1611 (Kopp-Beet)
                                                                                 ifthenelse(
                                                                                     fa_cr == 1611,
                                                                                     m_beet * m_gold * mask,
                                                                                     0 * mask)))))))))
app1 = ifthenelse(fa_cr == 1111, 1 * app_conc,
                               # 1111 (Friess, Beet), 1112 (Friess-Corn),
                               ifthenelse(fa_cr == 1112, 1 * app_conc,
                                          ifthenelse(fa_cr == 1711, 1 * app_conc,  # 1711 (Mathis-Beet)
                                                     0 * app_conc)))
# Pesticide applied (ug/m2) on Julian day 197 (April 14, 2016).
# April 14, Kopp and Burger
app2 = ifthenelse(fa_cr == 1511, 1 * app_conc,  # 1511 (Burger-Beet)
                       ifthenelse(fa_cr == 1611, 1 * app_conc,  # 1611 (Kopp-Beet),
                                  0 * app_conc))

# Pesticide applied (ug/m2) on Julian day 238 (May 25, 2016).
# May 25, Schmidt and Speich, and (out of transect): Friess and Mahler
# Note: Speich and Friess could be 1 week later.
app3 = ifthenelse(fa_cr == 1112, 1 * app_conc,  # 1112 (Friess-Corn)
                       ifthenelse(fa_cr == 1212, 1 * app_conc,  # 1212 (Speich-Corn),
                                  ifthenelse(fa_cr == 1412, 1 * app_conc,  # 1412 (Schmitt-Corn),
                                             ifthenelse(fa_cr == 1312, 1 * app_conc,  # 1312 (Mahler-Corn)
                                                        0 * app_conc))))
        
aguila(app1, app2, app3)
# Landscape analysis
# z0_tss = TimeoutputTimeseries("res_z0_mm", self, "outlet.map", noHeader=False)
# z1_tss = TimeoutputTimeseries("res_z1_mm", self, "outlet.map", noHeader=False)
# z2_tss = TimeoutputTimeseries("res_z2_mm", self, "outlet.map", noHeader=False)
# ztot_tss = TimeoutputTimeseries("res_ztot_mm", self, "outlet.map", noHeader=False)
# z0_tss.sample(z0)
# z1_tss.sample(z1)
# z2_tss.sample(z2)
# ztot_tss.sample(tot_depth)