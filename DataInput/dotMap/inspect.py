from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

landuse = readmap("landuse")
# aguila(landuse)
clone = readmap("clone")

create_outlet = False
if create_outlet:
    out_burn = readmap("dem_ldd_burn2")

    out_ditch_ldd = lddcreate(out_burn, 1E35, 1E35, 1E35, 1E35)
    out_ditch = pit(out_ditch_ldd)
    # Save
    report(out_ditch, 'outlet_multi') # 0 to 67 nominal outlets

    out_true = readmap('outlet_true')
    # aguila(out_ditch, out_true)
    multi_out_bool = boolean(out_ditch)  # True False outlets
    multi_out_nom = ifthenelse(multi_out_bool, nominal(1), nominal(0)) # 0 to 1 nominal outlets
    aguila(out_ditch, multi_out_bool, multi_out_nom)
    # Save
    report(multi_out_bool, 'out_multi_bool')
    report(multi_out_nom, 'out_multi_nom')

# Test for radians vs degrees
# aguila(slope_rad)
# angles = atan(max(slope(dem), 0.001))  # Results in a range of 0 - 360 (degrees)
# aguila(angles)
# slope_deg = slope_rad/0.0174533  # 1 deg = 0.0174533 rad
# aguila(slope_deg)


create_weekly = True
if create_weekly:
    weekly = readmap("weekly_smp")
    weekly = order(weekly)
    result = boolean(weekly)
    # report(weekly, 'weekly_ord.map')  # stores a ".map" file

    word = readmap("weekly_ord")
    north = ifthen(word < 31, word)
    north = order(north)
    # Talweg and South
    ts = ifthen(word > 30, word)
    south = ifthen(ts > 55, ts)
    south = order(south)
    valley = ifthen(ts < 56, ts)
    valley = order(valley)
    # valley -= 30
    # aguila(valley, landuse)
    report(north, 'north_wk.map')
    report(south, 'south_wk.map')
    report(valley, 'valley_wk.map')
    # detail = readmap("detailed_smp")


# ldd1 = readmap("ldd")
# Surface
create_ldd = False
if create_ldd:
    out_burn = readmap("dem_ldd_burn2")
    out_ditch_ldd = lddcreate(out_burn, 1E35, 1E35, 1E35, 1E35)
    # report(out_ditch_ldd, 'ldd_subs')

    # Burnt for Runoff routing:
    dem = readmap("dem_slope")
    dem_ldd = readmap("dem_ldd")
    slope_rad = sin(atan(max(slope(dem), 0.001)))
    ldd1 = lddcreate(dem_ldd, 1e31, 1e31, 1e31, 1e31)  # second param = "outflowdepth" ??
    # Subsurface
    ldd2 = lddcreate(dem, 1e31, 1e31, 1e31, 1e31)  # Now using the burn2_ldd


create_apps = False
if create_apps:
    dem = readmap("dem_slope")
    mask = dem / dem
    out_burn = readmap("dem_ldd_burn2")

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
    fa_cr = readmap("farmcrop")  # Contains codes to assign appropriate dosage
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

    # app2 = ifthenelse(fa_cr == 1511, 1*app_conc,  # 1511 (Burger-Beet)
    #                                ifthenelse(fa_cr == 1611, 1*app_conc,  # 1611 (Kopp-Beet),
    #                                           0*app_conc))
    # app3 = ifthenelse(fa_cr == 1112, 1 * app_conc,  # 1112 (Friess-Corn)
    #                                ifthenelse(fa_cr == 1212, 1 * app_conc,  # 1212 (Speich-Corn),
    #                                           ifthenelse(fa_cr == 1412, 1 * app_conc,  # 1412 (Schmitt-Corn),
    #                                                      ifthenelse(fa_cr == 1312, 1 * app_conc,  # 1312 (Mahler-Corn)
    #                                                                 0 * app_conc))))
    # app2ug = app2 * cellarea()
    # app2delta = ifthenelse(app2 > 0, scalar(-32.3), scalar(-23.7))
    new = defined(out_burn)
    burn_farmcrop = ifthen(new, fa_cr)
    # report(burn_farmcrop, 'crop_burn.map')
    aguila(fa_cr, app_conc, out_burn, new, burn_farmcrop)


