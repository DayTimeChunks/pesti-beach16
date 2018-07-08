from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

landuse = readmap("landuse2016")
# aguila(landuse)
clone = readmap("clone")


generate_maps = {
    "outlets": False,
    "weekly_soils": False,
    "applications": False,
    "create_ldd": False
}

create_outlet = generate_maps['outlets']
if create_outlet:
    # Burned to account for a connected ditch
    out_burn = readmap("dem_ldd_burn3")

    # Build the ldd that generates the outlet(s)
    out_ditch_ldd = lddcreate(out_burn, 1E35, 1E35, 1E35, 1E35)
    out_ditch = pit(out_ditch_ldd)

    # Get true outlet (to sum all outlet outputs)
    out_true = ifthenelse(out_ditch == 33, nominal(1), nominal(0))

    # Convert to nominal 0's and 2's for all outlets
    multi_out_bool = boolean(out_ditch)  # True False outlets
    multi_out_nom = ifthenelse(multi_out_bool, nominal(1), nominal(0))

    # View
    aguila(out_ditch, out_true, multi_out_nom)

    # Save
    report(out_ditch, 'outlet_multi_v3.map')  # 0 to 67 nominal outlets
    report(out_true, 'outlet_v3.map')
    report(multi_out_nom, 'out_multi_nom_v3.map')

# ldd1 = readmap("ldd")
# Surface
create_ldd = generate_maps['create_ldd']
if create_ldd:
    out_burn = readmap("dem_ldd_burn3")

    # Burnt for Runoff routing:
    # dem = readmap("dem_slope")
    # dem_ldd = readmap("dem_ldd")
    # slope_rad = sin(atan(max(slope(dem), 0.001)))
    # ldd1 = lddcreate(dem_ldd, 1e31, 1e31, 1e31, 1e31)

    # Subsurface
    ldd_subs = lddcreate(out_burn, 1E35, 1E35, 1E35, 1E35)  # Now using the burn3_ldd
    aguila(ldd_subs)
    # report(ldd_subs, 'ldd_subs_v3.map')

# Test for radians vs degrees
# aguila(slope_rad)
# angles = atan(max(slope(dem), 0.001))  # Results in a range of 0 - 360 (degrees)
# aguila(angles)
# slope_deg = slope_rad/0.0174533  # 2 deg = 0.0174533 rad
# aguila(slope_deg)


create_weekly = generate_maps['weekly_soils']
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
    # report(north, 'north_wk.map')
    # report(south, 'south_wk.map')
    # report(valley, 'valley_wk.map')
    # detail = readmap("detailed_smp")
    south_ave = ifthen(south == 1, south)
    valley_ave = ifthen(valley == 1, valley)
    north_ave = ifthen(north == 1, north)
    # report(order(north_ave), 'north_ave.map')
    # report(order(valley_ave), 'valley_ave.map')
    # report(order(south_ave), 'south_ave.map')
    north_nom = nominal(boolean(north))
    valley_nom = nominal(boolean(valley))
    south_nom = nominal(boolean(south))
    # report(north_nom, 'north_nom.map')
    # report(valley_nom, 'valley_nom.map')
    # report(south_nom, 'south_nom.map')
    # aguila(south_nom, south)

detailed_s11 = False
if detailed_s11:
    detail_ord = order(readmap("detailed_smp"))
    ts = ifthen(detail_ord > 34, detail_ord)
    ts = order(ts)
    # S
    sou = ifthen(ts > 30, ts)
    s11 = ifthen(sou < 39, sou)
    s11_ord = order(s11)
    report(s11_ord, 's11_ord.map')

    # report(s11_out, "s11_out.map")

detailed = False
if detailed:
    detail_ord = order(readmap("detailed_smp"))
    nor = ifthen(detail_ord < 35, detail_ord)
    nor = order(nor)

    ts = ifthen(detail_ord > 34, detail_ord)
    ts = order(ts)
    val = ifthen(ts < 31, ts)
    # val = order(val)

    # n1
    n1 = ifthen(nor < 5, nor)
    n1_out = ifthen(n1 == 4, n1)
    n1 = nominal(boolean(n1))
    n1_out = nominal(boolean(n1_out))

    # n2
    nor = ifthen(nor > 4, nor)
    n2 = ifthen(nor < 11, nor)
    n2_out = ifthen(n2 == 10, n2)
    n2 = nominal(boolean(n2))
    n2_out = nominal(boolean(n2_out))
    # n3
    nor = ifthen(nor > 10, nor)
    n3 = ifthen(nor < 17, nor)
    n3_out = ifthen(n3 == 16, n3)
    n3 = nominal(boolean(n3))
    n3_out = nominal(boolean(n3_out))

    # n4
    nor = ifthen(nor > 16, nor)
    n4 = ifthen(nor < 21, nor)
    n4_out = ifthen(n4 == 20, n4)
    n4 = nominal(boolean(n4))
    n4_out = nominal(boolean(n4_out))
    # n5
    nor = ifthen(nor > 20, nor)
    n5 = ifthen(nor < 25, nor)
    n5_out = ifthen(n5 == 24, n5)
    n5 = nominal(boolean(n5))
    n5_out = nominal(boolean(n5_out))
    # n7
    nor = ifthen(nor > 24, nor)
    n7 = ifthen(nor < 30, nor)
    n7_out = ifthen(n7 == 29, n7)
    n7 = nominal(boolean(n7))
    n7_out = nominal(boolean(n7_out))
    # n8
    nor = ifthen(nor > 29, nor)
    n8 = ifthen(nor < 35, nor)
    n8_out = ifthen(n8 == 34, n8)
    n8 = nominal(boolean(n8))
    n8_out = nominal(boolean(n8_out))
    report(n1, "n1_nom.map")
    report(n2, "n2_nom.map")
    report(n3, "n3_nom.map")
    report(n4, "n4_nom.map")
    report(n5, "n5_nom.map")
    report(n7, "n7_nom.map")
    report(n8, "n8_nom.map")
    report(n1_out, "n1_out.map")
    report(n2_out, "n2_out.map")
    report(n3_out, "n3_out.map")
    report(n4_out, "n4_out.map")
    report(n5_out, "n5_out.map")
    report(n7_out, "n7_out.map")
    report(n8_out, "n8_out.map")

    # S
    sou = ifthen(ts > 30, ts)
    s11 = ifthen(sou < 39, sou)
    s11_out = ifthen(s11 == 38, s11)
    s11 = nominal(boolean(s11))
    s11_out = nominal(boolean(s11_out))
    report(s11, "s11_nom.map")
    report(s11_out, "s11_out.map")

    sou = ifthen(sou > 38, sou)
    sou = order(sou)
    s12 = ifthen(sou < 8, sou)
    s12_out = ifthen(s12 == 7, s12)
    s13 = ifthen(sou > 7, sou)
    s13_out = ifthen(s13 == 8, s13)
    s12 = nominal(boolean(s12))
    s12_out = nominal(boolean(s12_out))
    s13 = nominal(boolean(s13))
    s13_out = nominal(boolean(s13_out))
    report(s12, "s12_nom.map")
    report(s13, "s13_nom.map")
    report(s12_out, "s12_out.map")
    report(s13_out, "s13_out.map")

    # T
    # Problem with ordering,
    # separate into west and east first
    val = ifthen(ts < 31, ts)
    val = order(val)
    west = ifthen(val > 16, val)

    t5 = ifthen(west < 21, west)
    t5_out = ifthen(t5 == 20, t5)
    t5 = nominal(boolean(t5))
    t5_out = nominal(boolean(t5_out))
    west = ifthen(west > 20, west)
    t7 = ifthen(west < 26, west)
    t7_out = ifthen(t7 == 25, t7)
    t7 = nominal(boolean(t7))
    t7_out = nominal(boolean(t7_out))

    t8 = ifthen(west > 25, west)
    t8_out = ifthen(t8 == 26, t8)
    t8 = nominal(boolean(t8))
    t8_out = nominal(boolean(t8_out))

    report(t5, "v5_nom.map")
    report(t7, "v7_nom.map")
    report(t8, "v8_nom.map")
    report(t5_out, "v5_out.map")
    report(t7_out, "v7_out.map")
    report(t8_out, "v8_out.map")

    east = ifthen(val < 17, val)
    east = order(east)
    east = ifthen(east != 13, east)
    #
    t10 = ifthen(east < 6, east)
    t10_out = ifthen(t10 == 5, t10)
    east = ifthen(east > 5, east)
    t9 = ifthen(east < 11, east)
    t9_out = ifthen(t9 == 10, t9)
    t4 = ifthen(east > 11, east)
    t4_out = ifthen(t4 == 12, t4)

    t10 = nominal(boolean(t10))
    t10_out = nominal(boolean(t10_out))
    t9 = nominal(boolean(t9))
    t9_out = nominal(boolean(t9_out))
    t4 = nominal(boolean(t4))
    t4_out = nominal(boolean(t4_out))
    report(t4, "v4_nom.map")
    report(t9, "v9_nom.map")
    report(t10, "v10_nom.map")
    report(t4_out, "v4_out.map")
    report(t9_out, "v9_out.map")
    report(t10_out, "v10_out.map")

    # report(valley_nom, 'valley_nom.map')
    # report(south_nom, 'south_nom.map')
    aguila(t8, t8_out, landuse)


# Have not run it yet!
create_apps = generate_maps['applications']
if create_apps:
    dem = readmap("dem_slope")
    mask = dem / dem
    out_burn = readmap("dem_ldd_burn3")

    # Applications Mass
    # Product concentration (active ing.)
    double = 2.0  # ~ Dosage for corn when growing beet
    d_gold = 915 * 10 ** 6  # ug/L S-met
    m_gold = 960 * 10 ** 6  # ug/L

    # Dosages # L/Ha * 1Ha/1000m2 = L/m2
    d_beet = None
    d_corn = 2.1 * 1 / 10 ** 4  # L/Ha * 2 Ha / 10000 m2
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

    # app2 = ifthenelse(fa_cr == 1511, 2*app_conc,  # 1511 (Burger-Beet)
    #                                ifthenelse(fa_cr == 1611, 2*app_conc,  # 1611 (Kopp-Beet),
    #                                           0*app_conc))
    # app3 = ifthenelse(fa_cr == 1112, 2 * app_conc,  # 1112 (Friess-Corn)
    #                                ifthenelse(fa_cr == 1212, 2 * app_conc,  # 1212 (Speich-Corn),
    #                                           ifthenelse(fa_cr == 1412, 2 * app_conc,  # 1412 (Schmitt-Corn),
    #                                                      ifthenelse(fa_cr == 1312, 2 * app_conc,  # 1312 (Mahler-Corn)
    #                                                                 0 * app_conc))))
    # app2ug = app2 * cellarea()
    # app2delta = ifthenelse(app2 > 0, scalar(-32.3), scalar(-23.7))
    new = defined(out_burn)
    burn_farmcrop = ifthen(new, fa_cr)
    # report(burn_farmcrop, 'crop_burn.map')
    # old_farm = readmap("farmcrop")
    report(burn_farmcrop, 'farm_burn_v3.map')


pcodes = readmap('plot_code16')
aguila(pcodes)