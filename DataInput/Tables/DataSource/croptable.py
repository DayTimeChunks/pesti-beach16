# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

# Which computer's directory?
PC = True
# What model version?
version = 'v4'

if PC:
    path = "D:/Documents/these_pablo/Models/BEACH2016/DataInput/Tables/DataSource/"
else:
    # path = "C:/Users/DayTimeChunks/Documents/PhD/HydrologicalMonitoring"
    path = "C:/Users/DayTimeChunks/Documents/Models/pesti-beach16/DataInput/Tables/DataSource/"

path += "croptable_start.csv"

croptable = pd.read_csv(path, sep=";")

# For landuse reference
landuse = {"Corn": 1, "Wheat": 2, "Oats": 3, "Alfalfa": 4, "Beet": 5,
           "Greenery": 6, "Dirt Road": 7, "Grass Road": 8, "Paved Road": 9, "Ditch": 10,
           "Fallow": 11, "Hedge": 12, "Orchard": 13}

# Corn, Wheat, Beet
crop_conditions = [
    (croptable.loc[1:, 'crop_type'] == 1),  # Corn
    (croptable.loc[1:, 'crop_type'] == 2),  # Wheat
    (croptable.loc[1:, 'crop_type'] == 5)  # Beet
]
sow_yy = [2016, 2015, 2016]
sow_mm = [4, 10, 3]  # Apr, Oct, March
sow_dd = [10, 15, 25]

# Corn, Wheat, Beet
len_grow_stage_ini = [28, 160, 50]  # (days)
len_dev_stage = [14, 75, 40]  # (days)
len_mid_stage = [14, 75, 50]  # (days)
len_end_stage = [14, 25, 40]  # (days)
kcb_ini = [0.3, 0.7, 0.35]  # (-)
kcb_mid = [1.2, 1.15, 1.2]  # (-)
kcb_end = [0.5, 0.25, 0.7]  # (-)
max_LAI = [7, 6.3, 4.5]  # (-)
mu = [3.5, 1.5, 1.8]  # (g biomass / MJ), not used, see: getBiomassCover()
max_height = [2.43, 0.79, 0.55]  # m  (SWAT needs 'm')
max_root_depth = [1, 1.4, 1]  # m (convert to mm on run())
p_tab = [0.55, 0.55, 0.55]  # Depletion coeff (-) assumed total from FAO.

croptable.loc[1:, 'sow_yy'] = np.select(crop_conditions, sow_yy, default=0)
croptable.loc[1:, 'sow_mm'] = np.select(crop_conditions, sow_mm, default=0)
croptable.loc[1:, 'sow_dd'] = np.select(crop_conditions, sow_dd, default=0)
croptable.loc[1:, 'len_grow_stage_ini'] = np.select(crop_conditions, len_grow_stage_ini, default=0)
croptable.loc[1:, 'len_dev_stage'] = np.select(crop_conditions, len_dev_stage, default=0)
croptable.loc[1:, 'len_mid_stage'] = np.select(crop_conditions, len_mid_stage, default=0)
croptable.loc[1:, 'len_end_stage'] = np.select(crop_conditions, len_end_stage, default=0)
croptable.loc[1:, 'kcb_ini'] = np.select(crop_conditions, kcb_ini, default=0.15)
croptable.loc[1:, 'kcb_mid'] = np.select(crop_conditions, kcb_mid, default=0.15)
croptable.loc[1:, 'kcb_end'] = np.select(crop_conditions, kcb_end, default=0.15)
croptable.loc[1:, 'max_LAI'] = np.select(crop_conditions, max_LAI, default=0)
croptable.loc[1:, 'mu'] = np.select(crop_conditions, mu, default=0)
croptable.loc[1:, 'max_height'] = np.select(crop_conditions, max_height, default=0)
croptable.loc[1:, 'max_root_depth'] = np.select(crop_conditions, max_root_depth, default=0)
croptable.loc[1:, 'p_tab'] = np.select(crop_conditions, p_tab, default=0)

theta_condition = [
    (croptable.loc[1:, 'crop_type'] == 10)  # Ditch
]
croptable.loc[1:, 'theta_sat_z0z1'] = np.select(theta_condition, [0.61], default=0.61)
croptable.loc[1:, 'theta_fcap_z0z1'] = 0.2
croptable.loc[1:, 'theta_sat_z2'] = np.select(theta_condition, [0.61], default=0.61)
croptable.loc[1:, 'theta_fcap_z2'] = 0.2
croptable.loc[1:, 'theta_wp'] = 0.1

ksat_condition = [
    (croptable.loc[1:, 'crop_type'] == 7),  # Dirt Road, assumed Group C, 1.3-3.8 mm/h -> 2.6 mm/h = 62mm/day
    (croptable.loc[1:, 'crop_type'] == 9),  # Paved Road, assumed Group D, < 1.3mm/h = 31 mm/day
    (croptable.loc[1:, 'crop_type'] == 10),  # Ditch

]
# Ksat mm/h in 1st rainfall event Leaching experiment: 134 mm/h -> 3216 mm/day
# Ksat mm/h in 2nd rainfall event Leaching experiment range:
# 0.13 mm/h -> 3.12 mm/day
# 1.8 mm/h (1.8*24h =43.2)
# 5.3 mm/h -> 127 mm/day   <- Group B!!!
# 26.8 mm/h -> 643.2 mm/day  <- Already Group A, unlikely based on soil characteristics!!!
croptable.loc[1:, 'k_sat_z0z1'] = np.select(ksat_condition, [62, 31, 643], default=127)  # mm/day
croptable.loc[1:, 'k_sat_z2'] = np.select(ksat_condition, [62, 31, 127], default=127)  # mm/day
# croptable.loc[1:, 'k_sat_z0z1'] = 43.2  # mm/day
# croptable.loc[1:, 'k_sat_z2'] = 43.2  # mm/day

# Curve Number guidelines:
# https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044171.pdf
# https://en.wikipedia.org/wiki/Runoff_curve_number
# Will assume:
# HSG Group B (final infiltration rate 3.8–7.6 mm per hour)
# HSG Group C (final infiltration rate 1.3-3.8 mm per hour)
cn_conditions = [
    (croptable.loc[1:, 'crop_type'] == 1),  # Corn
    (croptable.loc[1:, 'crop_type'] == 2),  # Wheat
    (croptable.loc[1:, 'crop_type'] == 5),  # Beet
    (croptable.loc[1:, 'crop_type'] == 6),  # Greenery
    (croptable.loc[1:, 'crop_type'] == 7),  # Dirt Road
    (croptable.loc[1:, 'crop_type'] == 8),  # Grass Road
    (croptable.loc[1:, 'crop_type'] == 9),  # Paved Road
    (croptable.loc[1:, 'crop_type'] == 10),  # Ditch
    (croptable.loc[1:, 'crop_type'] == 11),  # Fallow
    (croptable.loc[1:, 'crop_type'] == 12),  # Hedge
    (croptable.loc[1:, 'crop_type'] == 13)  # Orchard
]
# Assumed Poor hydrologic conditions (HC)
# Hydraulic condition is based on combination factors that affect infiltration and runoff, including
# (a) density and canopy of vegetative areas,
# (b) amount of year-round cover,
# (c) amount of grass or close-seeded legumes,
# (d) percent of residue cover on the land surface (good gt 20%),and
# (e) degree of surface roughness.
CN2 = {"Corn": {"A": 72, "B": 81, "C": 88, "D": 91},  # poor HC
       "Wheat": {"A": 72, "B": 81, "C": 88, "D": 91},  # poor HC
       "Beet": {"A": 72, "B": 81, "C": 88, "D": 91},  # poor HC
       "Greenery": {"A": 35, "B": 56, "C": 70, "D": 77},  # Brush, fair HC, # Table 2-2c
       "Dirt Road": {"A": 72, "B": 82, "C": 87, "D": 89},  # Table 2-2a
       "Grass Road": {"A": 59, "B": 74, "C": 82, "D": 86},  # Farmsteads, # Table 2-2c
       "Paved Road": {"A": 98, "B": 98, "C": 98, "D": 98},  # Table 2-2a
       "Ditch": {"A": 98, "B": 98, "C": 98, "D": 98},  # Paved -> should add to vol.
       "Fallow": {"A": 30, "B": 58, "C": 71, "D": 78},  # Assumed Meadow, Table 2-2c
       "Hedge": {"A": 35, "B": 56, "C": 70, "D": 77},  # Brush, fair HC, # Table 2-2c
       "Orchard": {"A": 43, "B": 65, "C": 76, "D": 82},  # Woods-grass, fair HC, # Table 2-2c
       "Bare Soil": {"A": 77, "B": 86, "C": 91, "D": 94}  # Fallow on Table, 2-2b, but Bare Soil treatment
       }

group = [CN2["Corn"]["B"],
         CN2["Wheat"]["B"],
         CN2["Beet"]["B"],
         CN2["Greenery"]["B"],
         CN2["Dirt Road"]["B"],
         CN2["Grass Road"]["B"],
         CN2["Paved Road"]["B"],
         CN2["Ditch"]["B"],
         CN2["Fallow"]["B"],
         CN2["Hedge"]["B"],
         CN2["Orchard"]["B"]
         ]

croptable.loc[1:, 'CN2'] = np.select(cn_conditions, group, default=98)

if PC:
    saved = "D:/Documents/these_pablo/Models/BEACH2016/DataInput/Tables/DataSource/croptable_end.csv"
    # saved = "D:/Documents/these_pablo/Models/BEACH2016/DataInput/Tables/croptable.csv"
    croptable.to_csv(saved, sep=';', index=False)
    model = "D:/Documents/these_pablo/Models/BEACH2016/model_" + version + '/croptable.tbl'
    np.savetxt(model, croptable.values, fmt='%10.5f', delimiter="\t")
else:
    print ("True")
    saved = "croptable_end.csv"
    croptable.to_csv(saved, sep=',', index=False)
    np.savetxt('C:/Users/DayTimeChunks/Documents/Models/pesti-beach16/model_' + version + '/croptable.tbl',
               croptable.values, fmt='%10.5f', delimiter="\t")  # header="X\tY\tZ\tValue")
