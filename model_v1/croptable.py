import numpy as np
import pandas as pd

# Crop table paramters
sow_mm = {"Corn": 4, "Wheat": 3, "Beet": 3}


path = "D:/Documents/these_pablo/Models/BEACH2016/DataInput/Tables/DataSource/"
path += "croptable_start.csv"

croptable = pd.read_csv(path, sep=";")

crop_conditions = [
    (croptable['croptype'] == 1),
    (croptable['croptype'] == 2),
    (croptable['croptype'] == 5)
]
sow_mm = [4, 3, 3]

croptable['sow_yy'] = 2016
croptable['sow_mm'] = np.select(crop_conditions, sow_mm, default=0)



landuse = {"Corn": 1, "Wheat": 2, "Oats": 3, "Alfalfa": 4, "Beet": 5,
               "Greenery": 6, "Dirt Road": 7, "Grass Road": 8, "Paved Road": 9, "Ditch": 10,
               "Fallow": 11, "Hedge": 12, "Orchard": 13}