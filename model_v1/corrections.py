
from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

"""
Converting clone.map to nominal type 
"""
# To assign ordinal, use "order(expression)"
# oldclone = readmap("clone")
# dem = readmap("dem_slope")
# res = boolean(dem)
# res = nominal(res)
# report(res, 'clone_nom.map')  # stores a ".map" file


"""
Correcting the outlet location
"""
# res = readmap("zzTest") # the result of the accuflux function (via ldd)
# outlet = ifthenelse(res == mapmaximum(res), nominal(1), nominal(0))
# report(outlet, 'outlet_true.map')

"""
fields_cover.map

Replacing missing values in landuse (due to ArcGis Errors"
with landuse = 5 (check again later)
"""
# fields = readmap("landuse")
# fields_true = cover(fields, 5)
# dem = readmap("dem_slope")
# res = boolean(dem)
# fields_fin = ifthen(res, fields_true)  # Reduce the map to the catchment extent
# report(fields_fin, 'fields_cover.map')  # store a ".map" file

# Checking how to get the value only at the outlet
outlet = readmap("outlet_true")
perc = readmap("wetness")
out = boolean(outlet)

q = ifthenelse(out, perc, scalar(0))
new = areaaverage(q, outlet)
cell = cellvalue(outlet, 30, 30)
print(cell)
aguila(new, outlet)
