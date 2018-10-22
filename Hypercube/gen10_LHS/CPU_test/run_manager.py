import math
import os
import shutil
import re
from copy import deepcopy
# Read in the file
lisem_runs = [183, 188, 199, 214, 225]
end_days_remain = deepcopy(lisem_runs)
event_type = [2, 1, 1, 2, 1]

# TODO: split three 2-day to 1-day events!
# Event types:
# 1 = 1 event in one day (no day break)
# 2 = event split into 2 days
# 3 = 2 or more events in 1 day

for i in range(len(lisem_runs)):
    directory = os.getcwd() + "\\res\\" + str(i)

    # Delete all previous results!
    if os.path.exists(directory):
        shutil.rmtree(directory)
        shutil.copytree(os.getcwd() + "\\res\\dummy", directory)

    # For first-time-ever events, create folder repo.
    if not os.path.exists(directory):
        shutil.copytree(os.getcwd() + "\\res\\dummy", directory)

    # Check tss/output folders exist. Need 1-extra folder.
    beach_directory = os.getcwd() + "\\tss\\output\\" + str(i)
    if not os.path.exists(directory):
        os.mkdir(beach_directory)
    # Extra BEACH tss/output folder
    if i == len(lisem_runs) - 1:
        beach_directory = os.getcwd() + "\\tss\\output\\" + str(i+1)
        if not os.path.exists(beach_directory):
            print("Created: " + beach_directory)
            os.mkdir(beach_directory)


# Preparing inputs for LISEM "run" file
target_time = [162., 290., 350., 360., 240.]
end_time = [202., 290., 350., 624., 240.]
# ts_old = [126, 126, 200, 200, 200]
# ts_new = [126, 126, 200, 200, 200]
ts_new = [9., 9., 9., 9., 9.]
ts_old = [10, 10, 10, 10, 10]
interval = 99.
map_list = []
hydro_row = []

# Following loop inputs new values to "run" file AND
# tells BEACH which (i) infiltration maps and (ii) runoff value
# to get from hydro.txt row
for t in lisem_runs:
    i = lisem_runs.index(t)
    # Infiltration Maps
    elapsed_min = ts_new[i]/60.*interval
    map_num = int(round(target_time[i]/elapsed_min))
    if target_time[i] == end_time[i]:
        map_list.append("\\infiltration.map")
    else:
        if map_num < 10:
            map_list.append("\\infilLM0.00" + str(map_num))
        else:
            map_list.append("\\infilLM0.0" + str(map_num))

    # Runoff rows
    row_ts = int(target_time[i]/(ts_new[i]/60.))
    if target_time[i] == end_time[i]:
        hydro_row.append(-1)
    else:
        hydro_row.append(row_ts)

    # Replace time step in run file
    with open('Alteck' + str(t) + '.run', 'r') as file:
        filedata = file.read()

    # Replace the target string
    target = r"Timestep=" + str(ts_old[i])
    target_new = r"Timestep=" + str(ts_new[i])

    filedata = filedata.replace(target, target_new)

    # Write the file out again
    with open('Alteck' + str(t) + '.run', 'w') as file:
        file.write(filedata)


