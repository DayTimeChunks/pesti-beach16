import shutil, os
import numpy as np
# Copy result TSS files
path = os.getcwd()
origin_free = path + "\\tss\\output\\0"
origin_coup = path + "\\tss\\output"
origin_res = path + "\\res"

dest_free = path + "\\visuals\\output_free\\0"
dest_coup = path + "\\visuals\\output_coup"
dest_res = path + "\\visuals\\res"


def copy_and_overwrite(from_path, to_path):
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree(from_path, to_path)


# Replace
""" Coupled """
copy_and_overwrite(origin_coup, dest_coup)
copy_and_overwrite(origin_res, dest_res)

""" Free """
# Replace non-coupled res.
# copy_and_overwrite(origin_free, dest_free)






# Test to change new values
# best_values = np.loadtxt("best_vector.txt")
#
# names = ['z3_factor',
#          'cZ0Z1',
#          'cZ',
#          'c_adr',
#          'k_g',
#          'gamma01', 'gammaZ',
#          'f_transp',
#          'f_evap',
#          'f_oc', 'k_oc',
#          'beta_runoff',
#          'dt_50_aged',
#          'dt_50_ab',
#          'dt_50_ref',
#          'epsilon_iso',
#          'beta_moisture'
#          ]
#
# print(best_values)
#
# new_values = []
# for i in range(len(names)):
#     if i != names.index("c_adr"):
#         new_values.append(best_values[i])
#     else:
#         new_values.append(best_values[i]*.5)
#
# print(new_values)
