import numpy as np
import pandas as pd
import pickle


def get_data(gen="Gen10"):
    filename = "results" + gen + '_copy.p'
    res = pickle.load(open(filename, "rb"))
    return pd.DataFrame.from_dict(res, orient='index')


gen1x = "Gen10"
df1 = get_data(gen=gen1x)
df1["Gen"] = gen1x
df = pd.concat([df1,
                ])
# Index names
df['Name'] = df['Gen'] + "/" + df.index
df = df.set_index('Name')

dfs = df[['z3_factor',
          'cZ0Z1',
          'cZ',
          'c_adr',
          'k_g',
          'gamma01',
          'gammaZ',
          'f_transp',
          'f_evap',
          'f_oc',
          'k_oc',
          'beta_runoff',
          'dt_50_aged',
          'dt_50_ab',
          'dt_50_ref',
          'epsilon_iso',
          'beta_moisture',
          'KGE-Q_out'
          ]]

# Selecting only Maximum Q-KGE
max_Q = dfs['KGE-Q_out'].max()
out = dfs.loc[(dfs['KGE-Q_out'] >= max_Q)]  # Soils
# out.iloc[:, 2:19]
select = out[['z3_factor',
              'cZ0Z1',
              'cZ',
              'c_adr',
              'k_g',
              'gamma01',
              'gammaZ',
              'f_transp',
              'f_evap',
              'f_oc',
              'k_oc',
              'beta_runoff',
              'dt_50_aged',
              'dt_50_ab',
              'dt_50_ref',
              'epsilon_iso',
              'beta_moisture',
              ]]

np.savetxt("best_vector.txt", select.values[0])
