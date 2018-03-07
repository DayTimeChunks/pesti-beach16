from obs_data import *
from plotFunctions import *
import math as m
import datetime
import jdcal

"""
Import simulated data 
Observations are already imported above.
"""
# Edit observations (select relevant columns)
obs['DayMoYr'] = pd.to_datetime(obs['DayMoYr'])  # Convert to Date object
obs = obs[['DayMoYr', 'VolTot.L']]  # keep desired columns

# Simulated Start time
yy = 2015
mm = 10
dd = 1
dt = datetime.date(yy, mm, dd)

version = "v1"
path = version + "/"

simQ = pd.read_table(path + "res_accuVol_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Volm3'], header=None)
simPerc = pd.read_table(path + "res_accuPercol_z2_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Percolm3'], header=None)
simLF = pd.read_table(path + "res_outLatflow_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'LFm3'], header=None)

sim = simQ.merge(simPerc, left_on='dt', right_on='dt')
sim = sim.merge(simLF, left_on='dt', right_on='dt')

sim['Q_sim.L'] = sim['Volm3'] * 10 ** 3
sim['Perc_sim.L'] = sim['Percolm3'] * 10 ** 3
sim['LF_sim.L'] = sim['LFm3'] * 10 ** 3

sim = sim.drop(['Volm3', 'Percolm3', 'LFm3'], axis=1)

# Convert model time step to Date
def convert_to_date(df):
    return dt + datetime.timedelta(days=df['dt']-1)


sim['Date'] = sim.apply(lambda df: convert_to_date(df), axis=1)  # axis 1 = rows
sim['Date'] = pd.to_datetime(sim['Date'])

# print("sim", sim.dtypes)
# print("obs", obs.dtypes)
# Merge with 'outer' to keep all rows without a match
mData = sim.merge(obs, left_on='Date', right_on='DayMoYr', how='outer')

hydro_plot(mData)

test = "test"

