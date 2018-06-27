import math as m
import datetime
import jdcal

# Make sure to update obs_data for the data you want...
from obs_data import *
from plotFunctions import *

"""
Import simulated data 
Observations are already imported above.
"""
# Edit observations (select relevant columns)
obs['DayMoYr'] = pd.to_datetime(obs['DayMoYr'])  # Convert to Date object
obs = obs[['DayMoYr', 'Jdays', 'VolTot.L']]  # keep desired columns

# Simulated Start time
yy = 2015
mm = 10
dd = 1
dt = datetime.date(yy, mm, dd)


# Convert model time step to Date
def convert_to_date(df):
    return dt + datetime.timedelta(days=df['dt'] - 1)


version = "v5"
path = version + "/"

if version == "v1":
    simQ = pd.read_table(path + "res_accuVol_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Volm3'],
                         header=None)
    simPerc = pd.read_table(path + "res_accuPercol_z2_m3.tss", skiprows=4, delim_whitespace=True,
                            names=['dt', 'Percolm3'], header=None)
    simLF = pd.read_table(path + "res_outLatflow_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'LFm3'],
                          header=None)

    sim = simQ.merge(simPerc, left_on='dt', right_on='dt')
    sim = sim.merge(simLF, left_on='dt', right_on='dt')
    sim['Q_sim.L'] = sim['Volm3'] * 10 ** 3
    sim['Perc_sim.L'] = sim['Percolm3'] * 10 ** 3
    sim['LF_sim.L'] = sim['LFm3'] * 10 ** 3

    sim = sim.drop(['Volm3', 'Percolm3', 'LFm3'], axis=1)

    sim['Date'] = sim.apply(lambda df: convert_to_date(df), axis=1)  # axis 2 = rows
    sim['Date'] = pd.to_datetime(sim['Date'])

    # print("sim", sim.dtypes)
    # print("obs", obs.dtypes)
    # Merge with 'outer' to keep all rows without a match
    mData = sim.merge(obs, left_on='Date', right_on='DayMoYr', how='outer')
    hydro_plot(mData)

elif version == "v5":
    simQ = pd.read_table(path + "2/res_accuVol_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Q_sim.m3'],
                         header=None)
    # print(obs.dtypes)
    # obs['Q_obs.m3'] = obs['VolTot.L'] * 2 / 10 ** 3
    # mData = obs.merge(simQ, left_on='Jdays', right_on='dt', how='outer')
    # hydro_dt(mData, m3=True)
    start = 1
    n_tests = 11
    for i in range(start, n_tests+1):
        folder = str(i)
        col = 'Nash.c' + str(i)
        nash = pd.read_table(path + folder + "/res_nash_q_m3.tss", skiprows=4, delim_whitespace=True,
                             names=['Jdays', col],
                             header=None)
        if i == 1:
            nash_all = nash
        else:
            nash_all = nash_all.merge(nash, left_on='Jdays', right_on='Jdays', how='outer')

    nash(nash_all, n_tests)
    # time = nash_all.ix[:, 0]


test = "test"
