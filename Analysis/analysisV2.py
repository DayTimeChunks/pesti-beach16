from obs_data import *
from plotFunctions import *
import math as m
import datetime
import jdcal
import ggplot


"""
Import simulated data 
Observations are already imported above.
"""
# Edit observations (select relevant columns)
obs['DayMoYr'] = pd.to_datetime(obs['DayMoYr'])  # Convert to Date object
obs_w = obs[['DayMoYr', 'VolTot.L']]  # keep desired columns


path1 = "v1/"
version = "v2"
running = False
if running:
    path = "/"
else:
    path = "../model_" + version + "/"

# C:/Users/DayTimeChunks/Documents/Models/pesti-beach16
simQ = pd.read_table(path + "res_accuVol_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Volm3'], header=None)
simLF = pd.read_table(path + "res_accuLatflow_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'LFm3'], header=None)
# simPerc = pd.read_table(path1 + "res_accuPercol_z2_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Percolm3'], header=None)
simStore = pd.read_table(path + "res_accuStorage_m3.tss", skiprows=4, delim_whitespace=True, names=['dt', 'Percolm3'], header=None)


sim = simQ.merge(simLF, left_on='dt', right_on='dt')
sim = sim.merge(simStore, left_on='dt', right_on='dt')


sim['Q_sim.L'] = sim['Volm3'] * 10 ** 3
# sim['Perc_sim.L'] = sim['Percolm3'] * 10 ** 3
sim['LF_sim.L'] = sim['LFm3'] * 10 ** 3

sim = sim.drop(['Volm3', 'LFm3'], axis=1)
# sim = sim.drop(['Percol_sim.L'], axis=1)

# Convert model time step to Date
def convert_to_date(df):
    # Simulated Start time
    yy = 2015
    mm = 10
    dd = 1
    dt = datetime.date(yy, mm, dd)
    return dt + datetime.timedelta(days=df['dt']-1)


sim['Date'] = sim.apply(lambda df: convert_to_date(df), axis=1)  # axis 1 = rows
sim['Date'] = pd.to_datetime(sim['Date'])

# print("sim", sim.dtypes)
# print("obs", obs.dtypes)
# Merge with 'outer' to keep all rows without a match
mData = sim.merge(obs_w, left_on='Date', right_on='DayMoYr', how='outer')

def hydro_v2(data):
    fig, ax1 = plt.subplots()
    # time = data['Date']
    time = data['dt']
    obs_vol = data['VolTot.L']
    sim_LF = data['LF_sim.L'] #+ data['Perc_sim.L']
    sim_vol = data['Q_sim.L'] #+ data['Perc_sim.L']

    legend = ["Q-obs (L/d)", "Q (L/d)", "LF (L/d)"]
    color_sequence = ['#d62728',
                      'darkviolet'
                      '#2ca02c',
                      '#1f77b4'
                      ]

    ax1.plot(time, obs_vol, label=legend[0])
    ax1.plot(time, sim_vol, linestyle='dashdot', label=legend[1])
    ax1.plot(time, sim_LF, linestyle='dashed', label=legend[2])
    # ax1.plot(time, obs_vol, color_sequence[0], linestyle='dashdot', label=legend[0])
    title = None
    if title:
        plt.title(title)
    # handles, labels = ax1.get_legend_handles_labels()
    # lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.3, 0.9), fancybox=True, framealpha=0.7,
    #                  ncol=2)
    # if save:
    #     plt.savefig('../Figures/' + str(title) + '.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

    # add_margin(ax1, x=0.01, y=0.01)
    plt.legend(loc='upper left', fancybox=True, framealpha=0.7)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Flux (L/day)')
    plt.show()

print("mData", mData.dtypes)
hydro_v2(mData)

test = "test"

