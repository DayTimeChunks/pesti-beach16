import numpy as np
import pandas as pd
from scipy.optimize import minimize
from datetime import datetime

import sys
sys.path.append('D:/Documents/these_pablo/Models/BEACH2016/Analysis')

from model_v2 import *


def getSoilDataCal(path, name):
    path += name + ".tss"
    obs = pd.read_table(path)
    return obs

# Function to add column ID, with first letter capitalized.
def newlabel(row, plot):
    return plot.capitalize() + '-' + str(int(row['Jdays']))


def get_dataframe(name_obs_detail, name_obs_comp, name_sim):
    """
    name_obs_detail = 'delta_det_cal' | 'conc_det_cal'
    name_obs_comp = 'delta_comp_cal' | 'conc_comp_cal'
    name_sim = 'd13C_real' | 'CONC_real'
    """
    # path = 'D:/Documents/these_pablo/Models/BEACH2016/Analysis/Data/BEACH_R/'
    path = 'C:/Users/pablo/Documents/pablo-models/pesti-beach16/Analysis/Data/BEACH_R/'

    detail_plots = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8',
                    'v4', 'v5', 'v7', 'v8', 'v9', 'v10',
                    's11', 's12', 's13']

    matches = []
    for sample in range(len(detail_plots)):
        plot = detail_plots[sample]

        # Simulated
        filename = "resM_" + plot + name_sim + ".tss"
        det_sim = pd.read_table(filename,
                                skiprows=4, delim_whitespace=True,
                                names=['Jdays', 'SIM'],
                                header=None)
        # Observed
        det_sim['IDcal'] = det_sim.apply(lambda row: newlabel(row, plot), axis=1)  # Add ID
        n = getSoilDataCal(path, name_obs_detail)  # Name observed detailed

        match = pd.merge(n, det_sim, how='inner', on='IDcal')
        matches.append(match)

    conc_det_dat = pd.concat(matches)

    # Prepare composites and concatenate
    sample_comp = ['nor', 'val', 'sou']
    matches = []

    for sample in range(len(sample_comp)):
        plot = sample_comp[sample]
        filename = "resM_" + plot + name_sim + ".tss"
        det_sim = pd.read_table(filename,
                                skiprows=4, delim_whitespace=True,
                                names=['Jdays', 'SIM'],
                                header=None)
        label = plot[0]
        det_sim['IDcal'] = det_sim.apply(lambda row: newlabel(row, label), axis=1)
        n = getSoilDataCal(path, name_obs_comp)  # Name observed composites
        match = pd.merge(n, det_sim, how='inner', on='IDcal')
        matches.append(match)
    # Concatenate composites
    conc_comp_dat = pd.concat(matches)
    conc_dat = pd.concat([conc_det_dat, conc_comp_dat])

    # Concatenate detailed AND COMPOSITES
    return conc_dat


def get_error(var):
    # Considers only soil data (detailed and composites)
    # Define observation dataframes to use
    if var == "ug.g":
        name_obs_detail = 'conc_det_cal'
        name_obs_comp = 'conc_comp_cal'
        name_sim = 'CONC_real'
    else:
        name_obs_detail = 'delta_det_cal'
        name_obs_comp = 'delta_comp_cal'
        name_sim = 'd13C_real'

    df1 = get_dataframe(name_obs_detail, name_obs_comp, name_sim)

    # Nash
    mean = df1[var].mean()
    # Diff sim vs. obs
    dfn = df1.assign(diff_sim=lambda row: (row['SIM'] - row[var]) ** 2)
    err_sim = dfn['diff_sim'].sum()
    # Variance
    dfn = dfn.assign(diff_obs=lambda row: (row[var] - mean) ** 2)
    err_obs = dfn['diff_obs'].sum()
    error = err_sim / err_obs # nash = 1 - error
    if var == 'ug.g':  # Log only for concentrations
        lnmean = np.log(df1[var]).mean()
        # Log Diff sim vs. obs
        dfn = dfn.assign(lndiff_sim=lambda row: (np.log(row['SIM']) - np.log(row[var])) ** 2)
        err_lnsim = dfn['lndiff_sim'].sum()
        # Log variance
        dfn = dfn.assign(lndiff_obs=lambda row: (np.log(row[var]) - lnmean) ** 2)
        err_lnobs = dfn['lndiff_obs'].sum()
        error += err_lnsim / err_lnobs
        error *= 0.5
    return error


def get_initial_array():
    # This section includes all non-stochastic parameters.
    # Get initial parameters, make a dictionary of the raw file.
    # import csv
    # ini_path = 'csv/initial.csv'
    # self.ini_param = {}  # Dictionary to store the values
    # with open(ini_path, 'r') as f:
    #     reader = csv.reader(f, delimiter=',')
    #     for row in reader:
    #         self.ini_param[row[0].strip()] = float(row[1])
    #
    # self.q_m3day_mean = scalar(self.ini_param.get("ave_outlet_q_m3day"))
    # self.conc_outlet_mean = scalar(self.ini_param.get("ave_outlet_conc_ugL"))
    # self.ln_conc_outlet_mean = scalar(self.ini_param.get("ave_outlet_lnconc_ugL"))
    # self.delta_outlet_mean = scalar(self.ini_param.get("ave_outlet_delta"))
    # self.conc_compNorth_mean = scalar(self.ini_param.get("ave_north_compConc_ugg"))
    # self.conc_compValley_mean = scalar(self.ini_param.get("ave_valley_compConc_ugg"))
    # self.conc_compSouth_mean = scalar(self.ini_param.get("ave_south_compConc_ugg"))
    # self.delta_compNorth_mean = scalar(self.ini_param.get("ave_north_compIso_delta"))
    # self.delta_compValley_mean = scalar(self.ini_param.get("ave_valley_compIso_delta"))
    # self.delta_compSouth_mean = scalar(self.ini_param.get("ave_south_compIso_delta"))

    # Constraint
    # dt50 (min, max), epsilon (min, max)
    # Initial, is reference for
    # dt50, epsilon
    x0 = [21, -1.37]
    return x0


b_det50 = (10, 37)
b_epsilon = (-2.0, -1.0)
bnds = (b_det50, b_epsilon)


def objective(x):
    """
    1) run the model with x
    2) get an error result, return
    :param x:
    :return:
    """
    nTimeSteps = 200  # 360
    myAlteck16 = BeachModel("mapInput/clone_nom.map", x)
    dynamicModel = DynamicFramework(myAlteck16,
                                    lastTimeStep=nTimeSteps,
                                    firstTimestep=firstTimeStep)
    dynamicModel.run()

    error_delta = get_error('d13C')
    error_conc = get_error('ug.g')
    error = 0.5*(error_delta + error_conc)
    return error


x0 = [5, -1]

# solution = minimize(objective, x0, method='SLSQP', bounds=bnds, jac=False)
# solution = minimize(objective, x0, method='L-BFGS-B', bounds=bnds, jac=False)
# solution = minimize(objective, x0) # , method='TNC') # , bounds=bnds, jac=False) # The default method is BFGS.
solution = minimize(objective, x0, method='Nelder-Mead') # , method='TNC') # , bounds=bnds, jac=False) # The default method is BFGS.
r = solution.x

# show final objective
print('Final Objective: ' + str(objective(r)))

# print solution
print('Solution')
print('x1 = ' + str(r[0]))
print('x2 = ' + str(r[1]))
