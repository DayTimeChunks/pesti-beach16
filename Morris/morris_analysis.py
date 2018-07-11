# -*- coding: utf-8 -*-
from SALib.analyze import morris as moa
from morris_test import *
# import numpy as np


# TODO: Evaluates only the last observation,
# so need to compute first the cumulative for each
def getTestOutput(variable, runs):
    """
    :param variable: name of the variable to extract
    :param runs: number of runs
    :return: the output vector, length = no. of runs
    """
    Y = np.zeros([runs])
    folder = 1
    for i in range(runs):
        path = str(folder) + "/" + variable + ".tss"
        # res = pd.read_table(path, header=None)
        res = np.loadtxt(path, float, skiprows=4, usecols=[1])
        Y[i] = res[-1]
        folder += 1
    return Y


def getSiList(ouput_vars, problem, param_values, grid_jump, p, runs):
    morris_results = []
    for indx in range(len(ouput_vars)):
        var_name = ouput_vars[indx]
        Y = getTestOutput(var_name, runs)
        Si = moa.analyze(problem, param_values, Y, num_resamples=1000, print_to_console=True,
                         grid_jump=grid_jump, num_levels=p)
        morris_results.append(Si)
    return morris_results


def getInputMatrix(path):
    path += "input_vectors.txt"
    return np.loadtxt(path, float)


# To test, uncomment:
# get_runs(get_param_values())
# First test of soil-water holding capacity on:
ouput_vars = [
    # Outlet - Hydro
    "resW_tot_accVol_m3",  # Tot Q-sim
    "resW_accEtp_m3",
    "resW_accRunoff_m3",
    "resW_o_accDrain_m3",
    "resW_o_cellLatflow_m3",
    "resW_accChStorage_m3",  # needed?
    "resW_accBaseflow_m3",
    # Outlet - Mass
    "resM_EXP_Smet_g",  # TOT = RO + ADR + LF
    "resM_oCONC_ugL",  # Tot concentration outlet
    "resM_oCONC_ROFF_ugL",
    "resM_oCONC_LF_ugL",
    "resM_oCONC_ADR_ugL",
    # Outlet - ISO
    "resM_outISO_d13C",  # TOT = RO + ADR + LF
    "resM_outISO_ROFF_d13C",
    "resM_outISO_LF_d13C",
    "resM_outISO_ADR_d13C",
    # Soils - Mass
    "resM_accDEG_L",
    "resM_accDEGz0_L",
    "resM_accAGED_L",
    "resM_accAGED_DEG_L",
    "resM_accVOLAT_L",
    "resM_accRO_L",
    "resM_accLCHz0_L",  # Leaching z0
    "resM_accDPz1_L",  # Leaching z1
    "resM_accADR_L",
    "resM_accLF_L",  # Outlet cell mass-loss
    # Nash - outlet
    "resNash_q_m3",  # Q
    "resNash_outConc_ugL",  # Conc
    "resNash_outIso_delta",  # Iso
    # Nash - soils
    "nashN_compConc_L",
    "nashV_compConc_L",
    "nashS_compConc_L",
    "nashN_compIso",
    "nashV_compIso",
    "nashS_compIso"
]

def run_test():
    p = 4.0
    grid_jump = 2
    r = 4
    params_test = get_param_values(p=p, grid_jump=grid_jump, r=r)
    morris_results = getSiList(ouput_vars, problem, params_test, 2, p, get_runs(params_test))
    saveMorris(ouput_vars, morris_results)

run_test()





