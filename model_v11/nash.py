# -*- coding: utf-8 -*-
from pcraster.framework import *


def reportNashHydro(model, q_obs, tot_vol_disch_m3):
    model.q_obs_cum += ifthenelse(q_obs >= 0, q_obs, 0)
    model.days_cum += ifthenelse(q_obs >= 0, scalar(1), scalar(0))  # Total days with data
    model.q_sim_cum += ifthenelse(q_obs >= 0, tot_vol_disch_m3, 0)
    model.q_sim_ave = model.q_sim_cum / model.days_cum
    model.q_diff += ifthenelse(q_obs >= 0, (tot_vol_disch_m3 - q_obs) ** 2, 0)
    model.q_var += ifthenelse(q_obs >= 0, (q_obs - 260.07) ** 2, 0)  # Mean discharge of data range = 260.07 m3/day
    nash_q = 1 - (model.q_diff / model.q_var)
    model.nash_q_tss.sample(nash_q)

    model.q_obs_cum_tss.sample(model.q_obs_cum)
    model.q_sim_cum_tss.sample(model.q_sim_cum)
    model.q_sim_ave_tss.sample(model.q_sim_ave)
    

def reportNashConcComposites(model, north_ave_conc, valley_ave_conc, south_ave_conc):
    """
    Nash Soil Concentrations
    1) Get the mean for each soil composite for entire year
    North = 1.909193 ug/g soil
    Talweg = 2.261839 ug/g soil
    South = 2.389668 ug/g soil

    2) The variance for each transect is
    var_north = (conc_north - mean_north)**2, if conc_north > 0

    3) Nash will be:
    1 - (conc_north_diff/var_north +  valley + south)
    """
    # Nash
    conc_north_obs = timeinputscalar('northConc.tss', ordinal("north_ave"))
    model.northConc_diff += ifthenelse(conc_north_obs > 0, (north_ave_conc - conc_north_obs) ** 2, scalar(0))
    model.northConc_var += ifthenelse(conc_north_obs > 0, (north_ave_conc - scalar(1.909193)) ** 2,
                                     scalar(0))  # ug/g
    
    # Nash
    conc_valley_obs = timeinputscalar('valleyConc.tss', ordinal("valley_ave"))
    model.valleyConc_diff += ifthenelse(conc_valley_obs > 0, (valley_ave_conc - conc_valley_obs) ** 2, scalar(0))
    model.valleyConc_var += ifthenelse(conc_valley_obs > 0, (valley_ave_conc - scalar(2.261839)) ** 2,
                                      scalar(0))  # ug/g

    conc_south_obs = timeinputscalar('southConc.tss', ordinal("south_ave"))
    model.southConc_diff += ifthenelse(conc_south_obs > 0, (south_ave_conc - conc_south_obs) ** 2, scalar(0))
    model.southConc_var += ifthenelse(conc_south_obs > 0, (south_ave_conc - scalar(2.389668)) ** 2,
                                     scalar(0))  # ug/g
    nash_compConc_L = 1 - ((model.northConc_diff / model.northConc_var) * 1 / 3 +
                           (model.valleyConc_diff / model.valleyConc_var) * 1 / 3 +
                           (model.southConc_diff / model.southConc_var) * 1 / 3)
    model.nash_compConc_L_tss.sample(nash_compConc_L)


def reportNashDeltaComposites():
    pass
