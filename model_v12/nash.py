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
    

def repNashOutConc(model, q_obs, outlet_light_export, 
                   catch_latflux_light, catch_drain_light, catch_runoff_light, 
                   catch_volat_light, catch_deg_light):
    # TODO:
    # Include nash computation

    # Outlet-specific, Cumulative masses
    model.cum_exp_L_g += ifthenelse(q_obs >= 0, outlet_light_export, scalar(0))
    model.resM_cumEXP_Smet_g_tss.sample(model.cum_exp_L_g)

    model.cum_latflux_L_g += ifthenelse(q_obs >= 0, catch_latflux_light, scalar(0))
    model.cum_latflux_L_g_tss.sample(model.cum_latflux_L_g)

    model.cum_adr_L_g += ifthenelse(q_obs >= 0, catch_drain_light, scalar(0)) 
    model.cum_adr_L_g_tss.sample(model.cum_adr_L_g)

    model.cum_roZ0_L_g += ifthenelse(q_obs >= 0, catch_runoff_light, scalar(0)) 
    model.cum_roZ0_L_g_tss.sample(model.cum_roZ0_L_g)
    
    # Other sinks
    model.cum_volatZ0_L_g += ifthenelse(q_obs >= 0, catch_volat_light, scalar(0)) 
    model.cum_volatZ0_L_g_tss.sample(model.cum_volatZ0_L_g)

    model.cum_deg_L_g += ifthenelse(q_obs >= 0, catch_deg_light, scalar(0)) 
    model.cum_deg_L_g_tss.sample(model.cum_deg_L_g)
    

def repNashOutIso(model, 
                  out_delta, 
                  roff_delta, latflux_delta, drain_delta):
    # TODO:
    # Include nash computation

    model.resM_outISO_d13C_tss.sample(out_delta)
    model.resM_outISO_ROFF_d13C_tss.sample(roff_delta)
    model.resM_outISO_LF_d13C_tss.sample(latflux_delta)
    model.resM_outISO_ADR_d13C_tss.sample(drain_delta)


def repNashOutCombined():
    pass


# Soils
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
