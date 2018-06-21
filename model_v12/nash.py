# -*- coding: utf-8 -*-
from pcraster.framework import *


def reportNashHydro(model, q_obs, tot_vol_disch_m3):
    # Global ave discharge of data range = 260.07 m3/day
    model.q_obs_cum += ifthenelse(q_obs >= 0, q_obs, 0)
    model.q_sim_cum += ifthenelse(q_obs >= 0, tot_vol_disch_m3, 0)
    model.q_diff += ifthenelse(q_obs >= 0, (tot_vol_disch_m3 - q_obs) ** 2, 0)
    model.q_var += ifthenelse(q_obs >= 0, (q_obs - model.q_m3day_mean) ** 2, 0)
    nash_q = 1 - (model.q_diff / model.q_var)
    model.nash_q_tss.sample(nash_q)

    model.q_obs_cum_tss.sample(model.q_obs_cum)
    model.q_sim_cum_tss.sample(model.q_sim_cum)

    # Dynamic mean leads to incorrect variance calc. upon cumulative addition (omitted)

    # model.q_sim_ave = model.q_sim_cum / model.days_cum
    # model.q_sim_ave_tss.sample(model.q_sim_ave)


def repNashOutConc(model, conc_outlet_obs, conc_ugL):
    # Nash computation consider normal and ln-transformed concentrations,
    # with the latter accounting for variance at low concentration ranges
    model.out_conc_diff += ifthenelse(conc_outlet_obs >= 0, (conc_ugL - conc_outlet_obs) ** 2, 0)
    model.out_conc_var += ifthenelse(conc_outlet_obs >= 0, (conc_ugL - model.conc_outlet_mean) ** 2, 0)
    model.out_lnconc_diff += ifthenelse(conc_outlet_obs >= 0, (ln(conc_ugL) - ln(conc_outlet_obs)) ** 2, 0)
    model.out_lnconc_var += ifthenelse(conc_outlet_obs >= 0, (ln(conc_ugL) - model.ln_conc_outlet_mean) ** 2, 0)
    normal_term = model.out_conc_diff / model.out_conc_var
    ln_term = model.out_lnconc_diff / model.out_lnconc_var
    nash_outlet_conc = 1 - 0.5 * (normal_term + ln_term)
    model.nash_outlet_conc_tss.sample(nash_outlet_conc)


def repNashOutIso(model, iso_outlet_obs, out_delta,
                  roff_delta, latflux_delta, drain_delta):

    model.out_iso_diff += ifthenelse(iso_outlet_obs < 1e6, (out_delta - iso_outlet_obs) ** 2, 0)
    model.out_iso_var += ifthenelse(iso_outlet_obs < 1e6, (out_delta - model.delta_outlet_mean) ** 2, 0)
    nash_outlet_iso = 1 - (model.out_iso_diff / model.out_iso_var)
    model.nash_outlet_iso_tss.sample(nash_outlet_iso)

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


def repCumOutMass(model, conc_outlet_obs, outlet_light_export,
                  catch_latflux_light, catch_drain_light, catch_runoff_light,
                  catch_volat_light, catch_deg_light):
    # Outlet-specific, Cumulative masses
    model.cum_exp_L_g += ifthenelse(conc_outlet_obs > 0, outlet_light_export, scalar(0))
    model.resM_cumEXP_Smet_g_tss.sample(model.cum_exp_L_g)

    model.cum_latflux_L_g += ifthenelse(conc_outlet_obs > 0, catch_latflux_light, scalar(0))
    model.cum_latflux_L_g_tss.sample(model.cum_latflux_L_g)

    model.cum_adr_L_g += ifthenelse(conc_outlet_obs > 0, catch_drain_light, scalar(0))
    model.cum_adr_L_g_tss.sample(model.cum_adr_L_g)

    model.cum_roZ0_L_g += ifthenelse(conc_outlet_obs > 0, catch_runoff_light, scalar(0))
    model.cum_roZ0_L_g_tss.sample(model.cum_roZ0_L_g)

    # Not in outlet, but relevant cumulative sinks
    model.cum_volatZ0_L_g += catch_volat_light
    model.cum_volatZ0_L_g_tss.sample(model.cum_volatZ0_L_g)

    model.cum_deg_L_g += catch_deg_light
    model.cum_deg_L_g_tss.sample(model.cum_deg_L_g)

    # TODO:
    # Report outlet mass (not cumulative), and components

