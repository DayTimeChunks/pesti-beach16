# -*- coding: utf-8 -*-
from pcraster.framework import *
from copy import deepcopy


def writeCumHydro(model, q_obs, out_runoff_m3, out_drain_m3, tot_rain_m3,
                  out_percol_m3, out_etp_m3, outlet_latflow_m3):

    model.rain_cum_m3 += tot_rain_m3
    model.rain_obs_cum_tss.sample(model.rain_cum_m3)

    # Mass balance components
    model.tot_runoff_m3 += ifthenelse(q_obs >= 0, out_runoff_m3, 0)
    model.tot_runoff_m3_tss.sample(model.tot_runoff_m3)
    model.tot_perc_z3_m3 += ifthenelse(q_obs >= 0, out_percol_m3, 0)  # <- == 0
    model.tot_perc_z3_m3_tss.sample(model.tot_perc_z3_m3)
    model.tot_etp_m3 += ifthenelse(q_obs >= 0, out_etp_m3, 0)
    model.tot_etp_m3_tss.sample(model.tot_etp_m3)
    # model.tot_baseflow_m3 += ifthenelse(q_obs >= 0, out_baseflow_m3, 0)
    # model.tot_baseflow_m3_tss.sample(model.tot_baseflow_m3)

    model.tot_drain_m3 += ifthenelse(q_obs >= 0, out_drain_m3, 0)  # o_drain_z1_m3
    model.tot_accu_drain_m3_tss.sample(model.tot_drain_m3)

    # model.tot_ilf_m3 += ifthenelse(q_obs >= 0, outlet_lat_inflow_m3, 0)
    # model.tot_accu_i_latflow_m3_tss.sample(model.tot_ilf_m3)
    model.cum_olf_m3 += ifthenelse(q_obs >= 0, outlet_latflow_m3, 0)
    model.resW_o_cumLatflow_m3_tss.sample(model.cum_olf_m3)
    # model.tot_nlf_m3 += ifthenelse(q_obs >= 0, n_latflow_m3, 0)
    # model.tot_accu_n_latflow_m3_tss.sample(model.tot_nlf_m3)
    # model.tot_of_m3 += ifthenelse(q_obs >= 0, of_latflow_m3, 0)
    # model.tot_accu_of_latflow_m3_tss.sample(model.tot_of_m3)


def computeGlobalWaterBalance(model, ch_storage_m3, percolation, etp_m3,
                              catch_n_latflow_m3,
                              cell_lat_outflow_m3,
                              outlet_latflow_m3,
                              tot_rain_m3, out_runoff_m3, q_obs, out_drain_m3):
    # Percolation
    percol_z3_m3 = percolation[-1] * cellarea() / 1000  # m3
    out_percol_m3 = accuflux(model.ldd_subs, percol_z3_m3)
    out_percol_m3 = areatotal(out_percol_m3, model.outlet_multi)

    # Evapotranspiration
    out_etp_m3 = accuflux(model.ldd_subs, etp_m3)
    out_etp_m3 = areatotal(out_etp_m3, model.outlet_multi)

    # Lateral flow
    # model.resW_o_cellLatflow_m3_tss.sample(outlet_latflow_m3)  # TODO: Odd place for this!
    # catch_lat_outflow_m3 = areatotal(cell_lat_outflow_m3, model.is_catchment)  # Sum catchment

    # Baseflow discharge (basement)
    # accu_baseflow = accuflux(model.ldd_subs, model.baseflow * cellarea() / 1000)  # m3
    # out_baseflow_m3 = areatotal(accu_baseflow, model.outlet_multi)
    # model.out_baseflow_m3_tss.sample(out_baseflow_m3)

    accu_ch_storage_m3 = accuflux(model.ldd_subs, ch_storage_m3)
    accu_ch_storage_m3 = areatotal(accu_ch_storage_m3, model.outlet_multi)

    # Cumulative
    writeCumHydro(model, q_obs, out_runoff_m3, out_drain_m3, tot_rain_m3,
                  out_percol_m3, out_etp_m3, outlet_latflow_m3)

    reportGlobalWaterBalance(model, tot_rain_m3, out_runoff_m3, out_drain_m3,
                             outlet_latflow_m3, out_percol_m3,
                             out_etp_m3, accu_ch_storage_m3)


def reportGlobalWaterBalance(model, tot_rain_m3, out_runoff_m3, out_drain_m3,
                             outlet_latflow_m3, out_percol_m3,
                             out_etp_m3, accu_ch_storage_m3):
    model.resW_accRain_m3_tss.sample(tot_rain_m3)
    model.resW_accEtp_m3_tss.sample(out_etp_m3)
    model.resW_accRunoff_m3_tss.sample(out_runoff_m3)  # save to outlet
    model.resW_accPercol_z3_m3_tss.sample(out_percol_m3)  # Basement
    model.resW_o_accDrain_m3_tss.sample(out_drain_m3)  # Outlet discharge - Drain
    model.resW_o_cellLatflow_m3_tss.sample(outlet_latflow_m3)  # Outlet discharge - LF
    model.resW_accChStorage_m3_tss.sample(accu_ch_storage_m3)
    # model.out_accu_o_latflow_m3_tss.sample(catch_lat_outflow_m3)
    # model.resW_accChStorage_m3_tss.sample(out_ch_storage_m3)

    # GLOBAL Water
    global_mb_water = (tot_rain_m3 -
                       out_etp_m3 - out_runoff_m3 -
                       out_percol_m3 -
                       outlet_latflow_m3 -
                       out_drain_m3 -
                       accu_ch_storage_m3
                       )

    model.global_mb_water_tss.sample(global_mb_water)

    # Discharge due to lateral flow (3 types: overflow, cell-outlet flow, catchment net)
    # Overflow
    # overflow_cellvol_z0 = overflow_height_z0 * cellarea() / 1000  # m3
    # overflow_cellvol_z1 = overflow_height_z1 * cellarea() / 1000  # m3
    # overflow_cellvol_z2 = overflow_height_z2 * cellarea() / 1000  # m3
    # column_overflow = overflow_cellvol_z2 + overflow_cellvol_z1 + overflow_cellvol_z0
    # accu_of_latflow_m3 = accuflux(model.ldd_subs, column_overflow)
    # of_latflow_m3 = areatotal(accu_of_latflow_m3, model.outlet_multi)
    # model.sat_accu_overflow_m3_tss.sample(of_latflow_m3)

    # In-flow (for each cell at the outlet)
    # in_latflow_z0_m3 = lat_inflow_z0 * cellarea() / 1000  # m3
    # in_latflow_z1_m3 = lat_inflow_z1 * cellarea() / 1000
    # in_latflow_z2_m3 = lat_inflow_z2 * cellarea() / 1000
    # in_latflow_z3_m3 = lat_inflow_z3 * cellarea() / 1000
    # cell_lat_inflow_m3 = in_latflow_z0_m3 + in_latflow_z1_m3 + in_latflow_z2_m3 + in_latflow_z3_m3
    # model.report(cell_lat_inflow_m3, 'cellIn')
    # outlet_cell_inflow_m3 = areatotal(cell_lat_inflow_m3, model.outlet_multi)  # Only outlet cells
    # model.out_cell_i_latflow_m3_tss.sample(outlet_cell_inflow_m3)
    # catch_lat_inflow_m3 = areatotal(cell_lat_inflow_m3, model.is_catchment)  # Sum catchment (needed for MB)
    # model.out_accu_i_latflow_m3_tss.sample(catch_lat_inflow_m3)


def getCatchmentStorage(model):
    cell_vol_tot_m3 = deepcopy(model.zero_map)
    for layer in range(model.num_layers):
        cell_vol_tot_m3 += model.theta[layer] * model.layer_depth[layer] * cellarea() / 1000

    vol_tot_m3 = accuflux(model.ldd_subs, cell_vol_tot_m3)
    multi_vol_tot_m3 = areatotal(vol_tot_m3, model.outlet_multi)
    model.storage_m3_tss.sample(multi_vol_tot_m3)


def getAverageMoisture(model):
    prop_sat_layers = []
    for layer in range(model.num_layers):
        prop_sat = model.theta[layer]/model.theta_sat[layer]
        prop_sat_ave = areaaverage(prop_sat, model.is_catchment)
        prop_sat_layers.append(prop_sat_ave)

    model.resW_z0_thetaPropSat.sample(prop_sat_layers[0])
    model.resW_z1_thetaPropSat.sample(prop_sat_layers[1])
    model.resW_z2_thetaPropSat.sample(prop_sat_layers[2])
    model.resW_z3_thetaPropSat.sample(prop_sat_layers[3])


def getRestitution():
    # rest_obs = tot_rain_m3 / q_obs
    # model.rest_obs_tss.sample(rest_obs)
    pass
