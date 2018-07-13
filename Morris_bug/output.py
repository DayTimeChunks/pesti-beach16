# -*- coding: utf-8 -*-
from pcraster.framework import *


# HYDRO
def defineHydroTSS(model):
    # Rain
    model.resW_accRain_m3_tss = TimeoutputTimeseries("resW_accRain_m3", model, nominal("outlet_v3"), noHeader=False)
    # Runoff
    model.resW_accRunoff_m3_tss = TimeoutputTimeseries("resW_accRunoff_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.tot_runoff_m3_tss = TimeoutputTimeseries("resW_totRunoff_m3", model, nominal("outlet_v3"), noHeader=False)
    # Percolation
    model.resW_accPercol_Bsmt_m3_tss = TimeoutputTimeseries("resW_accPercol_Bsmt_m3", model, nominal("outlet_v3"),
                                                           noHeader=False)
    # Deep percolation Basement
    model.tot_perc_z3_m3_tss = TimeoutputTimeseries("resW_totPercol_z3_m3", model, nominal("outlet_v3"),
                                                   noHeader=False)
    # ETP
    model.resW_accEtp_m3_tss = TimeoutputTimeseries("resW_accEtp_m3", model, nominal("outlet_v3"), noHeader=False)
    model.resW_accEvap_m3_tss = TimeoutputTimeseries("resW_accEvap_m3", model, nominal("outlet_v3"), noHeader=False)
    model.resW_accTransp_m3_tss = TimeoutputTimeseries("resW_accTransp_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.tot_etp_m3_tss = TimeoutputTimeseries("resW_totEtp_m3", model, nominal("outlet_v3"), noHeader=False)
    
    # Baseflow
    model.out_baseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", model, nominal("outlet_v3"),
                                                    noHeader=False)
    # model.tot_baseflow_m3_tss = TimeoutputTimeseries("resW_totBaseflow_m3", model, nominal("outlet_v3"),
    #                                                 noHeader=False)
    # LF Drainage
    model.resW_o_accDrain_m3_tss = TimeoutputTimeseries("resW_o_accDrain_m3", model, nominal("outlet_v3"),
                                                       noHeader=False)
    model.tot_accu_drain_m3_tss = TimeoutputTimeseries("resW_o_totDrain_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)  # Cumulative ADR
    # LF options
    model.sat_accu_overflow_m3_tss = TimeoutputTimeseries("resW_of_accLatflow_m3", model, nominal("outlet_v3"),
                                                         noHeader=False)
    
    model.tot_accu_of_latflow_m3_tss = TimeoutputTimeseries("resW_of_totLatflow_m3", model, nominal("outlet_v3"),
                                                           noHeader=False)
    # Inflow
    model.out_cell_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_cellLatflow_m3", model, nominal("outlet_v3"),
                                                          noHeader=False)
    model.out_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_accLatflow_m3", model, nominal("outlet_v3"),
                                                          noHeader=False)
    model.tot_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_totLatflow_m3", model, nominal("outlet_v3"),
                                                          noHeader=False)
    # Outflow
    model.resW_o_cellLatflow_m3_tss = TimeoutputTimeseries("resW_o_cellLatflow_m3", model, nominal("outlet_v3"),
                                                          noHeader=False)  # Outlet LF
    # model.out_accu_o_latflow_m3_tss = TimeoutputTimeseries("resW_o_accLatflow_m3", model, nominal("outlet_v3"),
    #                                                       noHeader=False)
    model.resW_o_cumLatflow_m3_tss = TimeoutputTimeseries("resW_o_cumLatflow_m3", model, nominal("outlet_v3"),
                                                         noHeader=False)
    # model.out_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_accLatflow_m3", model, nominal("outlet_v3"),
    #                                                       noHeader=False)
    # model.tot_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_totLatflow_m3", model, nominal("outlet_v3"),
    #                                                       noHeader=False)
    
    model.resW_accChStorage_m3_tss = TimeoutputTimeseries("resW_accChStorage_m3", model, nominal("outlet_v3"),
                                                         noHeader=False)
    model.global_mb_water_tss = TimeoutputTimeseries("resW_global_waterMB", model, nominal("outlet_v3"),
                                                    noHeader=False)
    model.storage_m3_tss = TimeoutputTimeseries("resW_accStorage_m3", model, nominal("outlet_v3"), noHeader=False)
    
    # Basement layer analysis
    # model.resW_accBAL_z3_m3_tss = TimeoutputTimeseries("resW_accBAL_z3_m3", model, nominal("outlet_v3"),
    #                                                   noHeader=False)
    # model.resW_accInfil_z3_m3_tss = TimeoutputTimeseries("resW_accInfil_z3_m3", model, nominal("outlet_v3"),
    #                                                     noHeader=False)
    # model.resW_accLF_z3_m3_tss = TimeoutputTimeseries("resW_accLF_z3_m3", model, nominal("outlet_v3"),
    #                                                  noHeader=False)
    # model.resW_accETP_z3_m3_tss = TimeoutputTimeseries("resW_accETP_z3_m3", model, nominal("outlet_v3"),
    #                                                   noHeader=False)
    
    model.resW_accBAL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accBAL_Bsmt_m3", model, nominal("outlet_v3"),
                                                        noHeader=False)
    model.resW_accInfil_Bsmt_m3_tss = TimeoutputTimeseries("resW_accInfil_Bsmt_m3", model, nominal("outlet_v3"),
                                                          noHeader=False)
    model.resW_accLF_Bsmt_m3_tss = TimeoutputTimeseries("resW_accLF_Bsmt_m3", model, nominal("outlet_v3"),
                                                       noHeader=False)
    model.resW_accETP_Bsmt_m3_tss = TimeoutputTimeseries("resW_accETP_Bsmt_m3", model, nominal("outlet_v3"),
                                                        noHeader=False)
    
    # Storage
    model.resW_accVOL_z0_m3_tss = TimeoutputTimeseries("resW_accVOL_z0_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.resW_accVOL_z1_m3_tss = TimeoutputTimeseries("resW_accVOL_z1_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.resW_accVOL_z2_m3_tss = TimeoutputTimeseries("resW_accVOL_z2_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.resW_accVOL_z3_m3_tss = TimeoutputTimeseries("resW_accVOL_z3_m3", model, nominal("outlet_v3"),
                                                      noHeader=False)
    model.resW_accVOL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accVOL_Bsmt_m3", model, nominal("outlet_v3"),
                                                        noHeader=False)

    # Analysis
    # This is 'q' as time series.
    model.i_Q_m3_tss = TimeoutputTimeseries("resW_i_accVol_m3", model, nominal("outlet_v3"), noHeader=False)
    model.o_Q_m3_tss = TimeoutputTimeseries("resW_o_accVol_m3", model, nominal("outlet_v3"), noHeader=False)
    model.tot_Q_m3_tss = TimeoutputTimeseries("resW_tot_accVol_m3", model, nominal("outlet_v3"), noHeader=False)
    model.q_obs_cum_tss = TimeoutputTimeseries("resW_cum_q_obs_m3", model, nominal("outlet_v3"),
                                              noHeader=False)  # Equivalent to net_Q
    model.rain_obs_cum_tss = TimeoutputTimeseries("resW_cum_rain_obs_m3", model, nominal("outlet_v3"),
                                                 noHeader=False)  # Equivalent to net_Q
    model.rest_obs_tss = TimeoutputTimeseries("resW_q_restit_obs_m3", model, nominal("outlet_v3"),
                                             noHeader=False)  # = rain/q_obs
    model.q_sim_cum_tss = TimeoutputTimeseries("resW_cum_q_sim_m3", model, nominal("outlet_v3"),
                                              noHeader=False)  # Sum sim discharge (if obs available).
    model.q_sim_ave_tss = TimeoutputTimeseries("resW_q_sim_ave_m3", model, nominal("outlet_v3"),
                                              noHeader=False)  # This is 'Nash_q' as time series.
    
    
def definePestTSS(model):
    # PESTI
    # Pesticide
    model.global_mb_pest_tss = TimeoutputTimeseries("resM_global_mb_pest", model, nominal("outlet_v3"),
                                                   noHeader=False)
    model.resM_accAPP_g_tss = TimeoutputTimeseries("resM_accAPP_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accVOLAT_L_tss = TimeoutputTimeseries("resM_accVOLAT_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accRO_L_tss = TimeoutputTimeseries("resM_accRO_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accDEG_L_tss = TimeoutputTimeseries("resM_accDEG_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accAGED_L_tss = TimeoutputTimeseries("resM_accAGED_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accAGED_DEG_L_tss = TimeoutputTimeseries("resM_accAGED_DEG_L", model, nominal("outlet_v3"),
                                                       noHeader=False)
    model.resM_accDEGz0_L_tss = TimeoutputTimeseries("resM_accDEGz0_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accLCHz0_L_tss = TimeoutputTimeseries("resM_accLCHz0_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accDPz1_L_tss = TimeoutputTimeseries("resM_accDPz1_L", model, nominal("outlet_v3"), noHeader=False)
    
    model.resM_accDP_L_tss = TimeoutputTimeseries("resM_accDP_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accADR_L_tss = TimeoutputTimeseries("resM_accADR_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accLF_L_tss = TimeoutputTimeseries("resM_accLF_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accBF_L_tss = TimeoutputTimeseries("resM_accBF_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accCHS_L_tss = TimeoutputTimeseries("resM_accCHS_L", model, nominal("outlet_v3"), noHeader=False)
    model.resM_accCHS_AGED_L_tss = TimeoutputTimeseries("resM_accCHS_AGED_L", model, nominal("outlet_v3"),
                                                       noHeader=False)
    
    model.resM_EXP_Smet_g_tss = TimeoutputTimeseries("resM_EXP_Smet_g", model, nominal("outlet_v3"),
                                                    noHeader=False)  # Total outlet mass (g) exports (light fraction only)
    # Concentrations outlet
    model.resM_oCONC_ugL_tss = TimeoutputTimeseries("resM_oCONC_ugL", model, nominal("outlet_v3"),
                                                   noHeader=False)  # Total outlet conc (ug/L)
    model.resM_oCONC_ROFF_ugL_tss = TimeoutputTimeseries("resM_oCONC_ROFF_ugL", model, nominal("outlet_v3"),
                                                        noHeader=False)  # Runoff outlet conc (ug/L)
    model.resM_oCONC_LF_ugL_tss = TimeoutputTimeseries("resM_oCONC_LF_ugL", model, nominal("outlet_v3"),
                                                      noHeader=False)  # Latflow outlet conc (ug/L)
    
    model.resM_oCONC_ADR_ugL_tss = TimeoutputTimeseries("resM_oCONC_ADR_ugL", model, nominal("outlet_v3"),
                                                       noHeader=False)  # Artificial drainage outlet conc (ug/L)
    # Isotopes outlet
    model.resM_outISO_d13C_tss = TimeoutputTimeseries("resM_outISO_d13C", model, nominal("outlet_v3"),
                                                     noHeader=False)  #
    model.resM_outISO_ROFF_d13C_tss = TimeoutputTimeseries("resM_outISO_ROFF_d13C", model, nominal("outlet_v3"),
                                                          noHeader=False)  # Runoff outlet
    model.resM_outISO_LF_d13C_tss = TimeoutputTimeseries("resM_outISO_LF_d13C", model, nominal("outlet_v3"),
                                                        noHeader=False)  # Latflow outlet
    model.resM_outISO_ADR_d13C_tss = TimeoutputTimeseries("resM_outISO_ADR_d13C", model, nominal("outlet_v3"),
                                                         noHeader=False)  # Artificial drainage outlet
    
    # Cumulative Pesticide
    model.cum_degZ0_L_g_tss = TimeoutputTimeseries("resM_cumDEGz0_L", model, nominal("outlet_v3"),
                                                  noHeader=False)  # Deg z0
    model.cum_deg_L_g_tss = TimeoutputTimeseries("resM_cumDEG_L", model, nominal("outlet_v3"),
                                                noHeader=False)  # Deg z0
    model.cum_aged_deg_L_g_tss = TimeoutputTimeseries("resM_cumAGE_DEG_L", model, nominal("outlet_v3"),
                                                     noHeader=False)
    model.resM_cumLCHz0_L_g_tss = TimeoutputTimeseries("resM_cumLCHz0_L", model, nominal("outlet_v3"),
                                                      noHeader=False)  # Leaching z0
    model.cum_roZ0_L_g_tss = TimeoutputTimeseries("resM_cumROz0_L", model, nominal("outlet_v3"),
                                                 noHeader=False)  # Runoff
    
    model.cum_volatZ0_L_g_tss = TimeoutputTimeseries("resM_cumVOLATz0_L", model, nominal("outlet_v3"),
                                                    noHeader=False)  # Runoff
    
    model.cum_adr_L_g_tss = TimeoutputTimeseries("resM_cumADR_L", model, nominal("outlet_v3"),
                                                noHeader=False)  # Art. drainage
    model.cum_latflux_L_g_tss = TimeoutputTimeseries("resM_cumLF_L", model, nominal("outlet_v3"),
                                                    noHeader=False)  # Soil column, outlet cells
    
    model.resM_cumEXP_Smet_g_tss = TimeoutputTimeseries("resM_cumEXP_Smet_g", model, nominal("outlet_v3"),
                                                       noHeader=False)  # Total cum. outlet mass (g) exports

    

def getThetaProportion(model):
    # Theta average proportion to saturation
    model.resW_z0_thetaPropSat = TimeoutputTimeseries("resW_z0_thetaPropSat", model, nominal("outlet_v3"),
                                                     noHeader=False)
    model.resW_z1_thetaPropSat = TimeoutputTimeseries("resW_z1_thetaPropSat", model, nominal("outlet_v3"),
                                                     noHeader=False)
    model.resW_z2_thetaPropSat = TimeoutputTimeseries("resW_z2_thetaPropSat", model, nominal("outlet_v3"),
                                                     noHeader=False)
    model.resW_z3_thetaPropSat = TimeoutputTimeseries("resW_z3_thetaPropSat", model, nominal("outlet_v3"),
                                                     noHeader=False)
    model.resW_Bsmt_thetaPropSat = TimeoutputTimeseries("resW_Bsmt_thetaPropSat", model, nominal("outlet_v3"),
                                                       noHeader=False)