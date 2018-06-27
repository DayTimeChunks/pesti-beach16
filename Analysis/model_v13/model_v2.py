# -*- coding: utf-8 -*-

# from time import *
import os
from copy import deepcopy
from datetime import datetime

from pcraster._pcraster import *
from pcraster.framework import *

from applications import getApplications
from hydro_v2 import *
from nash import *
from pesti_v2 import *
from soil_samples import *
from test_suite import *
from balance import *

print(os.getcwd())

global state
state = -1

morris = False
if morris:
    from morris_test import *
else:
    runs = 2

"""
model_v13
Model tested a non-linear basement reservoir, without constraint on Lateral flow
The model showed an increase in lateral flow at the start of the simulation (Day 2),
however, appear to reach a pseudo steady-state of low lateral flow during the observed-data period.
This indicates that the lateral-flow option for the basement layer is not the best solution.
"""

def get_state(old_state):
    new_state = old_state + 1
    global state
    state = new_state
    return state


def start_jday():
    start_sim = 1  # 213 # 166
    return start_sim


class BeachModel(DynamicModel, MonteCarloModel):
    def setDebug(self):
        pass

    def __init__(self, cloneMap):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        setclone(cloneMap)

        # dem = self.dem

    def premcloop(self):
        self.DEBUG = False
        self.TEST_depth = False
        self.TEST_roots = False
        self.TEST_Ksat = False
        self.TEST_thProp = False
        self.TEST_theta = False
        self.TEST_IR = False
        self.TEST_PERC = False

        # Hydro
        self.LF = True
        self.ETP = True
        self.f_transp = 0.40
        self.f_evap = 0.50
        self.linear_reservoir = False

        self.PEST = True
        self.TRANSPORT = True

        self.ROM = True
        self.LCH = True
        self.LFM = True
        self.DEG = True

        self.TEST_LCH = False
        self.TEST_LFM = False
        self.TEST_DEG = False
        # This section includes all non-stochastic parameters.
        # Get initial parameters, make a dictionary of the raw file.
        import csv
        ini_path = 'initial.csv'
        self.ini_param = {}  # Dictionary to store the values
        with open(ini_path, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                self.ini_param[row[0].strip()] = float(row[1])

        date_path = 'Time.csv'
        self.time_dict = {}
        with open(date_path, 'r') as file:
            reader = csv.reader(file, delimiter=';')
            for row in reader:
                self.time_dict[str(int(row[1]))] = str(row[0].strip())

        """
        Landscape Maps
        """
        self.dem = self.readmap("dem_slope")  # 192 - 231 m a.s.l
        self.datum_depth = (self.dem - mapminimum(self.dem)) * scalar(10 ** 3)  # mm

        # self.dem_route = self.readmap("dem_ldd")  # To route surface run-off
        # self.ldd_surf = lddcreate(self.dem_route, 1e31, 1e31, 1e31, 1e31)  # To route runoff
        out_burn = readmap("dem_ldd_burn3")
        self.is_catchment = defined(out_burn)

        # self.ldd_subs = lddcreate(self.dem, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI
        self.ldd_subs = readmap('ldd_subs_v3')  # To route lateral flow & build TWI

        self.zero_map = out_burn - out_burn  # Zero map to generate scalar maps
        self.mask = out_burn / out_burn
        self.report(self.mask, 'mask')
        self.aging = deepcopy(self.zero_map)  # Cumulative days after application on each pixel

        self.outlet_multi = self.readmap("out_multi_nom_v3")  # Multi-outlet with 0 or 2
        self.is_outlet = boolean(self.outlet_multi == 1)
        self.report(self.is_outlet, 'OUT')
        self.north_wk = nominal(self.readmap("north_nom"))
        self.valley_wk = nominal(self.readmap("valley_nom"))
        self.south_wk = nominal(self.readmap("south_nom"))

        self.n1_plot = nominal(self.readmap("n1_nom"))
        self.n2_plot = nominal(self.readmap("n2_nom"))
        self.n3_plot = nominal(self.readmap("n3_nom"))
        self.n4_plot = nominal(self.readmap("n4_nom"))
        self.n5_plot = nominal(self.readmap("n5_nom"))
        self.n7_plot = nominal(self.readmap("n7_nom"))
        self.n8_plot = nominal(self.readmap("n8_nom"))

        self.t4_plot = nominal(self.readmap("t4_nom"))
        self.t5_plot = nominal(self.readmap("t5_nom"))
        self.t7_plot = nominal(self.readmap("t7_nom"))
        self.t8_plot = nominal(self.readmap("t8_nom"))
        self.t9_plot = nominal(self.readmap("t9_nom"))
        self.t10_plot = nominal(self.readmap("t10_nom"))

        self.s11_plot = nominal(self.readmap("s11_nom"))
        self.s12_plot = nominal(self.readmap("s12_nom"))
        self.s13_plot = nominal(self.readmap("s13_nom"))

        self.landuse = self.readmap("landuse2016")

        # Topographical Wetness Index
        self.up_area = accuflux(self.ldd_subs, cellarea())
        self.slope = sin(atan(max(slope(self.dem), 0.001)))  # Slope in radians
        self.wetness = ln(self.up_area / tan(self.slope))

        """
        Output & Observations (tss and observation maps)
        """
        # Output time series (tss)
        # Outlet
        ###########

        # HYDRO
        # Rain
        self.resW_accRain_m3_tss = TimeoutputTimeseries("resW_accRain_m3", self, nominal("outlet_v3"), noHeader=False)
        # Runoff
        self.resW_accRunoff_m3_tss = TimeoutputTimeseries("resW_accRunoff_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.tot_runoff_m3_tss = TimeoutputTimeseries("resW_totRunoff_m3", self, nominal("outlet_v3"), noHeader=False)
        # Percolation
        self.resW_accPercol_Bsmt_m3_tss = TimeoutputTimeseries("resW_accPercol_Bsmt_m3", self, nominal("outlet_v3"),
                                                               noHeader=False)
        # Deep percolation Basement
        self.tot_perc_z3_m3_tss = TimeoutputTimeseries("resW_totPercol_z3_m3", self, nominal("outlet_v3"),
                                                       noHeader=False)
        # ETP
        self.resW_accEtp_m3_tss = TimeoutputTimeseries("resW_accEtp_m3", self, nominal("outlet_v3"), noHeader=False)
        self.resW_accEvap_m3_tss = TimeoutputTimeseries("resW_accEvap_m3", self, nominal("outlet_v3"), noHeader=False)
        self.resW_accTransp_m3_tss = TimeoutputTimeseries("resW_accTransp_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.tot_etp_m3_tss = TimeoutputTimeseries("resW_totEtp_m3", self, nominal("outlet_v3"), noHeader=False)

        # Baseflow
        self.resW_accBaseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", self, nominal("outlet_v3"),
                                                            noHeader=False)
        # self.tot_baseflow_m3_tss = TimeoutputTimeseries("resW_totBaseflow_m3", self, nominal("outlet_v3"),
        #                                                 noHeader=False)
        # LF Drainage
        self.resW_o_accDrain_m3_tss = TimeoutputTimeseries("resW_o_accDrain_m3", self, nominal("outlet_v3"),
                                                           noHeader=False)
        self.tot_accu_drain_m3_tss = TimeoutputTimeseries("resW_o_totDrain_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)  # Cumulative ADR
        # LF options
        self.sat_accu_overflow_m3_tss = TimeoutputTimeseries("resW_of_accLatflow_m3", self, nominal("outlet_v3"),
                                                             noHeader=False)

        self.tot_accu_of_latflow_m3_tss = TimeoutputTimeseries("resW_of_totLatflow_m3", self, nominal("outlet_v3"),
                                                               noHeader=False)
        # Inflow
        self.out_cell_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_cellLatflow_m3", self, nominal("outlet_v3"),
                                                              noHeader=False)
        self.out_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_accLatflow_m3", self, nominal("outlet_v3"),
                                                              noHeader=False)
        self.tot_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_totLatflow_m3", self, nominal("outlet_v3"),
                                                              noHeader=False)
        # Outflow
        self.resW_o_cellLatflow_m3_tss = TimeoutputTimeseries("resW_o_cellLatflow_m3", self, nominal("outlet_v3"),
                                                              noHeader=False)  # Outlet LF
        # self.out_accu_o_latflow_m3_tss = TimeoutputTimeseries("resW_o_accLatflow_m3", self, nominal("outlet_v3"),
        #                                                       noHeader=False)
        self.resW_o_cumLatflow_m3_tss = TimeoutputTimeseries("resW_o_cumLatflow_m3", self, nominal("outlet_v3"),
                                                             noHeader=False)
        # self.out_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_accLatflow_m3", self, nominal("outlet_v3"),
        #                                                       noHeader=False)
        # self.tot_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_totLatflow_m3", self, nominal("outlet_v3"),
        #                                                       noHeader=False)

        self.resW_accChStorage_m3_tss = TimeoutputTimeseries("resW_accChStorage_m3", self, nominal("outlet_v3"),
                                                             noHeader=False)
        self.global_mb_water_tss = TimeoutputTimeseries("resW_global_waterMB", self, nominal("outlet_v3"),
                                                        noHeader=False)
        self.storage_m3_tss = TimeoutputTimeseries("resW_accStorage_m3", self, nominal("outlet_v3"), noHeader=False)

        # Basement layer analysis
        # self.resW_accBAL_z3_m3_tss = TimeoutputTimeseries("resW_accBAL_z3_m3", self, nominal("outlet_v3"),
        #                                                   noHeader=False)
        # self.resW_accInfil_z3_m3_tss = TimeoutputTimeseries("resW_accInfil_z3_m3", self, nominal("outlet_v3"),
        #                                                     noHeader=False)
        # self.resW_accLF_z3_m3_tss = TimeoutputTimeseries("resW_accLF_z3_m3", self, nominal("outlet_v3"),
        #                                                  noHeader=False)
        # self.resW_accETP_z3_m3_tss = TimeoutputTimeseries("resW_accETP_z3_m3", self, nominal("outlet_v3"),
        #                                                   noHeader=False)

        self.resW_accBAL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accBAL_Bsmt_m3", self, nominal("outlet_v3"),
                                                            noHeader=False)
        self.resW_accInfil_Bsmt_m3_tss = TimeoutputTimeseries("resW_accInfil_Bsmt_m3", self, nominal("outlet_v3"),
                                                              noHeader=False)
        self.resW_accLF_Bsmt_m3_tss = TimeoutputTimeseries("resW_accLF_Bsmt_m3", self, nominal("outlet_v3"),
                                                           noHeader=False)
        self.resW_accETP_Bsmt_m3_tss = TimeoutputTimeseries("resW_accETP_Bsmt_m3", self, nominal("outlet_v3"),
                                                            noHeader=False)

        # Storage
        self.resW_accVOL_z0_m3_tss = TimeoutputTimeseries("resW_accVOL_z0_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.resW_accVOL_z1_m3_tss = TimeoutputTimeseries("resW_accVOL_z1_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.resW_accVOL_z2_m3_tss = TimeoutputTimeseries("resW_accVOL_z2_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.resW_accVOL_z3_m3_tss = TimeoutputTimeseries("resW_accVOL_z3_m3", self, nominal("outlet_v3"),
                                                          noHeader=False)
        self.resW_accVOL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accVOL_Bsmt_m3", self, nominal("outlet_v3"),
                                                            noHeader=False)

        # PESTI
        # Pesticide
        self.global_mb_pest_tss = TimeoutputTimeseries("resM_global_mb_pest", self, nominal("outlet_v3"),
                                                       noHeader=False)
        self.resM_accAPP_g_tss = TimeoutputTimeseries("resM_accAPP_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accVOLAT_L_tss = TimeoutputTimeseries("resM_accVOLAT_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accRO_L_tss = TimeoutputTimeseries("resM_accRO_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accDEG_L_tss = TimeoutputTimeseries("resM_accDEG_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accDEGz0_L_tss = TimeoutputTimeseries("resM_accDEGz0_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accLCHz0_L_tss = TimeoutputTimeseries("resM_accLCH_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accDPz1_L_tss = TimeoutputTimeseries("resM_accDPz1_L", self, nominal("outlet_v3"), noHeader=False)

        self.resM_accDP_L_tss = TimeoutputTimeseries("resM_accDP_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accADR_L_tss = TimeoutputTimeseries("resM_accADR_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accLF_L_tss = TimeoutputTimeseries("resM_accLF_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accBF_L_tss = TimeoutputTimeseries("resM_accBF_L", self, nominal("outlet_v3"), noHeader=False)
        self.resM_accCHS_L_tss = TimeoutputTimeseries("resM_accCHS_L", self, nominal("outlet_v3"), noHeader=False)

        self.resM_EXP_Smet_g_tss = TimeoutputTimeseries("resM_EXP_Smet_g", self, nominal("outlet_v3"),
                                                        noHeader=False)  # Total outlet mass (g) exports (light fraction only)
        # Concentrations outlet
        self.resM_oCONC_ugL_tss = TimeoutputTimeseries("resM_oCONC_ugL", self, nominal("outlet_v3"),
                                                       noHeader=False)  # Total outlet conc (ug/L)
        self.resM_oCONC_ROFF_ugL_tss = TimeoutputTimeseries("resM_oCONC_ROFF_ugL", self, nominal("outlet_v3"),
                                                            noHeader=False)  # Runoff outlet conc (ug/L)
        self.resM_oCONC_LF_ugL_tss = TimeoutputTimeseries("resM_oCONC_LF_ugL", self, nominal("outlet_v3"),
                                                          noHeader=False)  # Latflow outlet conc (ug/L)

        self.resM_oCONC_ADR_ugL_tss = TimeoutputTimeseries("resM_oCONC_ADR_ugL", self, nominal("outlet_v3"),
                                                           noHeader=False)  # Artificial drainage outlet conc (ug/L)
        # Isotopes outlet
        self.resM_outISO_d13C_tss = TimeoutputTimeseries("resM_outISO_d13C", self, nominal("outlet_v3"),
                                                         noHeader=False)  #
        self.resM_outISO_ROFF_d13C_tss = TimeoutputTimeseries("resM_outISO_ROFF_d13C", self, nominal("outlet_v3"),
                                                              noHeader=False)  # Runoff outlet
        self.resM_outISO_LF_d13C_tss = TimeoutputTimeseries("resM_outISO_LF_d13C", self, nominal("outlet_v3"),
                                                            noHeader=False)  # Latflow outlet
        self.resM_outISO_ADR_d13C_tss = TimeoutputTimeseries("resM_outISO_ADR_d13C", self, nominal("outlet_v3"),
                                                             noHeader=False)  # Artificial drainage outlet

        # Cumulative Pesticide
        self.cum_degZ0_L_g_tss = TimeoutputTimeseries("resM_cumDEGz0_L", self, nominal("outlet_v3"),
                                                      noHeader=False)  # Deg z0
        self.resM_cumLCHz0_L_g_tss = TimeoutputTimeseries("resM_cumLCHz0_L", self, nominal("outlet_v3"),
                                                          noHeader=False)  # Leaching z0
        self.cum_roZ0_L_g_tss = TimeoutputTimeseries("resM_cumROz0_L", self, nominal("outlet_v3"),
                                                     noHeader=False)  # Runoff
        self.cum_adr_L_g_tss = TimeoutputTimeseries("resM_cumADR_L", self, nominal("outlet_v3"),
                                                    noHeader=False)  # Art. drainage
        self.cum_latflux_L_g_tss = TimeoutputTimeseries("resM_cumLF_L", self, nominal("outlet_v3"),
                                                        noHeader=False)  # Soil column, outlet cells

        self.resM_cumEXP_Smet_g_tss = TimeoutputTimeseries("resM_cumEXP_Smet_g", self, nominal("outlet_v3"),
                                                           noHeader=False)  # Total cum. outlet mass (g) exports

        self.nash_compConc_L_tss = TimeoutputTimeseries("resNash_compConc_L", self, nominal("outlet_v3"),
                                                        noHeader=False)

        # Transects
        self.north_conc_tss = TimeoutputTimeseries("resM_norCONC", self, ordinal("north_ave"), noHeader=False)
        self.valley_conc_tss = TimeoutputTimeseries("resM_valCONC", self, ordinal("valley_ave"), noHeader=False)
        self.south_conc_tss = TimeoutputTimeseries("resM_souCONC", self, ordinal("south_ave"), noHeader=False)

        self.north_d13C_tss = TimeoutputTimeseries("resM_nord13C", self, ordinal("north_ave"), noHeader=False)
        self.valley_dC13_tss = TimeoutputTimeseries("resM_vald13C", self, ordinal("valley_ave"), noHeader=False)
        self.south_d13C_tss = TimeoutputTimeseries("resM_soud13C", self, ordinal("south_ave"), noHeader=False)

        # Detailed
        self.n1_d13C_tss = TimeoutputTimeseries("resM_n1d13C", self, ordinal("n1_out"), noHeader=False)
        self.n1_conc_tss = TimeoutputTimeseries("resM_n1CONC", self, ordinal("n1_out"), noHeader=False)
        self.n2_d13C_tss = TimeoutputTimeseries("resM_n2d13C", self, ordinal("n2_out"), noHeader=False)
        self.n2_conc_tss = TimeoutputTimeseries("resM_n2CONC", self, ordinal("n2_out"), noHeader=False)
        self.n3_d13C_tss = TimeoutputTimeseries("resM_n3d13C", self, ordinal("n3_out"), noHeader=False)
        self.n3_conc_tss = TimeoutputTimeseries("resM_n3CONC", self, ordinal("n3_out"), noHeader=False)
        self.n4_d13C_tss = TimeoutputTimeseries("resM_n4d13C", self, ordinal("n4_out"), noHeader=False)
        self.n4_conc_tss = TimeoutputTimeseries("resM_n4CONC", self, ordinal("n4_out"), noHeader=False)
        self.n5_d13C_tss = TimeoutputTimeseries("resM_n5d13C", self, ordinal("n5_out"), noHeader=False)
        self.n5_conc_tss = TimeoutputTimeseries("resM_n5CONC", self, ordinal("n5_out"), noHeader=False)
        self.n7_d13C_tss = TimeoutputTimeseries("resM_n7d13C", self, ordinal("n7_out"), noHeader=False)
        self.n7_conc_tss = TimeoutputTimeseries("resM_n7CONC", self, ordinal("n7_out"), noHeader=False)
        self.n8_d13C_tss = TimeoutputTimeseries("resM_n8d13C", self, ordinal("n8_out"), noHeader=False)
        self.n8_conc_tss = TimeoutputTimeseries("resM_n8CONC", self, ordinal("n8_out"), noHeader=False)

        self.t4_d13C_tss = TimeoutputTimeseries("resM_t4d13C", self, ordinal("t4_out"), noHeader=False)
        self.t4_conc_tss = TimeoutputTimeseries("resM_t4CONC", self, ordinal("t4_out"), noHeader=False)
        self.t5_d13C_tss = TimeoutputTimeseries("resM_t5d13C", self, ordinal("t5_out"), noHeader=False)
        self.t5_conc_tss = TimeoutputTimeseries("resM_t5CONC", self, ordinal("t5_out"), noHeader=False)
        self.t7_d13C_tss = TimeoutputTimeseries("resM_t7d13C", self, ordinal("t7_out"), noHeader=False)
        self.t7_conc_tss = TimeoutputTimeseries("resM_t7CONC", self, ordinal("t7_out"), noHeader=False)
        self.t8_d13C_tss = TimeoutputTimeseries("resM_t8d13C", self, ordinal("t8_out"), noHeader=False)
        self.t8_conc_tss = TimeoutputTimeseries("resM_t8CONC", self, ordinal("t8_out"), noHeader=False)
        self.t9_d13C_tss = TimeoutputTimeseries("resM_t9d13C", self, ordinal("t9_out"), noHeader=False)
        self.t9_conc_tss = TimeoutputTimeseries("resM_t9CONC", self, ordinal("t9_out"), noHeader=False)
        self.t10_d13C_tss = TimeoutputTimeseries("resM_t10d13C", self, ordinal("t10_out"), noHeader=False)
        self.t10_conc_tss = TimeoutputTimeseries("resM_t10CONC", self, ordinal("t10_out"), noHeader=False)

        self.s11_d13C_tss = TimeoutputTimeseries("resM_s11d13C", self, ordinal("s11_out"), noHeader=False)
        self.s11_conc_tss = TimeoutputTimeseries("resM_s11CONC", self, ordinal("s11_out"), noHeader=False)
        self.s12_d13C_tss = TimeoutputTimeseries("resM_s12d13C", self, ordinal("s12_out"), noHeader=False)
        self.s12_conc_tss = TimeoutputTimeseries("resM_s12CONC", self, ordinal("s12_out"), noHeader=False)
        self.s13_d13C_tss = TimeoutputTimeseries("resM_s13d13C", self, ordinal("s13_out"), noHeader=False)
        self.s13_conc_tss = TimeoutputTimeseries("resM_s13CONC", self, ordinal("s13_out"), noHeader=False)

        # Test s11 sample variation
        self.s11_sconc_tss = TimeoutputTimeseries("resM_s11_sconc", self, ordinal("s11_ord"), noHeader=False)
        self.s11_smass_tss = TimeoutputTimeseries("resM_s11_smass", self, ordinal("s11_ord"), noHeader=False)

        # Analysis
        # This is 'q' as time series.
        self.i_Q_m3_tss = TimeoutputTimeseries("resW_i_accVol_m3", self, nominal("outlet_v3"), noHeader=False)
        self.o_Q_m3_tss = TimeoutputTimeseries("resW_o_accVol_m3", self, nominal("outlet_v3"), noHeader=False)
        self.tot_Q_m3_tss = TimeoutputTimeseries("resW_tot_accVol_m3", self, nominal("outlet_v3"), noHeader=False)
        self.q_obs_cum_tss = TimeoutputTimeseries("resW_cum_q_obs_m3", self, nominal("outlet_v3"),
                                                  noHeader=False)  # Equivalent to net_Q
        self.rain_obs_cum_tss = TimeoutputTimeseries("resW_cum_rain_obs_m3", self, nominal("outlet_v3"),
                                                     noHeader=False)  # Equivalent to net_Q
        self.rest_obs_tss = TimeoutputTimeseries("resW_q_restit_obs_m3", self, nominal("outlet_v3"),
                                                 noHeader=False)  # = rain/q_obs
        self.q_sim_cum_tss = TimeoutputTimeseries("resW_cum_q_sim_m3", self, nominal("outlet_v3"),
                                                  noHeader=False)  # Sum sim discharge (if obs available).
        self.q_sim_ave_tss = TimeoutputTimeseries("resW_q_sim_ave_m3", self, nominal("outlet_v3"),
                                                  noHeader=False)  # This is 'Nash_q' as time series.
        self.nash_q_tss = TimeoutputTimeseries("resNash_q_m3", self, nominal("outlet_v3"),
                                               noHeader=False)  # This is 'Nash_q' as time series.

        # Theta average proportion to saturation
        self.resW_z0_thetaPropSat = TimeoutputTimeseries("resW_z0_thetaPropSat", self, nominal("outlet_v3"),
                                                         noHeader=False)
        self.resW_z1_thetaPropSat = TimeoutputTimeseries("resW_z1_thetaPropSat", self, nominal("outlet_v3"),
                                                         noHeader=False)
        self.resW_z2_thetaPropSat = TimeoutputTimeseries("resW_z2_thetaPropSat", self, nominal("outlet_v3"),
                                                         noHeader=False)
        self.resW_z3_thetaPropSat = TimeoutputTimeseries("resW_z3_thetaPropSat", self, nominal("outlet_v3"),
                                                         noHeader=False)
        self.resW_Bsmt_thetaPropSat = TimeoutputTimeseries("resW_Bsmt_thetaPropSat", self, nominal("outlet_v3"),
                                                           noHeader=False)

        # Transects and detailed soils
        ###########
        # self.obs_trans = self.readmap("weekly_smp")
        # self.obs_detail = self.readmap("detailed_smp")

        # self.obs_runoff_m3_tss = TimeoutputTimeseries("obs_runoff_m3", self, ordinal("weekly_ord.map"), noHeader=False)
        # self.obs_cum_runoff_m3_tss = TimeoutputTimeseries("obs_runoff_m3_cum", self, ordinal("weekly_ord.map"), noHeader=False)
        # self.obs_latflow_m3_tss = TimeoutputTimeseries("obs_latflow_m3", self, "weekly_smp.map", noHeader=False)
        # self.obs_percol_m3_tss = TimeoutputTimeseries("obs_percol_m3", self, "weekly_smp.map", noHeader=False)
        # self.obs_etp_m3_tss = TimeoutputTimeseries("obs_etp_m3", self, "weekly_smp.map", noHeader=False)
        # self.obs_ch_storage_m3_tss = TimeoutputTimeseries("obs_chStorage_m3", self, "weekly_smp.map", noHeader=False)

    def initial(self):
        self.num_layers = int(self.ini_param.get("layers"))
        # Hydrological scenarios
        self.bsmntIsPermeable = False  # basement percolation (DP)
        self.ADLF = True

        # Morris tests
        m_state = get_state(state)  # First run will return state = 0

        """ Physical parameters for each layer """
        self.gamma = []  # coefficient to calibrate Ksat1
        self.s = []  # Not used now <- see hydro_v2.py
        self.c_lf = []
        for layer in range(self.num_layers):
            self.gamma.append(scalar(self.ini_param.get("gamma" + str(layer))))  # percolation coefficient
            self.s.append(scalar(self.ini_param.get("s" + str(layer))))  # calibrate Ksat, mm/day
            self.c_lf.append(scalar(self.ini_param.get("c" + str(layer))))  # subsurface flow coefficient)

        self.fc_adj = scalar(self.ini_param.get("fc_adj"))  # Adjusting Paul's FC (equivalent so far)
        self.root_adj = scalar(self.ini_param.get("root_adj"))  # Adjust Root Depth factor
        self.c_adr = scalar(self.ini_param.get("c_adr"))
        self.drainage_layers = [False, False, True, False, False]  # z2 turned on!!
        z3_factor = scalar(self.ini_param.get("z3_factor"))  # Represents top proportion of bottom-most layer
        if m_state == 0:
            pass
        elif m_state == 1:
            # for layer in range(self.num_layers):
            #     self.c_lf[layer] = 0.50
            # z3_factor = scalar(0.7)
            self.c_adr = 0.10
            # self.root_adj *= 0.7
            # epsilon = -2 * 2.369  # high deg
            # self.bsmntIsPermeable = True
            # self.gamma[3] = 0.02
        elif m_state == 2:
            # for layer in range(self.num_layers):
            #     self.c_lf[layer] = 2
            # z3_factor = scalar(0.4)
            self.c_adr = 0.05
            # self.root_adj *= 0.5
            # epsilon = -2 * 2.476  # mid deg
            # self.gamma[1] *= 0.5
            # self.bsmntIsPermeable = True
            # self.gamma[3] = 0.1

        self.k = scalar(self.ini_param.get("k"))  # coefficient of declining LAI in end stage

        """ Soil Properties """
        self.p_b = scalar(self.ini_param.get("p_b"))  # Soil bulk density (g/cm^3)
        self.f_oc = scalar(self.ini_param.get("f_oc"))  # Organic Carbon in soil without grass (kg/kg)

        """
        Sorption parameters
        """
        # K_oc - S-metolachlor (K_oc in ml/g)
        # Marie's thesis: log(K_oc)= 2.8-1.6 [-] -> k_oc = 63 - 398 (Alletto et al., 2013).
        # Pesticide Properties Database: k_oc = 120 ml/g (range: 50-540 mL/g)
        self.k_oc = scalar(self.ini_param.get("k_oc"))  # ml/g
        self.k_d = self.k_oc * self.f_oc  # Dissociation coefficient K_d (mL/g = L/Kg)

        # Pesticide Properties Database states :
        # K_d=0.67; but here
        # K_d=120*0.021 = 1.52;  (L/kg)
        # Difference will lead to higher retardation factor (more sorption)

        """
        Volatilization parameters
        """
        #
        # Henry's constant @ 20 C (Metolachlor, Feigenbrugel et al., 2004)
        self.k_cp = max(scalar(self.ini_param.get("k_cp")), scalar(42600.00))  # mol/L atm
        # Henry, dimensionless conversion Hcc = Hcp*R*T
        self.k_h = self.k_cp * 0.0821 * 273.15  # scalar(self.ini_param.get("k_h"))
        self.molar = scalar(self.ini_param.get("molar"))  # g S-met/mol

        """
        Degradation parameters
        """
        # Half-lives
        # DT50 (typical) = 90
        # DT50 (lab at 20°C) = 15 (range: 8-38 days)
        # DT50 (field) = 21 (range 11-31 Switzerland and France)
        # DT50 (foliar) = 5
        # DT50 (water-sediment) = 365
        # DT50 (water phase only) = 88
        self.dt_50_ref = scalar(self.ini_param.get("dt_50_ref"))  # S-met (days)

        self.temp_ref = scalar(self.ini_param.get("temp_ref"))  # Temp.  reference

        self.beta_moisture = scalar(self.ini_param.get(
            "beta_moisture"))  # Need to find a correct value for 'B', exponent in moisture dependency (Dairon)
        self.alpha_temperature = scalar(self.ini_param.get(
            "alpha_temperature"))  # Need to confirm units Ea = 54000 KJ/mol; R = 8.314 J/mol/Kelvin

        self.act_e = scalar(self.ini_param.get("activation_e"))
        self.r_gas = scalar(self.ini_param.get("r_gas"))

        """
        Isotopes
        """
        self.r_standard = scalar(self.ini_param.get("r_standard"))  # VPDB
        # self.alpha_iso = scalar(self.ini_param.get("alpha_iso"))  # 2 is no fractionation

        # Degradation Scenarios
        epsilon = -1 * 1.743  # low deg
        self.alpha_iso = epsilon / 1000 + 1

        """
        Layer depths
        """
        # Lowest layer depth and baseflow
        if self.linear_reservoir:
            self.k_g = 40000  # [days]
        self.gw_factor = 1 - z3_factor  # Represents bottom-most portion of bottom layer

        self.layer_depth = []
        self.tot_depth = deepcopy(self.zero_map)
        bottom_depth = deepcopy(self.zero_map)
        for layer in range(self.num_layers):
            if layer < self.num_layers - 2:  # 5 - 1 = 3 (i.e. z0,z1,z2)
                self.layer_depth.append(self.zero_map +
                                        scalar(self.ini_param.get('z' + str(layer))))
                self.tot_depth += self.layer_depth[layer]
                self.report(self.layer_depth[layer], 'DepthZ' + str(layer))
            elif layer < self.num_layers - 1:  # 5 - 2 = 4 (i.e. z3)
                bottom_depth = (self.datum_depth +  # total height
                                scalar(self.ini_param.get('z' + str(layer))) + 100  # plus a min-depth
                                - self.tot_depth)
                self.layer_depth.append(bottom_depth * z3_factor)  # minus:(z0, z1, z2)*decreasing depth factor
                self.tot_depth += self.layer_depth[layer]
                self.report(self.layer_depth[layer], 'DepthZ' + str(layer))
            else:  # Basement Layer = n5  (z4)
                self.layer_depth.append(bottom_depth * self.gw_factor)  # minus:(z0, z1, ...)*decreasing depth factor
                self.tot_depth += self.layer_depth[layer]
                self.report(self.layer_depth[layer], 'DepthZ' + str(layer))

            if self.TEST_depth:
                checkLayerDepths(self, layer)

        self.smp_depth = self.layer_depth[0]
        self.report(self.tot_depth, 'zTot_mm')
        #  aguila --scenarios='{2}' DepthZ0 DepthZ1 DepthZ2 DepthZ3 DepthZ4 zTot_mm

        """
        Hydro Maps
        """
        # Initial moisture (Final from model v1, Sept 30, 2016)
        self.theta = []
        if start_jday() < 166:
            self.theta.append(readmap('d14_theta_z0'))  # map of initial soil moisture in top layer (-)
            self.theta.append(readmap('d14_theta_z1'))  # d14_theta_z1
            self.theta.append(readmap('d14_theta_z2'))
            self.theta.append(readmap("d14_theta_z3"))  # * z3_factor + scalar(0.6) * self.gw_factor
            self.theta.append(readmap("d14_theta_z4"))  # * z3_factor + scalar(0.6) * self.gw_factor
        else:
            self.theta.append(readmap('d166_theta_z0'))  # map of initial soil moisture in top layer (-)
            self.theta.append(readmap('d166_theta_z1'))
            self.theta.append(readmap('d166_theta_z2'))
            self.theta.append(readmap("d166_theta_z3"))
            self.theta.append(readmap("d166_theta_z4"))

        # Need initial states to compute change in storage after each run
        self.theta_ini = deepcopy(self.theta)

        """
        Pesticides Maps
        """
        # Application days
        self.app_days = [177, 196, 238]
        # self.app_days = [177, 180, 183] # test days
        self.aged_days = ifthen(boolean(self.is_catchment), scalar(365))

        # Mass
        # in ug = conc. (ug/g soil) * density (g/cm3) * (10^6 cm3/m3)*(2 m/10^3 mm)* depth_layer(mm) * cellarea(m2)
        # g = ug * 1e-06
        self.sm_background = []
        mean_back_conc = [0.06, 0.03, 0.001, 0.001, 0.001]
        for layer in range(self.num_layers):
            background = ((self.zero_map + mean_back_conc[layer]) * self.p_b * scalar(10 ** 6 / 10 ** 3) *
                          self.layer_depth[layer] * cellarea() * (10 ** -6))  # Based on detailed soils
            self.sm_background.append(background)

        # Fraction masses and Delta (Background)
        self.lightmass = []
        self.lightmass_ini = []
        self.heavymass = []
        self.heavymass_ini = []

        self.delta = []
        self.delta_ini = []
        self.light_back = []
        self.heavy_back = []

        # Initial Isotope Signature
        for layer in range(self.num_layers):
            # Initial deltas assume theoretical max @99% deg Streitwieser Semiclassical Limits
            self.delta.append(self.zero_map - 23.7)
            self.delta_ini.append(self.zero_map - 23.7)
            self.light_back.append(self.sm_background[layer] /
                                   (1 + self.r_standard * (self.delta[layer] / 1000 + 1)))
            self.heavy_back.append(self.sm_background[layer] - self.light_back[layer])

            # Set mass fractions <- background fractions
            self.lightmass.append(deepcopy(self.light_back[layer]))
            self.lightmass_ini.append(deepcopy(self.light_back[layer]))
            self.heavymass.append(deepcopy(self.heavy_back[layer]))
            self.heavymass_ini.append(deepcopy(self.heavy_back[layer]))

            if mapminimum(self.lightmass[layer]) < 0:
                print("Err INI, light")

            if mapminimum(self.heavymass[layer]) < 0:
                print("Err INI, heavy")

        # Assign dosages based on Farmer-Crop combinations [g/m2]
        self.fa_cr = readmap("farm_burn_v3")  # Contains codes to assign appropriate dosage
        self.apps = getApplications(self, self.fa_cr, massunit='g')  # returns list of applied masses

        # Applications delta
        # Use map algebra to produce a initial signature map,
        # ATT: Need to do mass balance on addition of new layer.
        # where app1 > 0, else background sig. (plots with no new mass will be 0)
        # where app1 > 0, else background sig. (plots with no new mass will be 0)
        self.appDelta = []
        self.appDelta.append(ifthenelse(self.apps[0] > 0, scalar(-32.3), scalar(-23.7)))
        self.appDelta.append(ifthenelse(self.apps[1] > 0, scalar(-32.3), scalar(-23.7)))
        self.appDelta.append(ifthenelse(self.apps[2] > 0, scalar(-32.3), scalar(-23.7)))

        # Cumulative maps
        # self.light_ini_storage_ug = self.lightmass_z0_ini + self.lightmass_z1_ini + self.lightmass_z2_ini
        # self.cum_appl_g = self.zero_map
        self.cum_runoff_ug = deepcopy(self.zero_map)
        self.cum_leached_ug_z0 = deepcopy(self.zero_map)
        self.cum_leached_ug_z1 = deepcopy(self.zero_map)
        self.cum_leached_ug_z2 = deepcopy(self.zero_map)
        self.cum_leached_ug_z3 = deepcopy(self.zero_map)

        self.cum_latflux_ug_z0 = deepcopy(self.zero_map)
        self.cum_latflux_ug_z1 = deepcopy(self.zero_map)
        self.cum_latflux_ug_z2 = deepcopy(self.zero_map)
        self.cum_latflux_ug_z3 = deepcopy(self.zero_map)

        self.cum_baseflx_ug_z3 = deepcopy(self.zero_map)

        self.cum_degZ0_L_g = deepcopy(self.zero_map)
        self.cum_roZ0_L_g = deepcopy(self.zero_map)
        self.cum_lchZ0_L_g = deepcopy(self.zero_map)
        self.cum_adr_L_g = deepcopy(self.zero_map)
        self.cum_latflux_L_g = deepcopy(self.zero_map)
        self.cum_exp_L_g = deepcopy(self.zero_map)
        self.northConc_diff = deepcopy(self.zero_map)
        self.northConc_var = deepcopy(self.zero_map)
        self.valleyConc_diff = deepcopy(self.zero_map)
        self.valleyConc_var = deepcopy(self.zero_map)
        self.southConc_diff = deepcopy(self.zero_map)
        self.southConc_var = deepcopy(self.zero_map)

        """
        Temperature maps and params
        """
        self.lag = scalar(0.8)  # lag coefficient (-), 0 < lag < 2; -> in SWAT, lag = 0.80
        # Generating initial surface temp map (15 deg is arbitrary)
        self.temp_fin = []
        self.temp_surf_fin = self.zero_map + 15
        for layer in range(self.num_layers):
            self.temp_fin.append(self.zero_map + 15)

        # Maximum damping depth (dd_max)
        # The damping depth (dd) is calculated daily and is a function of max. damping depth (dd_max), (mm):
        self.dd_max = (scalar(2500) * self.p_b) / (self.p_b + 686 * exp(-5.63 * self.p_b))

        # TODO
        # Average Annual air temperature (celcius - Layon!! Not Alteckendorf yet!!)
        self.temp_ave_air = scalar(12.2)  # 12.1 is for Layon

        """
        Simulation start time
        """
        start_day = start_jday()  # Returns initial timestep
        greg_date = self.time_dict[str(start_day)].split("/", 2)
        print("Date: ", greg_date[0], greg_date[1], greg_date[2])
        print("Sim Day: ", start_day)
        yy = int(greg_date[2])
        mm = int(greg_date[1])
        dd = int(greg_date[0])

        date_factor = 1
        if (100 * yy + mm - 190002.5) < 0:
            date_factor = -1

        # simulation start time in JD (Julian Day)
        self.jd_start = 367 * yy - rounddown(7 * (yy + rounddown((mm + 9) / 12)) / 4) + rounddown(
            (275 * mm) / 9) + dd + 1721013.5 - 0.5 * date_factor
        self.jd_cum = 0
        self.jd_dt = 1  # Time step size (days)

        # Analysis
        self.water_balance = []  # mm
        for layer in range(self.num_layers):
            self.water_balance.append(deepcopy(self.zero_map))
        self.days_cum = 0  # Track no. of days with data
        self.q_diff = 0
        self.q_var = 0
        self.q_obs_cum = 0
        self.q_sim_cum = 0  # net total disch
        self.q_sim_ave = 0

        self.rain_cum_m3 = 0  # Rainfall
        self.rain_cum_mm = self.zero_map + scalar(400.0)  # Cum Rainfall

        self.tot_drain_m3 = 0  # drainage z1
        # self.tot_nlf_m3 = 0
        # self.tot_ilf_m3 = 0  # upstream inflow
        self.cum_olf_m3 = 0  # downstream outflow

        self.tot_of_m3 = 0  # Overflow due to LF sat capacity reached

        self.tot_etp_m3 = 0
        self.tot_baseflow_m3 = 0
        self.tot_perc_z3_m3 = 0
        self.tot_runoff_m3 = 0

        # Stochastic / test parameters
        print("state:", m_state)
        self.theta_wp = scalar(self.ini_param.get("wp_zAll"))  # => 0.19
        theta_sat_z2 = self.zero_map + readmap("thSATz2")  # => 0.63
        # scalar(self.ini_param.get("sat_z2z3")) + mapnormal() * 0.04  # mean + 1SD(X)*0.04 = mean + (0.002)**0.5
        theta_fcap_z2 = self.zero_map + readmap("thFCz2")  # => 0.39
        # scalar(self.ini_param.get("fc_z2z3")) + mapnormal() * 0.04  # mean + 1SD(X)*0.04 = mean + (0.002)**0.5
        # self.report(self.theta_sat_z2, "thSATz2")
        # self.report(self.theta_fcap_z2, "thFCz2")

        self.theta_sat = []
        self.theta_fc = []
        for i in range(self.num_layers):
            if i < 2:
                self.theta_sat.append(deepcopy(self.zero_map))
                self.theta_fc.append(deepcopy(self.zero_map))
            else:
                self.theta_sat.append(deepcopy(theta_sat_z2))
                self.theta_fc.append(deepcopy(theta_fcap_z2))

                # Increase bottom layer's initial moisture by 20% saturation
                # self.theta[-2] += self.theta_sat[-2] * .10
                # self.theta[-2] = min(deepcopy(self.theta_sat[-2]), self.theta[-2])

    def dynamic(self):

        jd_sim = self.jd_start + self.jd_cum
        if self.PEST:
            self.aged_days += scalar(1)
        # timeinputscalar() gets the TSS's cell value of row (timestep) and TSS's column indexed by the landuse-map.
        # In other words, the value of the landuse-map pixel == column to to look for in landuse.tss
        # So currently becasue landuse does not change value in the year, this step is redundant
        # and we could simply use the landuse map to map the fields to the "Crop Parameters" below.
        # Mapping "landuse.map" to -> "fields map" (i.e. the latter is a dyanmic-landuse equivalent).
        fields = timeinputscalar('landuse.tss', nominal(self.landuse))  #
        # Note that the number of columns could still be reduced to 9 as, only 9 classes are considered in 2016.

        " Crop Parameters "
        # SEE: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_lookup.html?highlight=lookupscalar
        setglobaloption('matrixtable')  # allows lookupscalar to read more than 1 expressions.
        crop_type = lookupscalar('croptable.tbl', 1, fields)  # (table, col-value in crop table, column-value in fields)
        sow_yy = lookupscalar('croptable.tbl', 2, fields)
        sow_mm = lookupscalar('croptable.tbl', 3, fields)  # sowing or Greenup month
        sow_dd = lookupscalar('croptable.tbl', 4, fields)  # sowing day
        sow_dd = ifthenelse(self.fa_cr == 1111, sow_dd - 15,  # Beet Friess
                            sow_dd)
        len_grow_stage_ini = lookupscalar('croptable.tbl', 5,
                                          fields)  # old: Lini. length of initial crop growth stage
        len_dev_stage = lookupscalar('croptable.tbl', 6, fields)  # Ldev: length of development stage
        len_mid_stage = lookupscalar('croptable.tbl', 7, fields)  # Lmid: length of mid-season stage
        len_end_stage = lookupscalar('croptable.tbl', 8, fields)  # Lend: length of late season stage
        kcb_ini = lookupscalar('croptable.tbl', 9, fields)  # basal crop coefficient at initial stage
        kcb_mid = lookupscalar('croptable.tbl', 10, fields)  # basal crop coefficient at mid season stage
        kcb_end = lookupscalar('croptable.tbl', 11, fields)  # basal crop coefficient at late season stage
        max_LAI = lookupscalar('croptable.tbl', 12, fields)  # maximum leaf area index
        mu = lookupscalar('croptable.tbl', 13, fields)  # light use efficiency
        max_height = lookupscalar('croptable.tbl', 14, fields)  # maximum crop height

        max_root_depth = lookupscalar('croptable.tbl', 15, fields) * 1000  # max root depth converting from m to mm
        # Max RD (m) according to Allen 1998, Table 22 (now using FAO source)
        # Sugar beet = 0.7 - 2.1
        # Corn = 2.0 - 2.7
        # Grazing pasture 0.5 - 2.5
        # Spring Wheat = 2.0 -2.5
        # Winter Wheat = 2.5 -2.8
        # Apple trees = 2.0-1.0

        # depletable theta before water stress (Allen1998, Table no.22)
        p_tab = lookupscalar('croptable.tbl', 16, fields)
        # Sugar beet = 0.55
        # Corn = 0.55
        # Grazing Pasture = 0.6
        # Spring Wheat = 0.55
        # Winter Wheat = 0.55
        # Apple trees = 0.5

        """ Soil physical parameters
        # Moved theta properties to initial(), for Morris test.
        """
        # Basement layers defined under initial()
        self.theta_sat[0] = timeinputscalar('thetaSat_agr.tss', nominal(self.landuse))  # saturated moisture # [-]
        self.theta_sat[1] = deepcopy(self.theta_sat[0])
        self.theta_fc[0] = timeinputscalar('thetaFC_agr.tss', nominal(self.landuse)) * self.fc_adj  # field capacity
        self.theta_fc[1] = deepcopy(self.theta_fc[0])
        theta_wp = self.theta_wp  # scalar(0.19)  # lookupscalar('croptable.tbl', 21, fields)  # wilting point moisture

        if self.TEST_thProp:
            checkMoistureProps(self, self.theta_sat, 'aSATz')
            checkMoistureProps(self, self.theta_fc, 'aFCz')

        k_sat_z0z1 = timeinputscalar('ksats.tss', nominal(self.landuse))
        k_sat_z2z3 = lookupscalar('croptable.tbl', 17, fields)  # saturated conductivity of the second layer
        k_sat = []
        for i in range(self.num_layers):
            if i < 2:
                k_sat.append(deepcopy(k_sat_z0z1))
            else:
                k_sat.append(deepcopy(k_sat_z2z3))

        CN2_A = lookupscalar('croptable.tbl', 18, fields)  # curve number of moisture condition II
        CN2_B = lookupscalar('croptable.tbl', 19, fields)  # curve number of moisture condition II
        CN2_C = lookupscalar('croptable.tbl', 20, fields)  # curve number of moisture condition II
        CN2_D = lookupscalar('croptable.tbl', 21, fields)  # curve number of moisture condition II

        if self.TEST_Ksat:
            reportKsatEvolution(self, k_sat)

        """
        Time-series data to spatial location,
        map is implicitly defined as the clonemap.
        """
        precip = timeinputscalar('rain.tss', 1)  # daily precipitation data as time series (mm)
        # Precipitation total
        rain_m3 = self.mask * precip * cellarea() / 1000  # m3
        tot_rain_m3 = areatotal(rain_m3, self.is_catchment)

        temp_bare_soil = timeinputscalar('T_bare.tss', nominal('clone_nom'))  # SWAT, Neitsch2009, p.43.
        self.temp_air = timeinputscalar('airTemp.tss', nominal('clone_nom'))
        et0 = timeinputscalar('ET0.tss', 1)  # daily ref. ETP at Zorn station (mm)
        wind = timeinputscalar('U2.tss', 1)  # wind speed time-series at 1 meters height
        humid = timeinputscalar('RHmin.tss', 1)  # minimum relative humidity time-series # PA: (-)
        # precipVol = precip * cellarea() / 1000  # m3

        ################
        # Crop growth ##
        ################
        jd_sow = convertJulian(sow_yy, sow_mm, sow_dd)
        self.rain_cum_mm += precip
        self.rain_cum_mm = ifthenelse(jd_sim == jd_sow, scalar(0), self.rain_cum_mm)
        # self.report(self.rain_cum_mm, 'aCuRain')
        CN2 = ifthenelse(self.rain_cum_mm > 60, CN2_D, CN2_A)  # report_CN(self, CN2, self.rain_cum_mm)

        # soil_group = ifthenelse(self.rain_cum_mm > 90, ordinal(3), ordinal(2))
        all_stages = len_grow_stage_ini + len_dev_stage + len_mid_stage + len_end_stage

        # updating of sowing date by land use
        # TODO: Why update if sow date is already by land use?
        # sow_yy = ifthenelse(jd_sim < jd_sow + all_stages, sow_yy,
        #                     ifthenelse(jd_sim < jd_sow + all_stages + 365, sow_yy + 2,
        #                                ifthenelse(jd_sim < jd_sow + all_stages + 730, sow_yy + 1,
        #                                           ifthenelse(jd_sim < jd_sow + all_stages + 1095, sow_yy + 3,
        #                                                      ifthenelse(jd_sim < jd_sow + all_stages + 1460,
        #                                                                 sow_yy + 4,
        #                                                                 scalar(0))))))

        # Update sowing date / plant date
        jd_plant = convertJulian(sow_yy, sow_mm, sow_dd)

        jd_dev = jd_plant + len_grow_stage_ini
        jd_mid = jd_dev + len_dev_stage
        jd_late = jd_mid + len_mid_stage
        jd_end = jd_late + len_end_stage
        LAIful = max_LAI + 0.5

        height = timeinputscalar('height.tss', nominal(self.landuse))
        root_depth_tot = timeinputscalar('height.tss', nominal(self.landuse)) * self.root_adj
        root_depth_tot *= 10 ** 3  # Convert to mm

        root_depth = []
        for layer in range(self.num_layers):
            if layer == 0:
                root_depth_z0 = ifthenelse(root_depth_tot > self.layer_depth[0], self.layer_depth[0], root_depth_tot)
                root_depth.append(root_depth_z0)
            elif layer == 1:
                root_depth_z1 = ifthenelse(root_depth_tot < self.layer_depth[0], scalar(0),
                                           ifthenelse(root_depth_tot <= self.layer_depth[1] + self.layer_depth[0],
                                                      root_depth_tot - self.layer_depth[0], self.layer_depth[1]))
                root_depth.append(root_depth_z1)
            elif layer == 2:
                root_depth_z2 = ifthenelse(root_depth_tot <= self.layer_depth[0] + self.layer_depth[1], scalar(0),
                                           ifthenelse(
                                               root_depth_tot <= self.layer_depth[0] + self.layer_depth[1] +
                                               self.layer_depth[2],
                                               root_depth_tot - self.layer_depth[0] - self.layer_depth[1],
                                               self.layer_depth[2]))
                root_depth.append(root_depth_z2)
            elif layer == 3:
                root_depth_z3 = ifthenelse(
                    root_depth_tot <= self.layer_depth[0] + self.layer_depth[1] + self.layer_depth[2],
                    scalar(0),
                    ifthenelse(
                        root_depth_tot <= self.layer_depth[0] + self.layer_depth[1] + self.layer_depth[2] +
                        self.layer_depth[3],
                        root_depth_tot - self.layer_depth[0] - self.layer_depth[1] - self.layer_depth[2],
                        self.layer_depth[3]))
                root_depth.append(root_depth_z3)
            else:
                root_depth.append(scalar(0))

        if self.TEST_roots:
            checkRootDepths(root_depth)
        # calculation of leaf area index
        LAI = ifthenelse(jd_sim < jd_plant, scalar(0),
                         ifthenelse(jd_sim < jd_mid,
                                    max_LAI * (jd_sim - jd_plant) / (len_grow_stage_ini + len_dev_stage),
                                    ifthenelse(jd_sim < jd_mid + 0.5 * len_mid_stage,
                                               max_LAI + (LAIful - max_LAI) * (jd_sim - jd_mid) / (
                                                   0.5 * len_mid_stage),
                                               ifthenelse(jd_sim < jd_late, LAIful,
                                                          ifthenelse(jd_sim <= jd_end,
                                                                     LAIful * exp(
                                                                         -self.k * (jd_sim - jd_late)),
                                                                     scalar(0))))))
        # calculation of fraction of soil covered by vegetation
        # frac_soil_cover = 2 - exp(-mu * LAI)
        # \mu is a light-use efficiency parameter that
        # depends on land-use characteristics
        # (i.e. Grass: 0.35; Crops: 0.45; Trees: 0.5-0.77; cite: Larcher, 1975).

        # TODO: Check "f" definition by Allan et al., 1998 against previous (above)
        # fraction of soil cover is calculated inside the "getPotET" function.
        # frac_soil_cover = ((Kcb - Kcmin)/(Kcmax - Kcmin))**(2+0.5*mean_height)
        # self.fTss.sample(frac_soil_cover)

        # Get potential evapotranspiration for all layers
        etp_dict = getPotET(self, sow_yy, sow_mm, sow_dd,
                            jd_sim,
                            wind, humid,
                            et0,
                            kcb_ini, kcb_mid, kcb_end,
                            height,
                            len_grow_stage_ini, len_dev_stage, len_mid_stage, len_end_stage,
                            p_tab)
        pot_transpir = etp_dict["Tp"]
        pot_evapor = etp_dict["Ep"]
        depletable_water = etp_dict["P"]
        # TODO: printouts!
        # self.report(pot_transpir, 'aPotTRA')
        # self.report(pot_evapor, 'aPotEVA')

        # Not in use for water balance, but used to estimate surface temp due to bio-cover.
        frac_soil_cover = etp_dict["f"]
        # self.report(frac_soil_cover, 'aFracCV')

        bio_cover = getBiomassCover(self, frac_soil_cover)
        # bcv should range 0 (bare soil) to 2 (complete cover)
        # self.report(bio_cover, 'aBCV')

        # Applications
        mass_applied = deepcopy(self.zero_map)
        light_applied = deepcopy(self.zero_map)
        heavy_applied = deepcopy(self.zero_map)
        light_volat = deepcopy(self.zero_map)

        """
        Applications and Volatilization (on application days only)
        """
        if self.currentTimeStep() in self.app_days:
            mass_applied = ifthenelse(self.currentTimeStep() == self.app_days[0],
                                      self.apps[0],
                                      ifthenelse(self.currentTimeStep() == self.app_days[1], self.apps[1],
                                                 ifthenelse(self.currentTimeStep() == self.app_days[2], self.apps[2],
                                                            scalar(0))))

            self.aged_days = ifthenelse(mass_applied > 0, scalar(0), self.aged_days)

            light_applied = ifthenelse(self.currentTimeStep() == self.app_days[0],
                                       mass_applied / (1 + self.r_standard * (self.appDelta[0] / 1000 + 1)),
                                       ifthenelse(self.currentTimeStep() == self.app_days[1],
                                                  mass_applied / (
                                                      1 + self.r_standard * (self.appDelta[1] / 1000 + 1)),
                                                  ifthenelse(self.currentTimeStep() == self.app_days[2],
                                                             mass_applied / (
                                                                 1 + self.r_standard * (self.appDelta[2] / 1000 + 1)),
                                                             scalar(0))))
            heavy_applied = mass_applied - light_applied

            self.lightmass[0] += light_applied
            self.heavymass[0] += heavy_applied
            if mapminimum(self.lightmass[0]) < 0:
                print("Err APP, light")

            if mapminimum(self.heavymass[0]) < 0:
                print("Err APP, heavy")

            if self.TRANSPORT:
                # Mass volatilized
                light_volat = getVolatileMass(self, self.temp_air, self.lightmass[0],  # "LL",
                                              rel_diff_model='option-2', sorption_model="linear",
                                              gas=True, run=self.PEST)
                heavy_volat = getVolatileMass(self, self.temp_air, self.heavymass[0],  # "HH",
                                              rel_diff_model='option-2', sorption_model="linear",
                                              gas=True, run=self.PEST)

                light_volat = ifthenelse(mass_applied > 0, light_volat, scalar(0))
                heavy_volat = ifthenelse(mass_applied > 0, heavy_volat, scalar(0))

                self.lightmass[0] -= light_volat
                self.heavymass[0] -= heavy_volat
                if mapminimum(self.lightmass[0]) < 0:
                    print("Err Volat, light")
                    self.report(self.lightmass[0], 'VOL_Mz0')

                if mapminimum(self.heavymass[0]) < 0:
                    print("Err Volat, heavy")

        """
        Infiltration, runoff, & percolation (all layers)
        """
        infil_z0 = deepcopy(self.zero_map)
        runoff_z0 = deepcopy(self.zero_map)
        percolation = []
        mass_runoff = []  # 0 <- light, 2 <- heavy
        light_leached = []
        heavy_leached = []
        permeable = True
        for layer in range(self.num_layers):
            if layer == 0:  # Layer 0
                if mapminimum(self.lightmass[0]) < 0:
                    print("Err Start, light")
                if mapminimum(self.heavymass[0]) < 0:
                    print("Err Start, heavy")

                z0_IRO = getTopLayerInfil(self, precip, theta_wp, CN2, crop_type,
                                          jd_sim, jd_dev, jd_mid, jd_end, len_dev_stage)
                runoff_z0 = z0_IRO.get("roff")  # [mm]
                # Partition infiltration
                infil_z0 = z0_IRO.get("infil_z0")  # [mm]
                infil_z1 = z0_IRO.get("infil_z1")  # [mm]

                # Distribute to layer 0
                SW0 = self.theta[layer] * self.layer_depth[layer] + infil_z0
                self.theta[layer] = SW0 / self.layer_depth[layer]  # [-]

                # Distribution to layer 2 needed here, bc. getPercolation(layer = 0) will check below's capacity
                SW1 = self.theta[layer + 1] * self.layer_depth[layer + 1] + infil_z1
                self.theta[layer + 1] = SW1 / self.layer_depth[layer + 1]  # [-]

                if mapmaximum(self.theta[layer]) > mapmaximum(self.theta_sat[layer]):
                    val = float(mapmaximum(self.theta[layer])) - float(mapmaximum(self.theta_sat[layer]))
                    self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                                   self.theta[layer])
                    if float(val) > float(1e-06):
                        print("Error at Percolation(), SAT exceeded, layer " + str(layer) + ' by ' + str(val))

                if mapmaximum(self.theta[layer + 1]) > mapmaximum(self.theta_sat[layer + 1]):
                    val = float(mapmaximum(self.theta[layer + 1])) - float(mapmaximum(self.theta_sat[layer + 1]))
                    self.theta[layer + 1] = ifthenelse(self.theta[layer + 1] > self.theta_sat[layer + 1],
                                                       self.theta_sat[layer + 1],
                                                       self.theta[layer + 1])
                    if float(val) > float(1e-06):
                        print("Error at Percolation(), SAT exceeded, layer " + str(layer + 1))

                # RunOff Mass
                # Mass & delta run-off (RO)
                mass_runoff.append(getRunOffMass(self, precip, runoff_z0, self.lightmass[layer],
                                                 transfer_model="nu-mlm", sorption_model="linear",
                                                 gas=True, run=self.ROM))
                mass_runoff.append(getRunOffMass(self, precip, runoff_z0, self.heavymass[layer],
                                                 transfer_model="nu-mlm", sorption_model="linear",
                                                 gas=True, run=self.ROM))
                self.lightmass[layer] -= mass_runoff[0]  # light
                self.heavymass[layer] -= mass_runoff[1]  # heavy
                if mapminimum(self.lightmass[layer]) < 0:
                    print("Err RO, light")

                if mapminimum(self.heavymass[layer]) < 0:
                    print("Err RO, heavy")

                percolation.append(getPercolation(self, layer, k_sat[layer], isPermeable=permeable))  # [mm]
                water_flux = infil_z0 + infil_z1 + percolation[layer]

                light_leached.append(getLeachedMass(self, layer, water_flux, self.lightmass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.LCH))
                heavy_leached.append(getLeachedMass(self, layer, water_flux, self.heavymass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.LCH))

                SW1b = self.theta[layer] * self.layer_depth[layer] - percolation[layer]
                self.theta[layer] = SW1b / self.layer_depth[layer]
                self.lightmass[layer] -= light_leached[layer]
                self.heavymass[layer] -= heavy_leached[layer]

                if mapminimum(self.lightmass[layer]) < 0:
                    print("Err LCH, light")
                if mapminimum(self.heavymass[layer]) < 0:
                    print("Err LCH, heavy")

                # Discharge due to runoff at the outlet
                runoff_m3 = runoff_z0 * cellarea() / 1000  # m3
                accu_runoff_m3 = accuflux(self.ldd_subs, runoff_m3)
                out_runoff_m3 = areatotal(accu_runoff_m3, self.outlet_multi)

            else:  # Layers 2, 1, 3 & 4
                if layer == (self.num_layers - 1):
                    permeable = self.bsmntIsPermeable

                if mapminimum(self.lightmass[layer]) < 0:
                    print("Err Startz1, light")
                if mapminimum(self.heavymass[layer]) < 0:
                    print("Err Startz1, heavy")

                if mapmaximum(self.theta[layer]) > mapmaximum(self.theta_sat[layer]):
                    val = float(mapmaximum(self.theta[layer])) - float(mapmaximum(self.theta_sat[layer]))
                    if float(val) > float(1e-06):
                        print("Error at right before Percolation() layers: 2, 1, 3, SAT exceeded, layer " + str(layer))
                    self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                                   self.theta[layer])

                SW2 = self.theta[layer] * self.layer_depth[layer] + percolation[layer - 1]
                self.theta[layer] = SW2 / self.layer_depth[layer]
                # exceed = max(self.theta[layer] - self.theta_sat[layer], scalar(0))
                # self.report(exceed, 'outEXz' + str(layer))

                # recordInfiltration(self, percolation[layer - 2], layer)
                # self.report(SW, 'SWz' + str(layer))
                # if self.TEST_theta:
                #     checkMoisture(self, self.theta, 'athz')

                self.lightmass[layer] += light_leached[layer - 1]
                self.heavymass[layer] += heavy_leached[layer - 1]

                if mapmaximum(self.theta[layer]) > mapmaximum(self.theta_sat[layer]):
                    val = float(mapmaximum(self.theta[layer])) - float(mapmaximum(self.theta_sat[layer]))
                    if float(val) > float(1e-06):
                        print(
                            "Error at Percolation() layers: 2, 1, 3, SAT exceeded, layer " + str(layer) + ' by ' + str(
                                val))
                    self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                                   self.theta[layer])

                percolation.append(getPercolation(self, layer, k_sat[layer], isPermeable=permeable))

                if layer < (len(self.layer_depth) - 1):  # layers: 0,2,1,3
                    sw_check_bottom = self.theta[layer + 1] * self.layer_depth[layer + 1] + percolation[layer]
                    exceed_mm = max(sw_check_bottom - self.theta_sat[layer + 1] * self.layer_depth[layer + 1],
                                    scalar(0))
                    if float(mapmaximum(exceed_mm)) > float(1e-03):
                        self.report(exceed_mm, 'outEXz' + str(layer + 1))

                light_leached.append(getLeachedMass(self, layer, percolation[layer], self.lightmass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.LCH))
                heavy_leached.append(getLeachedMass(self, layer, percolation[layer], self.heavymass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.LCH))

                SW3 = self.theta[layer] * self.layer_depth[layer] - percolation[layer]
                self.theta[layer] = SW3 / self.layer_depth[layer]

                self.lightmass[layer] -= light_leached[layer]
                self.heavymass[layer] -= heavy_leached[layer]
                if mapminimum(self.lightmass[layer]) < 0:
                    print("Err LCHz1, light")
                if mapminimum(self.heavymass[layer]) < 0:
                    print("Err LCHz1, heavy")

            if mapminimum(self.theta[layer]) < 0:
                print('Error at Percolation, Layer ' + str(layer))

            if self.TEST_LCH:
                recordLCH(self, light_leached[layer], layer)
                self.report(self.lightmass[layer], 'LCH_Mz' + str(layer))

            if self.TEST_IR:
                if layer == 0:
                    recordInfiltration(self, infil_z0, layer)
                    recordRunOff(self, runoff_z0)
                else:
                    recordInfiltration(self, percolation[layer - 1], layer)

            if self.TEST_PERC:
                recordPercolation(self, percolation[layer], layer)

        # Artificial drainage (relevant layer)
        drained_layers = [n for n, x in enumerate(self.drainage_layers) if x is True]  # <- list of indexes
        adr_layer = int(drained_layers[0])  # <- 13.05.2018, implements only one layer (i.e. z2)!
        cell_drainge_outflow = getArtificialDrainage(self, adr_layer)  # mm
        light_drained = getDrainMassFlux(self, adr_layer, self.lightmass[adr_layer])  #
        heavy_drained = getDrainMassFlux(self, adr_layer, self.heavymass[adr_layer])  #
        self.lightmass[adr_layer] -= light_drained
        self.heavymass[adr_layer] -= heavy_drained

        SW4 = self.theta[adr_layer] * self.layer_depth[adr_layer] - cell_drainge_outflow
        self.theta[adr_layer] = SW4 / self.layer_depth[adr_layer]

        if mapminimum(self.theta[adr_layer]) < 0:
            print('Error at ADR, Layer ' + str(adr_layer))

        if mapminimum(self.lightmass[adr_layer]) < 0:
            print("Err DR, light")
            self.report(self.lightmass[adr_layer], 'ADR_Mz')

        if mapminimum(self.heavymass[adr_layer]) < 0:
            print("Err DR, heavy")

        # Artificial drainage (Outlet discharge)
        cell_drain_z2_m3 = cell_drainge_outflow * cellarea() / 1000  # m3
        accu_drain_m3 = accuflux(self.ldd_subs, cell_drain_z2_m3)  # m3
        out_drain_m3 = areatotal(accu_drain_m3, self.outlet_multi)

        # Evapotranspiration
        etp = []
        evap = []
        transp = []
        etp_m3 = deepcopy(self.zero_map)
        evap_m3 = deepcopy(self.zero_map)
        transp_m3 = deepcopy(self.zero_map)
        depth_evap = self.layer_depth[0] + self.layer_depth[1]
        for layer in range(self.num_layers):
            act_evaporation_layer = deepcopy(self.zero_map)
            if layer < 2:
                pot_evapor_layer = pot_evapor * self.layer_depth[layer] / depth_evap
                act_evaporation_layer = getActualEvap(self, layer, theta_wp, pot_evapor_layer, run=self.ETP)
            # Evaporation
            SW5 = self.theta[layer] * self.layer_depth[layer] - act_evaporation_layer
            self.theta[layer] = SW5 / self.layer_depth[layer]
            evap.append(act_evaporation_layer)
            evap_m3 += evap[layer] * cellarea() / 1000  # m3

            # Transpiration
            act_transpir_layer = getActualTransp(self, layer, root_depth_tot, root_depth[layer],
                                                 pot_transpir, theta_wp, depletable_water, run=self.ETP)
            act_transpir_layer *= self.f_transp
            SW6 = self.theta[layer] * self.layer_depth[layer] - act_transpir_layer
            self.theta[layer] = SW6 / self.layer_depth[layer]
            transp.append(act_transpir_layer)
            transp_m3 += transp[layer] * cellarea() / 1000  # m3

            etp.append(act_transpir_layer + act_evaporation_layer)
            etp_m3 += etp[layer] * cellarea() / 1000  # m3

            if mapminimum(self.theta[layer]) < 0:
                print('Error at ETP, Layer ' + str(layer))

        evap_m3 = areatotal(evap_m3, self.is_catchment)
        transp_m3 = areatotal(transp_m3, self.is_catchment)
        self.resW_accEvap_m3_tss.sample(evap_m3)
        self.resW_accTransp_m3_tss.sample(transp_m3)

        # Lateral flow
        catch_n_latflow_m3 = deepcopy(self.zero_map)  # Net
        cell_lat_outflow_m3 = deepcopy(self.zero_map)  # Only out

        latflow_outlet_mm = []
        latflow_cell_mm = []
        ligth_latflow = []
        heavy_latflow = []
        for layer in range(self.num_layers):
            # Check if basement layer
            if self.linear_reservoir and (layer == self.num_layers - 1):
                latflow_outlet_mm.append(deepcopy(self.zero_map))
                latflow_cell_mm.append(deepcopy(self.zero_map))
            else:
                # Outlet flux only first (needed to determine its own inflow capacity)
                # depth = self.layer_depth[layer]
                # c = self.c_lf[layer]
                # Outlet-cell flux Part I, independent of downstream capacity
                # outlet_flux1_mm = max(c * (depth * self.theta[layer] - depth * self.theta_fc[layer]),
                #                       scalar(0))  # [mm]
                # SW = self.theta[layer] * depth - outlet_flux1_mm
                # self.theta[layer] = ifthenelse(self.is_outlet, SW / depth, self.theta[layer])

                # Get lateral flow upstream cells
                latflow_dict = getLateralFlow(self, layer, run=self.LF)
                latflow_cell_mm.append(latflow_dict['cell_outflow'])  # flux map
                self.theta[layer] = latflow_dict['new_moisture']  # state map

                # Outlet-cell flux Part II: add flux map to Part I.
                outlet_flux2_mm = ifthenelse(self.is_outlet, latflow_dict['cell_outflow'], deepcopy(self.zero_map))
                # latflow_outlet_mm.append(outlet_flux1_mm + outlet_flux2_mm)
                latflow_outlet_mm.append(outlet_flux2_mm)

                if mapmaximum(self.theta[layer]) > mapmaximum(self.theta_sat[layer]):
                    val = mapmaximum(self.theta[layer]) - mapmaximum(self.theta_sat[layer])
                    self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                                   self.theta[layer])
                    if float(val) > float(1e-06):
                        print("Error at getLateralFlow(), SAT exceeded, layer " + str(layer) + ' by ' + str(val))

            # Pest. mass lateral flux
            ligth_latflow_dict = getLatMassFlux(self, layer, self.lightmass[layer], latflow_cell_mm[layer],
                                                debug=self.TEST_LFM, run=self.LFM)
            heavy_latflow_dict = getLatMassFlux(self, layer, self.heavymass[layer], latflow_cell_mm[layer],
                                                debug=self.TEST_LFM, run=self.LFM)
            ligth_latflow.append(ligth_latflow_dict['mass_loss'])
            heavy_latflow.append(heavy_latflow_dict['mass_loss'])

            self.lightmass[layer] = ligth_latflow_dict['new_mass']
            self.heavymass[layer] = heavy_latflow_dict['new_mass']
            cell_lat_outflow_m3 += latflow_outlet_mm[layer] * cellarea() / 1000  # m3, only out!

            if mapminimum(self.theta[layer]) < 0:
                print('Error at LF, Layer ' + str(layer))
            if mapminimum(self.lightmass[layer]) < 0:
                print("Err LF, light")
                self.report(self.lightmass[layer], 'LFM_Mz' + str(layer))
            if mapminimum(self.heavymass[layer]) < 0:
                print("Err LF, heavy")

        # Lateral flow (Outlet discharge)
        outlet_latflow_m3 = areatotal(cell_lat_outflow_m3, self.outlet_multi)  # Only outlet cells

        # Baseflow
        if self.linear_reservoir:
            SWbsmt = self.theta[-1] * self.layer_depth[-1]
            baseflow_mm = SWbsmt / self.k_g  # [mm/d]
            SWbsmt -= baseflow_mm
            if mapminimum(SWbsmt) < 0:
                val = float(mapminimum(SWbsmt))
                if val < float(-1e-06):
                    print("Negative Basement Soil Water")
            self.theta[-1] = max(SWbsmt / self.layer_depth[-1], scalar(0))

            accu_baseflow_m3 = accuflux(self.ldd_subs, baseflow_mm * cellarea() / 1000)
            out_baseflow_m3 = areatotal(accu_baseflow_m3, self.outlet_multi)
        else:
            out_baseflow_m3 = None

        light_deg = []
        ch_storage_light = []
        for layer in range(self.num_layers):
            # Temperature
            temp_dict = getLayerTemp(self, layer, bio_cover, temp_bare_soil)
            self.temp_surf_fin = temp_dict["temp_surface"]
            self.temp_fin[layer] = temp_dict["temp_layer"]

            # Degradation
            deg_light_dict = getMassDegradation(self, layer, theta_wp, self.lightmass[layer],
                                                frac="L", sor_deg_factor=1,
                                                debug=self.TEST_DEG, run=self.DEG)
            deg_heavy_dict = getMassDegradation(self, layer, theta_wp, self.heavymass[layer],
                                                frac="H", sor_deg_factor=1,
                                                debug=self.TEST_DEG, run=self.DEG)

            self.lightmass[layer] = deg_light_dict["mass_tot_new"]
            light_deg.append(deg_light_dict.get("mass_deg_aq") +
                             deg_heavy_dict.get("mass_deg_ads"))
            self.heavymass[layer] = deg_heavy_dict["mass_tot_new"]

            if mapminimum(self.lightmass[layer]) < 0:
                print("Err DEG, light")
                if layer == 0:
                    self.report(self.lightmass[layer], 'DEG_Mz' + str(layer))

            if mapminimum(self.heavymass[0]) < 0:
                print("Err DEG, heavy")

            # Change in mass storage after degradation - Pesticide Mass
            ch_storage_light.append(self.lightmass[layer] -
                                    self.lightmass_ini[layer])
            self.lightmass_ini[layer] = deepcopy(self.lightmass[layer])

            self.delta[layer] = ((self.heavymass[layer] / self.lightmass[layer] - self.r_standard) /
                                 self.r_standard) * 1000  # [permille]

        """ Layer analysis """
        for layer in range(self.num_layers):
            if layer == (self.num_layers - 1):
                getLayerAnalysis(self, layer, percolation, latflow_outlet_mm, evap, transp,
                                 root_depth, out_baseflow_m3=out_baseflow_m3)
            else:
                getLayerAnalysis(self, layer, percolation, latflow_outlet_mm, evap, transp,
                                 root_depth)

        # Update state variables
        # Change in storage - Moisture
        ch_storage = []
        ch_storage_m3 = deepcopy(self.zero_map)
        for layer in range(self.num_layers):
            ch_storage.append((self.theta[layer] * self.layer_depth[layer] * cellarea() / 1000) -
                              (self.theta_ini[layer] * self.layer_depth[layer] * cellarea() / 1000))
            self.theta_ini[layer] = deepcopy(self.theta[layer])
            ch_storage_m3 += ch_storage[layer]  # Reservoir storage (m3)

        getCatchmentStorage(self)
        getAverageMoisture(self)

        # Get Transect concentrations
        # Observed conc. is betw 2 and 8 ug/g dry soil (on transect)
        # Observed conc. can reach 20 ug/g dry soil (on single plot on application)

        # Record soil concentrations and isotopes
        if self.PEST:
            cell_mass = self.lightmass[0] + self.heavymass[0]
            cell_massXdelta = cell_mass * self.delta[0]

            north_sampling_pts = 30
            soils_north = reportNorthSoils(self, cell_mass, cell_massXdelta, north_sampling_pts)
            valley_sampling_pts = 25
            soils_valley = reportValleySoils(self, cell_mass, cell_massXdelta, valley_sampling_pts)
            south_sampling_pts = 26
            soils_south = reportSouthSoils(self, cell_mass, cell_massXdelta, south_sampling_pts)

        #####################
        # End of Model Loop #
        self.jd_cum += self.jd_dt  # updating JDcum, currently dt = 2 day

        ###################
        # Water Balance  ##
        ###################
        q_obs = timeinputscalar('q_obs_m3day.tss', nominal("outlet_v3"))

        # Total discharge
        tot_vol_disch_m3 = getTotalDischarge(out_runoff_m3,
                                             outlet_latflow_m3,
                                             out_drain_m3, baseflow=out_baseflow_m3)
        self.tot_Q_m3_tss.sample(tot_vol_disch_m3)

        computeGlobalWaterBalance(self, ch_storage_m3, percolation, etp_m3,
                                  catch_n_latflow_m3,
                                  cell_lat_outflow_m3,
                                  outlet_latflow_m3,
                                  tot_rain_m3, out_runoff_m3, q_obs, out_drain_m3,
                                  out_baseflow_m3=out_baseflow_m3)

        # Analysis (NASH Discharge)
        reportNashHydro(self, q_obs, tot_vol_disch_m3)

        if self.TEST_theta and self.currentTimeStep() % 2 == 0:
            checkMoisture(self, self.theta, 'athz')

        ######################
        # Pesticide Balance ##
        ######################
        # Applied ug on catchment
        catch_app_light = areatotal(light_applied, self.is_catchment)  #
        self.resM_accAPP_g_tss.sample(catch_app_light)

        if self.TRANSPORT:
            # Degradation
            light_deg_tot = deepcopy(self.zero_map)
            for layer in range(self.num_layers):
                light_deg_tot += light_deg[layer]
            catch_deg_light = areatotal(light_deg_tot, self.is_catchment)
            z0_light_deg_catch = areatotal(light_deg[0], self.is_catchment)
            self.cum_degZ0_L_g += z0_light_deg_catch
            self.resM_accDEG_L_tss.sample(catch_deg_light)
            self.resM_accDEGz0_L_tss.sample(z0_light_deg_catch)
            self.cum_degZ0_L_g_tss.sample(self.cum_degZ0_L_g)

            # Volatilized
            catch_volat_light = areatotal(light_volat, self.is_catchment)
            self.resM_accVOLAT_L_tss.sample(catch_volat_light)

            # Mass loss to run-off
            # Index: 0 <- light, Index: 2 <- heavy
            catch_runoff_light = areatotal(mass_runoff[0], self.is_catchment)
            catch_runoff_heavy = areatotal(mass_runoff[1], self.is_catchment)
            self.resM_accRO_L_tss.sample(catch_runoff_light)
            self.cum_roZ0_L_g += catch_runoff_light
            self.cum_roZ0_L_g_tss.sample(self.cum_roZ0_L_g)

            # z0-mass leached
            catch_leach_light_z0 = areatotal(light_leached[0], self.is_catchment)
            self.cum_lchZ0_L_g += catch_leach_light_z0
            self.resM_accLCHz0_L_tss.sample(catch_leach_light_z0)
            self.resM_cumLCHz0_L_g_tss.sample(self.cum_lchZ0_L_g)

            # z1-mass leached
            catch_leach_light_z1 = areatotal(light_leached[1], self.is_catchment)
            self.resM_accDPz1_L_tss.sample(catch_leach_light_z1)

            # Basement-mass leached = zero, if no basement percolation
            catch_leach_light_Bsmt = areatotal(light_leached[-1], self.is_catchment)
            self.resM_accDP_L_tss.sample(catch_leach_light_Bsmt)

            # Artificial drained mass (layer z2)
            catch_drain_light = areatotal(light_drained, self.is_catchment)
            catch_drain_heavy = areatotal(heavy_drained, self.is_catchment)
            self.resM_accADR_L_tss.sample(catch_drain_light)
            # catch_drain_heavy = areatotal(z1_heavy_drain, self.is_catchment)
            self.cum_adr_L_g += catch_drain_light
            self.cum_adr_L_g_tss.sample(self.cum_adr_L_g)

            # Lateral flux at outlet cells
            # outlet_cell_lightflux = (ligth_latflow_outlet[0] + ligth_latflow_outlet[2] +
            #                          ligth_latflow_outlet[1] + ligth_latflow_outlet[3])
            # outlet_cell_lightflux = areatotal(outlet_cell_lightflux, self.outlet_multi)  # Sum only outlet cells
            # self.cum_latflux_L_g += outlet_cell_lightflux
            # self.cum_latflux_L_g_tss.sample(self.cum_latflux_L_g)

            # outlet_cell_heavyflux = (heavy_latflow_outlet[0] + heavy_latflow_outlet[2] +
            #                          heavy_latflow_outlet[1] + heavy_latflow_outlet[3])
            # outlet_cell_heavyflux = areatotal(outlet_cell_heavyflux, self.outlet_multi)  # Sum only outlet cells

            # Lateral flux (Required in MB)
            latflux_light_catch = deepcopy(self.zero_map)
            latflux_heavy_catch = deepcopy(self.zero_map)
            for layer in range(self.num_layers):
                latflux_light_catch += ligth_latflow[layer]
                latflux_heavy_catch += heavy_latflow[layer]

            catch_latflux_light = areatotal(latflux_light_catch, self.outlet_multi)  # Needed for MB
            catch_latflux_heavy = areatotal(latflux_heavy_catch, self.outlet_multi)  # Needed for MB
            self.resM_accLF_L_tss.sample(catch_latflux_light)  # Reports the outlet-only loss

            # Baseflow flux
            # out_baseflow_light = areatotal(baseflow_light, self.is_catchment)
            # self.resM_accBF_L_tss.sample(out_baseflow_light)

            # Change in mass storage
            ch_storage_light_catch = deepcopy(self.zero_map)
            for layer in range(self.num_layers):
                ch_storage_light_catch += ch_storage_light[layer]

            catch_ch_storage_light = areatotal(ch_storage_light_catch, self.is_catchment)
            self.resM_accCHS_L_tss.sample(catch_ch_storage_light)

            ####################
            # Outlet Pesticide #
            ####################
            # Total mass export
            outlet_light_export = (catch_runoff_light + catch_drain_light + catch_latflux_light)
            outlet_heavy_export = (catch_runoff_heavy + catch_drain_heavy + catch_latflux_heavy)

            conc_ugL = (outlet_light_export + outlet_heavy_export) * 1e6 / (tot_vol_disch_m3 * 1e3)
            conc_ROFF_ug_L = (catch_runoff_light + catch_runoff_heavy) * 1e6 / (tot_vol_disch_m3 * 1e3)
            conc_LF_ug_L = (catch_latflux_light + catch_latflux_heavy) * 1e6 / (tot_vol_disch_m3 * 1e3)
            conc_ADR_ug_L = (catch_drain_light + catch_drain_heavy) * 1e6 / (tot_vol_disch_m3 * 1e3)
            #
            self.cum_exp_L_g += outlet_light_export
            self.resM_EXP_Smet_g_tss.sample(outlet_light_export)  # grams
            self.resM_oCONC_ugL_tss.sample(conc_ugL)  # ug/L
            self.resM_oCONC_ROFF_ugL_tss.sample(conc_ROFF_ug_L)  # ug/L
            self.resM_oCONC_LF_ugL_tss.sample(conc_LF_ug_L)  # ug/L
            self.resM_oCONC_ADR_ugL_tss.sample(conc_ADR_ug_L)  # ug/L
            self.resM_cumEXP_Smet_g_tss.sample(self.cum_exp_L_g)
            #
            # out_delta = ((outlet_heavy_export / outlet_light_export - self.r_standard) /
            #              self.r_standard) * 1000  # [permille]
            # roff_delta = ((catch_runoff_heavy / catch_runoff_light - self.r_standard) /
            #               self.r_standard) * 1000  # [permille]
            # latflux_delta = ((outlet_cell_heavyflux / outlet_cell_lightflux - self.r_standard) /
            #                  self.r_standard) * 1000  # [permille]
            # drain_delta = ((catch_drain_heavy / catch_drain_light - self.r_standard) /
            #                self.r_standard) * 1000  # [permille]

            # self.resM_outISO_d13C_tss.sample(out_delta)
            # self.resM_outISO_ROFF_d13C_tss.sample(roff_delta)
            # self.resM_outISO_LF_d13C_tss.sample(latflux_delta)
            # self.resM_outISO_ADR_d13C_tss.sample(drain_delta)

        reportGlobalPestBalance(self,
                                catch_app_light,
                                catch_deg_light,
                                catch_volat_light,
                                catch_runoff_light,
                                catch_leach_light_Bsmt,
                                catch_drain_light,
                                catch_latflux_light,
                                catch_ch_storage_light)

        # Analysis Pest
        if self.PEST:
            # TODO:
            # repNashOutConc()
            # repNashOutIso()
            # repNashOutCombined()

            # TODO:
            # Combine Composites and detailed
            reportNashConcComposites(self,
                                     soils_north['ave_conc'],
                                     soils_valley['ave_conc'],
                                     soils_south['ave_conc'])

            # reportNashDeltaComposites()

            # repNashSoilConc()
            # repNashSoilIso()

    def postmcloop(self):
        pass
        # names = ["q"]  # Discharge, Nash_Discharge
        # mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        # aguila --timesteps=[170,280,2] q-ave q-var outlet_v3.map
        # percentiles = [0.25, 0.5, 0.75]
        # mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())
        # aguila --quantiles=[0.25,0.75,0.25] --timesteps=[170,280,2] q


# Visualization
# aguila 2\at0dC000.177 2\at1dC000.177
# aguila --scenarios='{2,1}' --multi=1x4  --timesteps=[175,179,2] aLEACH aLEACHz aLF aLFz
# aguila --scenarios='{2}'  --timesteps=[100,280,2] az0dC az1dC az2dC
# aguila --scenarios='{2}'  --timesteps=[2,280,2] aHeight aRDtot aCrop aPotETP akcb akcb1 akcmax
#  aguila --scenarios='{2}'  --timesteps=[2,360,2] aHeight aRDtot aCrop akcb aPotTRA aPotEVA
#  aguila --scenarios='{2,1,3}' --multi=1x4 --timesteps=[1,300,1] athz0 athz1 athz2 athz3
#  aguila --scenarios='{2,1,3}' --multi=1x4  --timesteps=[2,300,2] athz0 athz1 athz2 athz3
#  aguila --scenarios='{2,1,3}' --multi=1x4 --timesteps=[2,300,2] aObj1 aObj2
#  aguila --scenarios='{2,1,3}' thFCz2 thSATz2
#  aguila --scenarios='{2}' --timesteps=[1,300,1] aROm3 athz0 athz1 athz2 athz3
#  aguila --scenarios='{2}' --timesteps=[1,300,1] aROm3 aZ1LCH
# aguila --scenarios='{2}' --timesteps=[2,300,2] z3EVA z3TRA z3ROOT

# Time series
# aguila 2\res_nash_q_m3.tss 6\res_nash_q_m3.tss
# aguila 2\resW_accStorage_m3.tss
# aguila 2\resM_norCONC.tss 2\resM_valCONC.tss 2\resM_souCONC.tss

nrOfSamples = int(runs)  # Samples are each a MonteCarlo realization
firstTimeStep = start_jday()  # 166 -> 14/03/2016
nTimeSteps = 360  # 360
myAlteck16 = BeachModel("clone_nom.map")  # an instance of the model, which inherits from class: DynamicModel
dynamicModel = DynamicFramework(myAlteck16, lastTimeStep=nTimeSteps,
                                firstTimestep=firstTimeStep)  # an instance of the Dynamic Framework
mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)

t0 = datetime.now()
# dynamicModel.run()
mcModel.run()
t1 = datetime.now()

duration = t1 - t0
tot_min = duration.total_seconds() / 60
print("Total minutes: ", tot_min)
print("Minutes/monte carlo", tot_min / int(runs))
print("Minutes/Yr: ", (duration.total_seconds() / 60) / (nTimeSteps - firstTimeStep) * 365)

if morris:
    # First test of soil-water holding capacity on:
    ouput_vars = ["resNash_q_m3",
                  "resW_q_sim_ave_m3",
                  "resW_cum_q_sim_m3",
                  "resM_cumDEGz0_L",  # Mass degraded in z0
                  "resM_cumLCHz0_L",  # Mass leached out of z0 (loss to subsurface)
                  "resM_cumEXP_L",  # Accounts for mass runoff, drainage, outlet-cell-flux and baseflow
                  "resNash_compConc_L"
                  # "resN_nash_compDelta"
                  ]
    morris_results = getSiList(ouput_vars, problem, param_values, grid_jump, p, runs)
    saveMorris(ouput_vars, morris_results)
