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

print(os.getcwd())

global state
state = -1

morris = False
if morris:
    from morris_test import *
else:
    runs = 3

"""
Testing decreases in gamma2 and gamma3,
with SCN accounting ofr 3 layers depth
"""

# 1st
# - Generate the set of input parameters (see: morris_analysis.py)
# - Make input parameters a global numpy array (for model access)


# def getTimeStamp(timestep, sep=','):
#     path = "Data/Time.csv"
#     obs = read_csv(path, sep=sep)
#     obs = obs[['Jdays', 'Date']]
#     obs_dict = obs.to_dict(orient='split')
#     obs_dict['data']
#     return str(obs_dict['data'][timestep-1][1])


def get_state(old_state):
    new_state = old_state + 1
    global state
    state = new_state
    return state


def start_jday():
    start_sim = 166
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
        self.TEST = True
        self.TEST_depth = False
        self.TEST_roots = False
        self.TEST_Ksat = False
        self.TEST_thProp = False
        self.TEST_theta = True
        self.TEST_IR = True
        self.TEST_PERC = False

        self.PEST = False
        self.TRANSPORT = True
        self.DEGRADE = True
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
        self.ldd_subs = lddcreate(out_burn, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI

        self.zero_map = out_burn - out_burn  # Zero map to generate scalar maps
        self.mask = out_burn / out_burn
        self.aging = deepcopy(self.zero_map)  # Cumulative days after application on each pixel

        # self.outlet = self.readmap("outlet_multi")
        self.outlet_multi = self.readmap("out_multi_nom")  # Multi-outlet with 0 or 1
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
        self.tot_rain_m3_tss = TimeoutputTimeseries("resW_accRain_m3", self, nominal("outlet_true"), noHeader=False)
        # Runoff
        self.out_runoff_m3_tss = TimeoutputTimeseries("resW_accRunoff_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_runoff_m3_tss = TimeoutputTimeseries("resW_totRunoff_m3", self, nominal("outlet_true"), noHeader=False)
        # Percolation
        self.out_percol_bsmt_m3_tss = TimeoutputTimeseries("resW_accPercol_z2_m3", self, nominal("outlet_true"),
                                                           noHeader=False)
        self.tot_perc_z3_m3_tss = TimeoutputTimeseries("resW_totPercol_z2_m3", self, nominal("outlet_true"),
                                                       noHeader=False)
        # ETP
        self.out_etp_m3_tss = TimeoutputTimeseries("resW_accEtp_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_etp_m3_tss = TimeoutputTimeseries("resW_totEtp_m3", self, nominal("outlet_true"), noHeader=False)

        self.out_baseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", self, nominal("outlet_true"),
                                                        noHeader=False)
        self.tot_baseflow_m3_tss = TimeoutputTimeseries("resW_totBaseflow_m3", self, nominal("outlet_true"),
                                                        noHeader=False)
        # LF Drainage
        self.out_accu_drain_m3_tss = TimeoutputTimeseries("resW_o_accDrain_m3", self, nominal("outlet_true"),
                                                          noHeader=False)
        self.tot_accu_drain_m3_tss = TimeoutputTimeseries("resW_o_totDrain_m3", self, nominal("outlet_true"),
                                                          noHeader=False)  # Cumulative ADR
        # LF options
        self.sat_accu_overflow_m3_tss = TimeoutputTimeseries("resW_of_accLatflow_m3", self, nominal("outlet_true"),
                                                             noHeader=False)

        self.tot_accu_of_latflow_m3_tss = TimeoutputTimeseries("resW_of_totLatflow_m3", self, nominal("outlet_true"),
                                                               noHeader=False)
        # Inflow
        self.out_cell_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_cellLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.out_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_accLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        # Outflow
        self.out_cell_o_latflow_m3_tss = TimeoutputTimeseries("resW_o_cellLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.out_accu_o_latflow_m3_tss = TimeoutputTimeseries("resW_o_accLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_o_latflow_m3_tss = TimeoutputTimeseries("resW_o_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.out_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_accLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_n_latflow_m3_tss = TimeoutputTimeseries("resW_n_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)

        self.out_ch_storage_m3_tss = TimeoutputTimeseries("resW_accChStorage_m3", self, nominal("outlet_true"),
                                                          noHeader=False)
        self.global_mb_water_tss = TimeoutputTimeseries("resW_global_waterMB", self, nominal("outlet_true"),
                                                        noHeader=False)
        self.storage_m3_tss = TimeoutputTimeseries("resW_accStorage_m3", self, nominal("outlet_true"), noHeader=False)

        # PESTI
        # Pesticide
        self.global_mb_pest_tss = TimeoutputTimeseries("resM_global_mb_pest", self, nominal("outlet_true"),
                                                       noHeader=False)
        self.out_app_L_tss = TimeoutputTimeseries("resM_accAPP_L", self, nominal("outlet_true"), noHeader=False)
        self.out_volat_L_tss = TimeoutputTimeseries("resM_accVOL_L", self, nominal("outlet_true"), noHeader=False)
        self.out_runoff_L_tss = TimeoutputTimeseries("resM_accRO_L", self, nominal("outlet_true"), noHeader=False)
        self.out_degZ0_L_tss = TimeoutputTimeseries("resM_accDEG_L", self, nominal("outlet_true"), noHeader=False)
        self.out_leachZ0_L_tss = TimeoutputTimeseries("resM_accLCH_L", self, nominal("outlet_true"), noHeader=False)
        self.out_leachZ1_L_tss = TimeoutputTimeseries("resM_accDPz1_L", self, nominal("outlet_true"), noHeader=False)

        self.out_leach_L_tss = TimeoutputTimeseries("resM_accDP_L", self, nominal("outlet_true"), noHeader=False)
        self.out_drain_L_tss = TimeoutputTimeseries("resM_accADR_L", self, nominal("outlet_true"), noHeader=False)
        self.out_latflux_L_tss = TimeoutputTimeseries("resM_accLF_L", self, nominal("outlet_true"), noHeader=False)
        self.out_baseflow_L_tss = TimeoutputTimeseries("resM_accBF_L", self, nominal("outlet_true"), noHeader=False)
        self.out_chstorage_L_tss = TimeoutputTimeseries("resM_accCHS_L", self, nominal("outlet_true"), noHeader=False)

        # Cumulative Pesticide
        self.cum_degZ0_L_g_tss = TimeoutputTimeseries("resM_cumDEGz0_L", self, nominal("outlet_true"),
                                                      noHeader=False)  # Deg z0
        self.cum_lchZ0_L_g_tss = TimeoutputTimeseries("resM_cumLCHz0_L", self, nominal("outlet_true"),
                                                      noHeader=False)  # Leaching z0
        self.cum_roZ0_L_g_tss = TimeoutputTimeseries("resM_cumROz0_L", self, nominal("outlet_true"),
                                                     noHeader=False)  # Runoff
        self.cum_adr_L_g_tss = TimeoutputTimeseries("resM_cumADR_L", self, nominal("outlet_true"),
                                                    noHeader=False)  # Art. drainage
        self.cum_latflux_L_g_tss = TimeoutputTimeseries("resM_cumLF_L", self, nominal("outlet_true"),
                                                        noHeader=False)  # Soil column, outlet cells
        self.cum_exp_L_g_tss = TimeoutputTimeseries("resM_cumEXP_L", self, nominal("outlet_true"),
                                                    noHeader=False)  # Total cum. exports
        self.nash_compConc_L_tss = TimeoutputTimeseries("resNash_compConc_L", self, nominal("outlet_true"),
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
        self.i_Q_m3_tss = TimeoutputTimeseries("resW_i_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.o_Q_m3_tss = TimeoutputTimeseries("resW_o_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_Q_m3_tss = TimeoutputTimeseries("resW_tot_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.q_obs_cum_tss = TimeoutputTimeseries("resW_cum_q_obs_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # Equivalent to net_Q
        self.rain_obs_cum_tss = TimeoutputTimeseries("resW_cum_rain_obs_m3", self, nominal("outlet_true"),
                                                     noHeader=False)  # Equivalent to net_Q
        self.rest_obs_tss = TimeoutputTimeseries("resW_q_restit_obs_m3", self, nominal("outlet_true"),
                                                 noHeader=False)  # = rain/q_obs
        self.q_sim_cum_tss = TimeoutputTimeseries("resW_cum_q_sim_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # Sum sim discharge (if obs available).
        self.q_sim_ave_tss = TimeoutputTimeseries("resW_q_sim_ave_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # This is 'Nash_q' as time series.
        self.nash_q_tss = TimeoutputTimeseries("resNash_q_m3", self, nominal("outlet_true"),
                                               noHeader=False)  # This is 'Nash_q' as time series.
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
        self.bsmntIsPermeable = True  # basement percolation (DP)
        self.ADLF = True

        # Morris tests
        m_state = get_state(state)  # First run will return state = 0

        """ Physical parameters for each layer """
        self.gamma = []  # coefficient to calibrate Ksat1
        self.s = []
        self.c_lf = []
        for layer in range(self.num_layers):
            self.gamma.append(scalar(self.ini_param.get("gamma" + str(layer))))  # percolation coefficient
            self.s.append(scalar(self.ini_param.get("s" + str(layer))))  # calibrate Ksat, mm/day
            self.c_lf.append(scalar(self.ini_param.get("c" + str(layer))))  # subsurface flow coefficient)

        self.fc_adj = scalar(self.ini_param.get("fc_adj"))  # Adjusting Paul's FC (equivalent so far)
        self.root_adj = scalar(self.ini_param.get("root_adj"))  # Adjust Root Depth factor
        self.c_adr = scalar(self.ini_param.get("c_adr"))
        self.drainage_layers = [False, False, True, False]  # z2 turned on!!
        z3_factor = scalar(self.ini_param.get("z3_factor"))  # Represents top proportion of bottom-most layer
        if m_state == 0:
            pass
        elif m_state == 1:
            # for layer in range(self.num_layers):
            #     self.c_lf[layer] += 0.1
            # z3_factor = scalar(0.7)
            self.c_adr /= 2
            # self.root_adj *= 0.7
            # epsilon = -1 * 1.369  # high deg
            # self.gamma[3] *= 0.5
        elif m_state == 2:
            # for layer in range(self.num_layers):
            #     self.c_lf[layer] += 0.2
            # z3_factor = scalar(0.4)
            self.c_adr /= 3
            # self.root_adj *= 0.5
            # epsilon = -1 * 1.476  # mid deg
            # self.gamma[2] *= 0.5
            # self.gamma[3] *= 0.25

        self.k = scalar(self.ini_param.get("k"))  # coefficient of declining LAI in end stage

        """ Soil Properties """
        self.p_b = scalar(self.ini_param.get("p_b"))  # Soil bulk density (g/cm^3)
        self.f_oc = scalar(self.ini_param.get("f_oc"))  # Organic Carbon in soil without grass (kg/kg)

        """
        Sorption parameters
        """
        # K_oc - S-metolachlor (K_oc in ml/g)
        # Marie's thesis: log(K_oc)= 1.8-2.6 [-] -> k_oc = 63 - 398 (Alletto et al., 2013).
        # Pesticide Properties Database: k_oc = 120 ml/g (range: 50-540 mL/g)
        self.k_oc = scalar(self.ini_param.get("k_oc"))  # ml/g
        self.k_d = self.k_oc * self.f_oc  # Dissociation coefficient K_d (mL/g = L/Kg)

        # Pesticide Properties Database states :
        # K_d=0.67; but here
        # K_d=120*0.021 = 2.52;  (L/kg)
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
        # self.alpha_iso = scalar(self.ini_param.get("alpha_iso"))  # 1 is no fractionation

        # Degradation Scenarios
        epsilon = -1 * 1.743  # low deg
        self.alpha_iso = epsilon / 1000 + 1

        """
        Layer depths
        """
        # Lowest layer depth and baseflow
        self.k_g = 1500  # [days] = 4 yrs
        self.gw_factor = 1 - z3_factor  # Represents bottom-most portion of bottom layer

        self.layer_depth = []
        self.tot_depth = deepcopy(self.zero_map)
        for layer in range(self.num_layers):
            if layer < self.num_layers - 1:
                self.layer_depth.append(self.zero_map +
                                        scalar(self.ini_param.get('z' + str(layer))))
                self.tot_depth += self.layer_depth[layer]
            else:
                self.layer_depth.append((self.datum_depth +  # total height
                                         scalar(self.ini_param.get('z' + str(layer)))  # plus a min-depth
                                         - self.tot_depth) * z3_factor)  # minus:(z0, z1, z2)*decreasing depth factor
                self.tot_depth += self.layer_depth[layer]
            if self.TEST_depth:
                checkLayerDepths(self, layer)

        self.smp_depth = self.layer_depth[0]

        """
        Hydro Maps
        """
        # Initial moisture (Final from model v1, Sept 30, 2016)
        self.theta = []
        if start_jday() < 166:
            self.theta.append(readmap('d14_theta_z0'))  # map of initial soil moisture in top layer (-)
            self.theta.append(readmap('d14_theta_z0'))  # d14_theta_z1
            self.theta.append(readmap('d14_theta_z2'))
            self.theta.append(readmap("d14_theta_z3"))  # * z3_factor + scalar(0.6) * self.gw_factor
        else:
            self.theta.append(readmap('d166_theta_z0'))  # map of initial soil moisture in top layer (-)
            self.theta.append(readmap('d166_theta_z1'))
            self.theta.append(readmap('d166_theta_z2'))
            self.theta.append(readmap("d166_theta_z3"))  # * z3_factor + scalar(0.6) * self.gw_factor

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
        # in ug = conc. (ug/g soil) * density (g/cm3) * (10^6 cm3/m3)*(1 m/10^3 mm)* depth_layer(mm) * cellarea(m2)
        # g = ug * 10e-06
        self.sm_background = []
        mean_back_conc = [0.06, 0.03, 0.001, 0.001]
        for layer in range(self.num_layers):
            background = ((self.zero_map + mean_back_conc[layer]) * self.p_b * scalar(10 ** 6 / 10 ** 3) *
                          self.layer_depth[layer] * cellarea() * 10e-06)  # Based on detailed soils
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

        # Assign dosages based on Farmer-Crop combinations [g/m2]
        fa_cr = readmap("crop_burn")  # Contains codes to assign appropriate dosage
        self.apps = getApplications(self, fa_cr)  # returns list of applied masses

        # Applications delta
        # Use map algebra to produce a initial signature map,
        # ATT: Need to do mass balance on addition of new layer.
        # where app1 > 0, else background sig. (plots with no new mass will be 0)
        # where app1 > 0, else background sig. (plots with no new mass will be 0)
        self.appDelta = []
        self.appDelta.append( ifthenelse(self.apps[0] > 0, scalar(-32.3), scalar(-23.7)) )
        self.appDelta.append( ifthenelse(self.apps[1] > 0, scalar(-32.3), scalar(-23.7)) )
        self.appDelta.append( ifthenelse(self.apps[2] > 0, scalar(-32.3), scalar(-23.7)) )

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
        self.lag = scalar(0.8)  # lag coefficient (-), 0 < lag < 1; -> in SWAT, lag = 0.80
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
        self.temp_ave_air = scalar(12.2)  # 12.2 is for Layon

        """
        Simulation start time: Oct 1st, 2015
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
        self.tot_nlf_m3 = 0
        self.tot_ilf_m3 = 0  # upstream inflow
        self.tot_olf_m3 = 0  # downstream outflow

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

        self.theta_sat = [deepcopy(self.zero_map), deepcopy(self.zero_map),
                          deepcopy(theta_sat_z2), deepcopy(theta_sat_z2)]
        self.theta_fc = [deepcopy(self.zero_map), deepcopy(self.zero_map),
                         deepcopy(theta_fcap_z2), deepcopy(theta_fcap_z2)]

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
        setglobaloption('matrixtable')  # allows lookupscalar to read more than 2 expressions.
        crop_type = lookupscalar('croptable.tbl', 1, fields)  # (table, col-value in crop table, column-value in fields)
        sow_yy = lookupscalar('croptable.tbl', 2, fields)
        sow_mm = lookupscalar('croptable.tbl', 3, fields)  # sowing or Greenup month
        sow_dd = lookupscalar('croptable.tbl', 4, fields)  # sowing day
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
        # Sugar beet = 0.7 - 1.2
        # Corn = 1.0 - 1.7
        # Grazing pasture 0.5 - 1.5
        # Spring Wheat = 1.0 -1.5
        # Winter Wheat = 1.5 -1.8
        # Apple trees = 1.0-2.0

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
        k_sat = [k_sat_z0z1, deepcopy(k_sat_z0z1),
                 k_sat_z2z3, deepcopy(k_sat_z2z3)]
        CN2_A = lookupscalar('croptable.tbl', 18, fields)  # curve number of moisture condition II
        CN2_B = lookupscalar('croptable.tbl', 19, fields)  # curve number of moisture condition II
        CN2_C = lookupscalar('croptable.tbl', 20, fields)  # curve number of moisture condition II

        if self.TEST_Ksat:
            reportKsatEvolution(self, k_sat)

        """
        Time-series data to spatial location,
        map is implicitly defined as the clonemap.
        """
        precip = timeinputscalar('rain.tss', 1)  # daily precipitation data as time series (mm)

        temp_bare_soil = timeinputscalar('T_bare.tss', nominal('clone_nom'))  # SWAT, Neitsch2009, p.43.
        self.temp_air = timeinputscalar('airTemp.tss', nominal('clone_nom'))
        et0 = timeinputscalar('ET0.tss', 1)  # daily ref. ETP at Zorn station (mm)
        wind = timeinputscalar('U2.tss', 1)  # wind speed time-series at 2 meters height
        humid = timeinputscalar('RHmin.tss', 1)  # minimum relative humidity time-series # PA: (-)
        # precipVol = precip * cellarea() / 1000  # m3

        ################
        # Crop growth ##
        ################
        jd_sow = convertJulian(sow_yy, sow_mm, sow_dd)
        self.rain_cum_mm += precip
        self.rain_cum_mm = ifthenelse(jd_sim == jd_sow, scalar(0), self.rain_cum_mm)
        # self.report(self.rain_cum_mm, 'aCuRain')
        CN2 = ifthenelse(self.rain_cum_mm > 90, CN2_C, CN2_A)
        # soil_group = ifthenelse(self.rain_cum_mm > 90, ordinal(3), ordinal(1))
        all_stages = len_grow_stage_ini + len_dev_stage + len_mid_stage + len_end_stage

        # updating of sowing date by land use
        # TODO: Why update if sow date is already by land use?
        # sow_yy = ifthenelse(jd_sim < jd_sow + all_stages, sow_yy,
        #                     ifthenelse(jd_sim < jd_sow + all_stages + 365, sow_yy + 1,
        #                                ifthenelse(jd_sim < jd_sow + all_stages + 730, sow_yy + 2,
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
        # frac_soil_cover = 1 - exp(-mu * LAI)
        # \mu is a light-use efficiency parameter that
        # depends on land-use characteristics
        # (i.e. Grass: 0.35; Crops: 0.45; Trees: 0.5-0.77; cite: Larcher, 1975).

        # TODO: Check "f" definition by Allan et al., 1998 against previous (above)
        # fraction of soil cover is calculated inside the "getPotET" function.
        # frac_soil_cover = ((Kcb - Kcmin)/(Kcmax - Kcmin))**(1+0.5*mean_height)
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
        # bcv should range 0 (bare soil) to 1 (complete cover)
        # self.report(bio_cover, 'aBCV')

        # Non-hydrological, non-reactive transport on z0
        #########################
        # Applications
        mass_applied = deepcopy(self.zero_map)
        light_applied = deepcopy(self.zero_map)
        heavy_applied = deepcopy(self.zero_map)
        light_volat = deepcopy(self.zero_map)

        # Volatilization (on application days only)
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

            if self.TRANSPORT:
                # Mass volatilized
                light_volat = getVolatileMass(self, self.temp_air, self.lightmass[0],  # "LL",
                                              rel_diff_model='option-1', sorption_model="linear",
                                              gas=True, run=self.PEST)
                heavy_volat = getVolatileMass(self, self.temp_air, self.heavymass[0],  # "HH",
                                              rel_diff_model='option-1', sorption_model="linear",
                                              gas=True, run=self.PEST)

                light_volat = ifthenelse(mass_applied > 0, light_volat, scalar(0))
                heavy_volat = ifthenelse(mass_applied > 0, heavy_volat, scalar(0))

                self.lightmass[0] -= light_volat
                self.heavymass[0] -= heavy_volat

        #########################
        infil_z0 = deepcopy(self.zero_map)
        runoff_z0 = deepcopy(self.zero_map)
        percolation = []
        mass_runoff = []  # 0 <- light, 1 <- heavy
        light_leached = []
        heavy_leached = []
        permeable = True
        # # Infiltration, runoff, & percolation (all layers)
        for layer in range(self.num_layers):
            if layer == 0:  # Layer 0
                z0_IRO = getTopLayerInfil(self, precip, theta_wp, CN2, crop_type,
                                          jd_sim, jd_dev, jd_mid, jd_end, len_dev_stage)
                # Partition infiltration
                infil_z0 = z0_IRO.get("infil_z0")  # [mm]
                infil_z1 = z0_IRO.get("infil_z1")  # [mm]
                runoff_z0 = z0_IRO.get("roff")  # [mm]

                # Distribute to layer 0
                self.theta[layer] += infil_z0 / self.layer_depth[layer]  # [-]
                # Distribution to layer 1 needed here, bc. getPercolation(layer = 0) will check below's capacity
                self.theta[layer + 1] += infil_z1 / self.layer_depth[layer + 1]  # [-]

                # RunOff Mass
                # Mass & delta run-off (RO)
                mass_runoff.append(getRunOffMass(self, precip, runoff_z0, self.lightmass[layer],
                                                 transfer_model="simple-mt", sorption_model="linear",
                                                 gas=True, run=self.PEST))
                mass_runoff.append(getRunOffMass(self, precip, runoff_z0, self.heavymass[layer],
                                                 transfer_model="simple-mt", sorption_model="linear",
                                                 gas=True, run=self.PEST))
                self.lightmass[layer] -= mass_runoff[0]  # light
                self.heavymass[layer] -= mass_runoff[1]  # heavy

                percolation.append(getPercolation(self, layer, k_sat[layer], isPermeable=permeable))  # [mm]
                water_flux = infil_z0 + infil_z1 + percolation[layer]

                light_leached.append(getLeachedMass(self, layer, water_flux, self.lightmass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.PEST))
                heavy_leached.append(getLeachedMass(self, layer, water_flux, self.heavymass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.PEST))
                self.theta[layer] -= percolation[layer] / self.layer_depth[layer]
                self.lightmass[layer] -= light_leached[layer]
                self.heavymass[layer] -= heavy_leached[layer]

            else:  # Layers 1, 2 & 3
                if layer == (self.num_layers - 1):
                    permeable = self.bsmntIsPermeable
                self.theta[layer] += percolation[layer - 1] / self.layer_depth[layer - 1]
                self.lightmass[layer] += light_leached[layer - 1]
                self.heavymass[layer] += heavy_leached[layer - 1]
                percolation.append(getPercolation(self, layer, k_sat[layer], isPermeable=permeable))
                light_leached.append(getLeachedMass(self, layer, percolation[layer], self.lightmass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.PEST))
                heavy_leached.append(getLeachedMass(self, layer, percolation[layer], self.heavymass[layer],
                                                    sorption_model="linear", leach_model="mcgrath", gas=True,
                                                    debug=self.DEBUG, run=self.PEST))
                self.theta[layer] -= percolation[layer] / self.layer_depth[layer]
                self.lightmass[layer] -= light_leached[layer]
                self.heavymass[layer] -= heavy_leached[layer]

            if self.TEST_IR:
                if layer == 0:
                    recordInfiltration(self, infil_z0, layer)
                    recordRunOff(self, runoff_z0)
                else:
                    recordInfiltration(self, percolation[layer-1], layer)

            if self.TEST_PERC:
                recordPercolation(self, percolation[layer], layer)


        # Artificial drainage (relevant layer)
        drained_layers = [n for n, x in enumerate(self.drainage_layers) if x is True]  # <- list of indexes
        adr_layer = int(drained_layers[0])  # <- 13.05.2018, implements only one layer (i.e. z2)!
        cell_drainge_outflow = getArtificialDrainage(self, drained_layers)  # mm
        light_drained = getDrainMassFlux(self, adr_layer, self.lightmass[adr_layer])  #
        heavy_drained = getDrainMassFlux(self, adr_layer, self.heavymass[adr_layer])  #
        self.lightmass[adr_layer] -= light_drained
        self.heavymass[adr_layer] -= heavy_drained
        self.theta[adr_layer] -= cell_drainge_outflow / self.layer_depth[adr_layer]

        # Evapotranspiration
        etp = []
        for layer in range(self.num_layers):
            act_transpir_layer = deepcopy(self.zero_map)
            act_evaporation_layer = deepcopy(self.zero_map)
            act_transpir_layer = getActualTransp(self, layer, root_depth_tot, root_depth[layer],
                                                 pot_transpir, theta_wp, depletable_water)
            if layer < 2:
                act_evaporation_layer = getActualEvap(self, layer, theta_wp, pot_evapor)
            etp.append(act_transpir_layer + act_evaporation_layer)
            self.theta[layer] -= etp[layer] / self.layer_depth[layer]

        # Lateral flow
        latflow_net = []  # every cell
        latflow_out = []
        latflow_in = []
        ligth_latflow = []
        ligth_latflow_outlet = []
        heavy_latflow = []
        for layer in range(self.num_layers):
            latflow_dict = getLateralFlow(self, layer)
            latflow_net.append(latflow_dict['lateral_flow_layer'])
            latflow_out.append(latflow_dict['cell_outflow'])
            latflow_in.append(latflow_dict['upstream_cell_inflow'])

            ligth_latflow_dict = getLatMassFlux(self, layer, self.lightmass[layer],
                                                latflow_out[layer], latflow_in[layer])
            heavy_latflow_dict = getLatMassFlux(self, layer, self.heavymass[layer],
                                                latflow_out[layer], latflow_in[layer])
            ligth_latflow.append(ligth_latflow_dict['net_mass_latflux'])
            ligth_latflow_outlet.append(ligth_latflow_dict['mass_loss'])

            self.theta[layer] += latflow_dict['lateral_flow_layer'] / self.layer_depth[layer]
            self.lightmass[layer] += ligth_latflow_dict['net_mass_latflux']
            self.heavymass[layer] += heavy_latflow_dict['net_mass_latflux']

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
                                                debug=self.DEBUG, run=self.PEST)
            deg_heavy_dict = getMassDegradation(self, layer, theta_wp, self.heavymass[layer],
                                                frac="H", sor_deg_factor=1,
                                                debug=self.DEBUG, run=self.PEST)

            self.lightmass[layer] = deg_light_dict["mass_tot_new"]
            light_deg.append(deg_light_dict.get("mass_deg_aq") +
                             deg_heavy_dict.get("mass_deg_ads"))
            self.heavymass[layer] = deg_heavy_dict["mass_tot_new"]

            # Change in mass storage after degradation - Pesticide Mass
            ch_storage_light.append(self.lightmass[layer] -
                                    self.lightmass_ini[layer])
            self.lightmass_ini[layer] = self.lightmass[layer]

            self.delta[layer] = ((self.heavymass[layer] / self.lightmass[layer] - self.r_standard) /
                                 self.r_standard) * 1000  # [permille]

        # Update state variables
        # Change in storage - Moisture
        ch_storage = []
        for layer in range(self.num_layers):
            ch_storage.append((self.theta[layer] * self.layer_depth[layer] * cellarea() / 1000) -
                              (self.theta_ini[layer] * self.layer_depth[layer] * cellarea() / 1000))
            self.theta_ini[layer] = deepcopy(self.theta[layer])

        # Get Transect concentrations
        # Observed conc. is betw 1 and 8 ug/g dry soil (on transect)
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

        ###########################################################################
        # FINAL MASS BALANCE #
        ######################
        # 'Sample' the time-series associated to each component (e.g. runoff)
        # at the outlet or due to accuflux()

        ######################
        # Water Balance
        ######################
        if self.TEST_theta and self.currentTimeStep() % 10 == 0:
            checkMoisture(self, self.theta, 'athz')

        # Precipitation total
        rain_m3 = self.mask * precip * cellarea() / 1000  # m3
        tot_rain_m3 = areatotal(rain_m3, self.is_catchment)
        self.tot_rain_m3_tss.sample(tot_rain_m3)

        # Discharge due to runoff at the outlet
        runoff_m3 = runoff_z0 * cellarea() / 1000  # m3
        accu_runoff_m3 = accuflux(self.ldd_subs, runoff_m3)
        out_runoff_m3 = areatotal(accu_runoff_m3, self.outlet_multi)
        self.out_runoff_m3_tss.sample(out_runoff_m3)  # save to outlet
        # self.obs_cum_runoff_m3_tss.sample(out_runoff_m3)  # save to sample locations

        # Discharge due to z1 artificial drainage
        # cell_drain_z1_m3 = drain_outflow_z1 * cellarea() / 1000  # m3
        # accu_o_drain_z1_m3 = accuflux(self.ldd_subs, cell_drain_z1_m3)  # m3
        # out_drain_z1_m3 = areatotal(accu_o_drain_z1_m3, self.outlet_multi)
        # self.out_accu_drain_m3_tss.sample(out_drain_z1_m3)

        # Discharge due to z2 artificial drainage
        cell_drain_z2_m3 = cell_drainge_outflow * cellarea() / 1000  # m3
        accu_drain_m3 = accuflux(self.ldd_subs, cell_drain_z2_m3)  # m3
        out_drain_m3 = areatotal(accu_drain_m3, self.outlet_multi)
        self.out_accu_drain_m3_tss.sample(out_drain_m3)

        # Discharge due to lateral flow (3 types: overflow, cell-outlet flow, catchment net)
        # Overflow
        # overflow_cellvol_z0 = overflow_height_z0 * cellarea() / 1000  # m3
        # overflow_cellvol_z1 = overflow_height_z1 * cellarea() / 1000  # m3
        # overflow_cellvol_z2 = overflow_height_z2 * cellarea() / 1000  # m3
        # column_overflow = overflow_cellvol_z2 + overflow_cellvol_z1 + overflow_cellvol_z0
        # accu_of_latflow_m3 = accuflux(self.ldd_subs, column_overflow)
        # of_latflow_m3 = areatotal(accu_of_latflow_m3, self.outlet_multi)
        # self.sat_accu_overflow_m3_tss.sample(of_latflow_m3)

        # In-flow (for each cell at the outlet)
        # in_latflow_z0_m3 = lat_inflow_z0 * cellarea() / 1000  # m3
        # in_latflow_z1_m3 = lat_inflow_z1 * cellarea() / 1000
        # in_latflow_z2_m3 = lat_inflow_z2 * cellarea() / 1000
        # in_latflow_z3_m3 = lat_inflow_z3 * cellarea() / 1000
        # cell_lat_inflow_m3 = in_latflow_z0_m3 + in_latflow_z1_m3 + in_latflow_z2_m3 + in_latflow_z3_m3
        # self.report(cell_lat_inflow_m3, 'cellIn')
        # outlet_cell_inflow_m3 = areatotal(cell_lat_inflow_m3, self.outlet_multi)  # Only outlet cells
        # self.out_cell_i_latflow_m3_tss.sample(outlet_cell_inflow_m3)
        # catch_lat_inflow_m3 = areatotal(cell_lat_inflow_m3, self.is_catchment)  # Sum catchment (needed for MB)
        # self.out_accu_i_latflow_m3_tss.sample(catch_lat_inflow_m3)

        # Out-flow (for each cell at the outlet)
        out_latflow_z0_m3 = latflow_out[0] * cellarea() / 1000  # m3
        out_latflow_z1_m3 = latflow_out[1] * cellarea() / 1000
        out_latflow_z2_m3 = latflow_out[2] * cellarea() / 1000
        out_latflow_z3_m3 = latflow_out[3] * cellarea() / 1000
        cell_lat_outflow_m3 = out_latflow_z0_m3 + out_latflow_z1_m3 + out_latflow_z2_m3 + out_latflow_z3_m3
        # self.report(cell_lat_outflow_m3, 'cellOut')
        outlet_cell_outflow_m3 = areatotal(cell_lat_outflow_m3, self.outlet_multi)  # Only outlet cells
        self.out_cell_o_latflow_m3_tss.sample(outlet_cell_outflow_m3)
        catch_lat_outflow_m3 = areatotal(cell_lat_outflow_m3, self.is_catchment)  # Sum catchment
        self.out_accu_o_latflow_m3_tss.sample(catch_lat_outflow_m3)

        # Netflow
        latflow_net_m3 = []
        catch_n_latflow_m3 = deepcopy(self.zero_map)
        for layer in range(self.num_layers):
            vol_by_cell = latflow_net[layer] * cellarea() / 1000  # m3
            latflow_net_m3.append(vol_by_cell)  # <- volume on each layer's cell
            catch_n_latflow_m3 += vol_by_cell  # <- all layers

        # catch_n_latflow_m3 = catch_lat_inflow_m3 - catch_lat_outflow_m3
        # self.out_accu_n_latflow_m3_tss.sample(catch_n_latflow_m3)

        # Percolation
        percol_z0_m3 = percolation[0] * cellarea() / 1000  # m3
        percol_z1_m3 = percolation[1] * cellarea() / 1000  # m3
        percol_z2_m3 = percolation[2] * cellarea() / 1000  # m3
        percol_z3_m3 = percolation[3] * cellarea() / 1000  # m3
        accu_out_percol_m3 = accuflux(self.ldd_subs, percol_z3_m3)
        out_percol_m3 = areatotal(accu_out_percol_m3, self.outlet_multi)
        self.out_percol_bsmt_m3_tss.sample(out_percol_m3)

        # Evapotranspiration
        etp_z0_m3 = etp[0] * cellarea() / 1000  # m3
        # self.report(etp_z0_m3, 'az0ETP')
        etp_z1_m3 = etp[1] * cellarea() / 1000
        # self.report(etp_z1_m3, 'az1ETP')
        etp_z2_m3 = etp[2] * cellarea() / 1000
        etp_z3_m3 = etp[3] * cellarea() / 1000
        # self.report(etp_z2_m3, 'az2ETP')
        etp_m3 = etp_z0_m3 + etp_z1_m3 + etp_z2_m3 + etp_z3_m3
        accu_out_etp_m3 = accuflux(self.ldd_subs, etp_m3)
        out_etp_m3 = areatotal(accu_out_etp_m3, self.outlet_multi)
        self.out_etp_m3_tss.sample(out_etp_m3)

        # Baseflow discharge (basement)
        # accu_baseflow = accuflux(self.ldd_subs, self.baseflow * cellarea() / 1000)  # m3
        # out_baseflow_m3 = areatotal(accu_baseflow, self.outlet_multi)
        # self.out_baseflow_m3_tss.sample(out_baseflow_m3)

        # Change in storage
        ch_storage_m3 = ch_storage[0] + ch_storage[1] + ch_storage[2] + ch_storage[3]
        accu_out_ch_storage_m3 = accuflux(self.ldd_subs, ch_storage_m3)
        out_ch_storage_m3 = areatotal(accu_out_ch_storage_m3, self.outlet_multi)
        self.out_ch_storage_m3_tss.sample(out_ch_storage_m3)

        # GLOBAL Water
        global_mb_water = (tot_rain_m3 -
                           out_etp_m3 - out_runoff_m3 -
                           out_percol_m3 + catch_n_latflow_m3 -  # of_latflow_m3 -
                           out_drain_m3 -  # o_drain_z1_m3 -
                           out_ch_storage_m3  # - out_baseflow_m3
                           )

        self.global_mb_water_tss.sample(global_mb_water)

        cell_vol_z0 = self.theta[0] * self.layer_depth[0] * cellarea() / 1000
        cell_vol_z1 = self.theta[1] * self.layer_depth[1] * cellarea() / 1000
        cell_vol_z2 = self.theta[2] * self.layer_depth[2] * cellarea() / 1000
        cell_vol_z3 = self.theta[3] * self.layer_depth[3] * cellarea() / 1000
        cell_vol_tot_m3 = cell_vol_z0 + cell_vol_z1 + cell_vol_z2 + cell_vol_z3
        vol_tot_m3 = accuflux(self.ldd_subs, cell_vol_tot_m3)
        multi_vol_tot_m3 = areatotal(vol_tot_m3, self.outlet_multi)
        self.storage_m3_tss.sample(multi_vol_tot_m3)

        tot_vol_disch_m3 = getTotalDischarge(out_runoff_m3, outlet_cell_outflow_m3, out_drain_m3)

        # self.i_Q_m3_tss.sample(in_vol_disch_m3)
        # self.o_Q_m3_tss.sample(out_vol_disch_m3)
        self.tot_Q_m3_tss.sample(tot_vol_disch_m3)


        ######################
        # Pesticide Balance
        ######################

        # Applied ug on catchment
        accu_app_light = accuflux(self.ldd_subs, light_applied)  # ug
        # out_app_light = areatotal(accu_app_light, self.outlet_multi)  # ug
        out_app_light = areatotal(light_applied, self.is_catchment)  # ug
        self.out_app_L_tss.sample(out_app_light)

        # cum_appl_catch_light = accuflux(self.ldd_subs, self.cum_appl_g)
        # Todo: Check if below same result
        # cum_appl_catch_mg = upstream(self.ldd_subs, self.cum_appl_g)
        if self.TRANSPORT:
            # Degradation
            light_deg_tot = light_deg[0] + light_deg[1] + light_deg[2] + light_deg[3]
            out_deg_light = areatotal(light_deg_tot, self.is_catchment)
            z0_light_deg_catch = areatotal(light_deg_tot, self.is_catchment)
            self.cum_degZ0_L_g += z0_light_deg_catch
            self.out_degZ0_L_tss.sample(z0_light_deg_catch)
            self.cum_degZ0_L_g_tss.sample(self.cum_degZ0_L_g)

            # Volatilized
            out_volat_light = areatotal(light_volat, self.is_catchment)
            self.out_volat_L_tss.sample(out_volat_light)

            # Mass loss to run-off
            # Index: 0 <- light, Index: 1 <- heavy
            out_runoff_light = areatotal(mass_runoff[0], self.is_catchment)
            self.out_runoff_L_tss.sample(out_runoff_light)
            self.cum_roZ0_L_g += out_runoff_light
            self.cum_roZ0_L_g_tss.sample(self.cum_roZ0_L_g)

            # z0-mass leached
            out_leach_light_z0 = areatotal(light_leached[0], self.is_catchment)
            self.cum_lchZ0_L_g += out_leach_light_z0
            self.out_leachZ0_L_tss.sample(out_leach_light_z0)
            self.cum_lchZ0_L_g_tss.sample(self.cum_lchZ0_L_g)

            # z1-mass leached
            out_leach_light_z1 = areatotal(light_leached[1], self.is_catchment)
            self.out_leachZ1_L_tss.sample(out_leach_light_z1)

            # z3-mass leached = zero, if no basement percolation
            out_leach_light = areatotal(light_leached[2], self.is_catchment)
            self.out_leach_L_tss.sample(out_leach_light)

            # Artificial drained mass (layer z2)
            out_drain_light = areatotal(light_drained, self.is_catchment)
            self.out_drain_L_tss.sample(out_drain_light)
            # out_drain_heavy = areatotal(z1_heavy_drain, self.is_catchment)
            self.cum_adr_L_g += out_drain_light
            self.cum_adr_L_g_tss.sample(self.cum_adr_L_g)

            # Lateral flux at outlet cells
            outlet_cell_lightflux = (ligth_latflow_outlet[0] + ligth_latflow_outlet[1] +
                                     ligth_latflow_outlet[2] + ligth_latflow_outlet[3])
            outlet_cell_lightflux = areatotal(outlet_cell_lightflux, self.outlet_multi)  # Sum only outlet cells
            self.cum_latflux_L_g += outlet_cell_lightflux
            self.cum_latflux_L_g_tss.sample(self.cum_latflux_L_g)

            # Net loss to lateral flux (Required in MB)
            mass_latflux_light = (ligth_latflow[0] + ligth_latflow[1] +
                                  ligth_latflow[2] + ligth_latflow[3])

            net_latflux_light = areatotal(mass_latflux_light, self.is_catchment)  # Needed for MB
            # self.out_latflux_L_tss.sample(net_latflux_light)

            # Baseflow flux
            # out_baseflow_light = areatotal(baseflow_light, self.is_catchment)
            # self.out_baseflow_L_tss.sample(out_baseflow_light)

            # Total mass export
            outlet_light_export = (out_runoff_light + out_drain_light + outlet_cell_lightflux  # +
                                   # out_baseflow_light
                                   )
            self.cum_exp_L_g += outlet_light_export
            self.cum_exp_L_g_tss.sample(self.cum_exp_L_g)

            # Change in mass storage
            ch_storage_light = (ch_storage_light[0] + ch_storage_light[1] +
                                ch_storage_light[2] + ch_storage_light[3])

            out_ch_storage_light = areatotal(ch_storage_light, self.is_catchment)
            self.out_chstorage_L_tss.sample(out_ch_storage_light)

            # TODO: still need to add outlet cells to MB
            light_mb_pest = (out_app_light - out_deg_light -
                             out_volat_light - out_runoff_light -
                             out_leach_light -
                             out_drain_light + net_latflux_light -  #
                             # out_baseflow_light -
                             out_ch_storage_light)
            self.global_mb_pest_tss.sample(light_mb_pest)

        self.jd_cum += self.jd_dt  # updating JDcum, currently dt = 1 day

        # Analysis (NASH Discharge)
        q_obs = timeinputscalar('q_obs_m3day.tss', nominal("outlet_true"))
        rest_obs = tot_rain_m3 / q_obs
        self.rest_obs_tss.sample(rest_obs)

        self.q_obs_cum += ifthenelse(q_obs >= 0, q_obs, 0)
        self.days_cum += ifthenelse(q_obs >= 0, scalar(1), scalar(0))  # Total days with data
        self.q_sim_cum += ifthenelse(q_obs >= 0, tot_vol_disch_m3, 0)
        self.q_sim_ave = self.q_sim_cum / self.days_cum
        self.q_diff += ifthenelse(q_obs >= 0, (tot_vol_disch_m3 - q_obs) ** 2, 0)
        self.q_var += ifthenelse(q_obs >= 0, (q_obs - 260.07) ** 2, 0)  # Mean discharge of data range = 260.07 m3/day
        nash_q = 1 - (self.q_diff / self.q_var)
        self.nash_q_tss.sample(nash_q)

        self.q_obs_cum_tss.sample(self.q_obs_cum)
        self.q_sim_cum_tss.sample(self.q_sim_cum)
        self.q_sim_ave_tss.sample(self.q_sim_ave)

        self.rain_cum_m3 += tot_rain_m3
        self.rain_obs_cum_tss.sample(self.rain_cum_m3)

        # Mass balance components
        self.tot_runoff_m3 += ifthenelse(q_obs >= 0, out_runoff_m3, 0)
        self.tot_runoff_m3_tss.sample(self.tot_runoff_m3)
        self.tot_perc_z3_m3 += ifthenelse(q_obs >= 0, out_percol_m3, 0)  # <- == 0
        self.tot_perc_z3_m3_tss.sample(self.tot_perc_z3_m3)
        self.tot_etp_m3 += ifthenelse(q_obs >= 0, out_etp_m3, 0)
        self.tot_etp_m3_tss.sample(self.tot_etp_m3)
        # self.tot_baseflow_m3 += ifthenelse(q_obs >= 0, out_baseflow_m3, 0)
        # self.tot_baseflow_m3_tss.sample(self.tot_baseflow_m3)

        self.tot_drain_m3 += ifthenelse(q_obs >= 0, out_drain_m3, 0)  # o_drain_z1_m3
        self.tot_accu_drain_m3_tss.sample(self.tot_drain_m3)

        # self.tot_ilf_m3 += ifthenelse(q_obs >= 0, outlet_lat_inflow_m3, 0)
        # self.tot_accu_i_latflow_m3_tss.sample(self.tot_ilf_m3)
        self.tot_olf_m3 += ifthenelse(q_obs >= 0, outlet_cell_outflow_m3, 0)
        self.tot_accu_o_latflow_m3_tss.sample(self.tot_olf_m3)
        # self.tot_nlf_m3 += ifthenelse(q_obs >= 0, n_latflow_m3, 0)
        # self.tot_accu_n_latflow_m3_tss.sample(self.tot_nlf_m3)
        # self.tot_of_m3 += ifthenelse(q_obs >= 0, of_latflow_m3, 0)
        # self.tot_accu_of_latflow_m3_tss.sample(self.tot_of_m3)

        # Analysis Pest
        if self.PEST:
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
            reportNashConcComposites(self,
                                     soils_north['ave_conc'],
                                     soils_valley['ave_conc'],
                                     soils_south['ave_conc'])

            # reportNashDeltaComposites()

    def postmcloop(self):
        pass
        # names = ["q"]  # Discharge, Nash_Discharge
        # mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        # aguila --timesteps=[170,280,1] q-ave q-var outlet_true.map
        # percentiles = [0.25, 0.5, 0.75]
        # mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())
        # aguila --quantiles=[0.25,0.75,0.25] --timesteps=[170,280,1] q


# Visualization
# aguila 1\at0dC000.177 1\at1dC000.177
# aguila --scenarios='{1,2}' --multi=1x4  --timesteps=[175,179,1] aLEACH aLEACHz aLF aLFz
# aguila --scenarios='{1}'  --timesteps=[100,280,1] az0dC az1dC az2dC
# aguila --scenarios='{1}'  --timesteps=[1,280,1] aHeight aRDtot aCrop aPotETP akcb akcb1 akcmax
#  aguila --scenarios='{1}'  --timesteps=[1,360,1] aHeight aRDtot aCrop akcb aPotTRA aPotEVA
#  aguila --scenarios='{1,2,3}' --multi=1x4 --timesteps=[2,300,2] athz0 athz1 athz2 athz3
#  aguila --scenarios='{1,2,3}' --multi=1x4  --timesteps=[1,300,1] athz0 athz1 athz2 athz3
#  aguila --scenarios='{1,2,3}' --multi=1x4 --timesteps=[1,300,1] aObj1 aObj2
#  aguila --scenarios='{1,2,3}' thFCz2 thSATz2
#  aguila --scenarios='{1}' --timesteps=[2,300,2] aROm3 athz0 athz1 athz2 athz3


# Time series
# aguila 1\res_nash_q_m3.tss 6\res_nash_q_m3.tss
# aguila 1\resW_accStorage_m3.tss
# aguila 1\resM_norCONC.tss 1\resM_valCONC.tss 1\resM_souCONC.tss

nrOfSamples = int(runs)  # Samples are each a MonteCarlo realization
firstTimeStep = start_jday()  # 166 -> 14/03/2016
nTimeSteps = 300  # 360
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
