# -*- coding: utf-8 -*-

# from time import *
import time
from datetime import datetime
from hydro import *
from pesti import *
from pesti_v2 import *

from pcraster._pcraster import *
from pcraster.framework import *
import os

print(os.getcwd())

global state
state = 0


def get_state(old_state):
    new_state = old_state + 1
    global state
    state = new_state
    return state


class BeachModel(DynamicModel, MonteCarloModel):
    def setDebug(self):
        pass

    def __init__(self, cloneMap):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        setclone(cloneMap)

        # dem = self.dem

    def premcloop(self):
        self.PEST = True
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

        """
        Landscape Maps
        """
        self.dem = self.readmap("dem_slope")  # 192 - 231 m a.s.l
        self.dem_route = self.readmap("dem_ldd")  # To route surface run-off
        self.datum_depth = (self.dem - mapminimum(self.dem)) * scalar(10 ** 3)  # mm

        # self.ldd_surf = lddcreate(self.dem_route, 1e31, 1e31, 1e31, 1e31)  # To route runoff
        out_burn = readmap("dem_ldd_burn2")
        self.is_catchment = defined(out_burn)
        # self.ldd_subs = lddcreate(self.dem, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI
        self.ldd_subs = lddcreate(out_burn, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI

        self.zero_map = out_burn - out_burn  # Zero map to generate scalar maps
        self.mask = out_burn / out_burn

        # self.outlet = self.readmap("outlet_multi")
        self.outlet_multi = self.readmap("out_multi_nom")  # Multi-outlet with 0 or 1
        self.north_wk = nominal(self.readmap("north_nom"))
        self.valley_wk = nominal(self.readmap("valley_nom"))
        self.south_wk = nominal(self.readmap("south_nom"))

        self.landuse = self.readmap("landuse")

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
        self.out_percol_z2_m3_tss = TimeoutputTimeseries("resW_accPercol_z2_m3", self, nominal("outlet_true"),
                                                         noHeader=False)
        self.tot_perc_z2_m3_tss = TimeoutputTimeseries("resW_totPercol_z2_m3", self, nominal("outlet_true"),
                                                       noHeader=False)
        # ETP
        self.out_etp_m3_tss = TimeoutputTimeseries("resW_accEtp_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_etp_m3_tss = TimeoutputTimeseries("resW_totEtp_m3", self, nominal("outlet_true"), noHeader=False)

        self.out_baseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", self, nominal("outlet_true"),
                                                        noHeader=False)
        self.tot_baseflow_m3_tss = TimeoutputTimeseries("resW_totBaseflow_m3", self, nominal("outlet_true"),
                                                        noHeader=False)
        # LF Drainage
        self.out_accu_o_drain_m3_tss = TimeoutputTimeseries("resW_o_accDrain_m3", self, nominal("outlet_true"),
                                                            noHeader=False)
        self.tot_accu_drain_m3_tss = TimeoutputTimeseries("resW_o_totDrain_m3", self, nominal("outlet_true"),
                                                          noHeader=False)
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
        self.out_delta_tss = TimeoutputTimeseries("out_delta", self, nominal("outlet_true"), noHeader=False)
        self.out_mass_tss = TimeoutputTimeseries("out_mass", self, nominal("outlet_true"), noHeader=False)

        self.out_app_L_tss = TimeoutputTimeseries("resM_accAPP_L", self, nominal("outlet_true"), noHeader=False)
        self.out_volat_L_tss = TimeoutputTimeseries("resM_accVOL_L", self, nominal("outlet_true"), noHeader=False)
        self.out_runoff_L_tss = TimeoutputTimeseries("resM_accRO_L", self, nominal("outlet_true"), noHeader=False)
        self.out_leach_L_tss = TimeoutputTimeseries("resM_accDP_L", self, nominal("outlet_true"), noHeader=False)
        self.out_drain_L_tss = TimeoutputTimeseries("resM_accADR_L", self, nominal("outlet_true"), noHeader=False)
        self.out_latflux_L_tss = TimeoutputTimeseries("resM_accLF_L", self, nominal("outlet_true"), noHeader=False)
        self.out_baseflow_L_tss = TimeoutputTimeseries("resM_accBF_L", self, nominal("outlet_true"), noHeader=False)
        self.out_chstorage_L_tss = TimeoutputTimeseries("resM_accCHS_L", self, nominal("outlet_true"), noHeader=False)

        # Transects
        self.north_conc_tss = TimeoutputTimeseries("resM_norCONC", self, ordinal("north_ave"), noHeader=False)
        self.valley_conc_tss = TimeoutputTimeseries("resM_valCONC", self, ordinal("valley_ave"), noHeader=False)
        self.south_conc_tss = TimeoutputTimeseries("resM_souCONC", self, ordinal("south_ave"), noHeader=False)

        self.north_dC13_tss = TimeoutputTimeseries("resM_nordC13", self, ordinal("north_ave"), noHeader=False)
        self.valley_dC13_tss = TimeoutputTimeseries("resM_valdC13", self, ordinal("valley_ave"), noHeader=False)
        self.south_dC13_tss = TimeoutputTimeseries("resM_soudC13", self, ordinal("south_ave"), noHeader=False)

        # Analysis
        # This is 'q' as time series.
        self.i_Q_m3_tss = TimeoutputTimeseries("resN_i_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.o_Q_m3_tss = TimeoutputTimeseries("resN_o_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.n_Q_m3_tss = TimeoutputTimeseries("resN_n_accVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.q_obs_tot_tss = TimeoutputTimeseries("resN_q_obs_tot_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # Equivalent to net_Q
        self.rain_obs_tot_tss = TimeoutputTimeseries("resN_rain_obs_tot_m3", self, nominal("outlet_true"),
                                                     noHeader=False)  # Equivalent to net_Q
        self.rest_obs_tss = TimeoutputTimeseries("resN_restit_obs_m3", self, nominal("outlet_true"),
                                                 noHeader=False)  # = rain/q_obs
        self.q_sim_tot_tss = TimeoutputTimeseries("resN_q_sim_tot_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # This is 'Nash_q' as time series.
        self.nash_q_tss = TimeoutputTimeseries("resN_nash_q_m3", self, nominal("outlet_true"),
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
        # coefficient to calibrate/test Ksat1
        self.s1 = scalar(self.ini_param.get("s1"))  # mm/day
        self.s2 = scalar(self.ini_param.get("s2"))  # coefficient to calibrate Ksat2

        # First-sensitivity:
        """ Physical parameters """
        self.c1 = scalar(self.ini_param.get("c1"))  # subsurface flow coefficient
        self.c2 = scalar(self.ini_param.get("c2"))  # not used (second layer)
        self.drain_coef = scalar(self.ini_param.get("drain_coef"))  # drainage coefficient
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
        self.alpha_iso = scalar(self.ini_param.get("alpha_iso"))  # 1 is no fractionation

        # First run will return state = 1
        m_state = get_state(state)
        self.PERCOL = False  # z2 deep percolation (DP)
        self.ADLF = True
        self.c_adr = 0.50
        self.k_g = 1500  # [days] = 4 yrs
        z2_factor = 0.9
        if m_state == 1:
            pass  # c1, c2 = 0.15, d1 = 0.8
        elif m_state == 2:
            pass
        self.gw_factor = 1 - z2_factor

        self.z0 = self.zero_map + 10  # mm -> 1cm
        self.z1 = self.zero_map + 300  # mm
        self.z2 = (self.datum_depth + 610 - self.z0 - self.z1)  # * z2_factor  # mm (300mm at outlet) * z2_factor
        self.tot_depth = self.z0 + self.z1 + self.z2
        self.smp_depth = self.z0

        """
        Hydro Maps
        """
        # Initial moisture (Final from model v1, Sept 30, 2016)
        self.theta_z0 = readmap('ini_theta_z0')  # map of initial soil moisture in top layer (-)
        self.theta_z1 = readmap('ini_theta_z1')
        self.theta_z2 = max(readmap("ini_theta_z2"), scalar(0.1)) * z2_factor + scalar(0.6) * self.gw_factor

        # Need initial states to compute change in storage after each run
        self.theta_z0_ini = self.theta_z0
        self.theta_z1_ini = self.theta_z1
        self.theta_z2_ini = self.theta_z2

        """
        Pesticides Maps
        """
        if (self.PEST):
            # Application days
            self.app_days = [177, 197, 238]
            # Mass
            # in ug = conc. (ug/g soil) * density (g/cm3) * (10^6 cm3/m3)*(1 m/10^3 mm)* depth_layer(mm) * cellarea(m2)
            # g = ug * 10e-06
            self.smback_z0 = (self.zero_map + 0.06) * self.p_b * scalar(
                10 ** 6 / 10 ** 3) * self.z0 * cellarea() * 10e-06  # Based on detailed soils
            self.smback_z1 = (self.zero_map + 0.03) * self.p_b * scalar(
                10 ** 6 / 10 ** 3) * self.z1 * cellarea() * 10e-06  # Based on detailed soils
            self.smback_z2 = (self.zero_map + 0.00001) * self.p_b * scalar(
                10 ** 6 / 10 ** 3) * self.z2 * cellarea() * 10e-06  # Assumed

            # Carbon Delta (Background)
            # Assumed theoretical max @99% deg Streitwieser Semiclassical Limits
            self.delta_z0 = self.zero_map - 23.7
            self.delta_z0_ini = self.delta_z0
            self.delta_z1 = self.zero_map - 23.7
            self.delta_z1_ini = self.delta_z1
            self.delta_z2 = self.zero_map - 23.7
            self.delta_z2_ini = self.delta_z2

            self.lightback_z0 = self.smback_z0 / (1 + self.r_standard * (self.delta_z0 / 1000 + 1))
            self.lightback_z1 = self.smback_z1 / (1 + self.r_standard * (self.delta_z1 / 1000 + 1))
            self.lightback_z2 = self.smback_z2 / (1 + self.r_standard * (self.delta_z2 / 1000 + 1))

            self.heavyback_z0 = self.smback_z0 - self.lightback_z0
            self.heavyback_z1 = self.smback_z1 - self.lightback_z1
            self.heavyback_z2 = self.smback_z2 - self.lightback_z2

            # self.report(self.smback_z0, 'a_tback')
            # self.report(self.lightback_z0, "a_lback")
            # self.report(self.heavyback_z0, "a_hback")

            # Applications Mass
            # Product concentration (active ing.)
            double = scalar(2.0)  # ~ Dosage for corn when growing beet
            d_gold = scalar(915)  # g/L S-met # * 10 ** 6  # ug/L S-met
            m_gold = scalar(960)  # g/L S-met # * 10 ** 6  # ug/L

            # Dosages # L/Ha * 1Ha/1000m2 = L/m2
            d_beet = None
            d_corn = scalar(2.1) * 1 / 10 ** 4  # 2.1 L/Ha * 1 Ha / 10000 m2
            m_beet = scalar(0.6) * 1 / 10 ** 4 * double
            m_corn = scalar(2.0) * 1 / 10 ** 4
            # (Likely larger dosage, early in the season)
            m_beet_Friess = scalar(0.6) * 1 / 10 ** 4 * (double + 1)  # 0.6 L/Ha * 1 Ha / 10000 m2
            m_beet_Mathis = scalar(0.6) * 1 / 10 ** 4 * (double + 1)  #

            # Assign dosages based on Farmer-Crop combinations [ug/m2]
            fa_cr = readmap("crop_burn")  # Contains codes to assign appropriate dosage
            app_conc = (  # [ug/m2] -> [g/m2]
                ifthenelse(fa_cr == 1111,  # 1111 (Friess, Beet)
                           m_beet_Friess * m_gold * self.mask,
                           ifthenelse(fa_cr == 1122,  # 1112 (Friess-Corn),
                                      m_corn * m_gold * self.mask,
                                      ifthenelse(fa_cr == 1212,  # 1212 (Speich-Corn),
                                                 m_corn * m_gold * self.mask,
                                                 ifthenelse(fa_cr == 1312,  # 1312 (Mahler-Corn),
                                                            m_corn * m_gold * self.mask,
                                                            ifthenelse(fa_cr == 1412,  # 1412 (Schmitt-Corn)
                                                                       d_corn * d_gold * self.mask,
                                                                       ifthenelse(fa_cr == 1511,  # 1511 (Burger-Beet)
                                                                                  m_beet * m_gold * self.mask,
                                                                                  # 1711 (Mathis-Beet),
                                                                                  ifthenelse(fa_cr == 1711,
                                                                                             m_beet_Mathis * m_gold * self.mask,
                                                                                             # 1611 (Kopp-Beet)
                                                                                             ifthenelse(
                                                                                                 fa_cr == 1611,
                                                                                                 m_beet * m_gold * self.mask,
                                                                                                 0 * self.mask))))))))
            )
            # Pesticide applied (ug->g) on Julian day 177 (March 25, 2016).
            # March 26th, Friess and Mathis
            self.app1 = ifthenelse(fa_cr == 1111, 1 * app_conc * cellarea(),
                                   # 1111 (Friess, Beet), 1112 (Friess-Corn),
                                   ifthenelse(fa_cr == 1112, 1 * app_conc * cellarea(),
                                              ifthenelse(fa_cr == 1711, 1 * app_conc * cellarea(),  # 1711 (Mathis-Beet)
                                                         0 * app_conc * cellarea())))
            # Pesticide applied (ug->g) on Julian day 197 (April 14, 2016).
            # April 14, Kopp and Burger
            self.app2 = ifthenelse(fa_cr == 1511, 1 * app_conc * cellarea(),  # 1511 (Burger-Beet)
                                   ifthenelse(fa_cr == 1611, 1 * app_conc * cellarea(),  # 1611 (Kopp-Beet),
                                              0 * app_conc * cellarea()))

            # Pesticide applied (ug->g) on Julian day 238 (May 25, 2016).
            # May 25, Schmidt and Speich, and (out of transect): Friess and Mahler
            # Note: Speich and Friess could be 1 week later.
            self.app3 = ifthenelse(fa_cr == 1112, 1 * app_conc * cellarea(),  # 1112 (Friess-Corn)
                                   ifthenelse(fa_cr == 1212, 1 * app_conc * cellarea(),  # 1212 (Speich-Corn),
                                              ifthenelse(fa_cr == 1412, 1 * app_conc * cellarea(),
                                                         # 1412 (Schmitt-Corn),
                                                         ifthenelse(fa_cr == 1312, 1 * app_conc * cellarea(),
                                                                    # 1312 (Mahler-Corn)
                                                                    0 * app_conc * cellarea()))))

            # Applications delta
            # Use map algebra to produce a initial signature map,
            # ATT: Need to do mass balance on addition of new layer.
            # where app1 > 0, else background sig. (plots with no new mass will be 0)
            # where app1 > 0, else background sig. (plots with no new mass will be 0)
            self.app1delta = ifthenelse(self.app1 > 0, scalar(-32.3), scalar(-23.7))
            self.app2delta = ifthenelse(self.app2 > 0, scalar(-32.3), scalar(-23.7))
            self.app3delta = ifthenelse(self.app3 > 0, scalar(-32.3), scalar(-23.7))

            # Assign background as initial
            self.pestmass_z0 = self.smback_z0  # ug
            self.pestmass_z0_ini = self.pestmass_z0
            self.pestmass_z1 = self.smback_z1  # ug
            self.pestmass_z1_ini = self.pestmass_z1
            self.pestmass_z2 = self.smback_z2  # ug
            self.pestmass_z2_ini = self.pestmass_z2

            self.lightmass_z0 = self.lightback_z0  # ug
            self.lightmass_z0_ini = self.lightmass_z0
            self.lightmass_z1 = self.lightback_z1  # ug
            self.lightmass_z1_ini = self.lightmass_z1
            self.lightmass_z2 = self.lightback_z2  # ug
            self.lightmass_z2_ini = self.lightmass_z2

            self.heavymass_z0 = self.heavyback_z0  # ug
            self.heavymass_z0_ini = self.heavymass_z0
            self.heavymass_z1 = self.heavyback_z1  # ug
            self.heavymass_z1_ini = self.heavymass_z1
            self.heavymass_z2 = self.heavyback_z2  # ug
            self.heavymass_z2_ini = self.heavymass_z2

            # Cumulative maps
            # self.light_ini_storage_ug = self.lightmass_z0_ini + self.lightmass_z1_ini + self.lightmass_z2_ini
            self.cum_appl_ug = self.zero_map
            self.cum_runoff_ug = self.zero_map
            self.cum_leached_ug_z0 = self.zero_map
            self.cum_leached_ug_z1 = self.zero_map
            self.cum_leached_ug_z2 = self.zero_map
            self.cum_latflux_ug_z0 = self.zero_map
            self.cum_latflux_ug_z1 = self.zero_map
            self.cum_latflux_ug_z2 = self.zero_map
            self.cum_baseflx_ug_z2 = self.zero_map

        """
        Temperature maps and params
        """
        self.lag = scalar(0.8)  # lag coefficient (-), 0 < lag < 1; -> in SWAT, lag = 0.80
        # Generating initial surface temp map (15 deg is arbitrary)
        self.temp_z0_fin = self.zero_map + 15
        self.temp_z1_fin = self.zero_map + 15
        self.temp_z2_fin = self.zero_map + 15
        self.temp_surf_fin = self.zero_map + 15

        # Maximum damping depth (dd_max)
        # The damping depth (dd) is calculated daily and is a function of max. damping depth (dd_max), (mm):
        self.dd_max = (scalar(2500) * self.p_b) / (self.p_b + 686 * exp(-5.63 * self.p_b))

        # TODO
        # Average Annual air temperature (celcius - Layon!! Not Alteckendorf yet!!)
        self.temp_ave_air = scalar(12.2)  # 12.2 is for Layon

        """ 
        Simulation start time: Oct 1st, 2015 
        """
        yy = scalar(2015)
        mm = scalar(10)
        dd = scalar(1)

        date_factor = 1
        if (100 * yy + mm - 190002.5) < 0:
            date_factor = -1

        # simulation start time in JD (Julian Day)
        self.jd_start = 367 * yy - rounddown(7 * (yy + rounddown((mm + 9) / 12)) / 4) + rounddown(
            (275 * mm) / 9) + dd + 1721013.5 - 0.5 * date_factor
        self.jd_cum = scalar(0)
        self.jd_dt = scalar(1)  # Time step size (days)

        # Analysis
        self.q_diff = 0
        self.q_var = 0
        self.q_obs_tot = 0
        self.q_sim_tot = 0  # net total disch

        self.rain_cum_m3 = 0  # Rainfall

        self.tot_drain_m3 = 0  # drainage z1
        self.tot_nlf_m3 = 0
        self.tot_ilf_m3 = 0  # upstream inflow
        self.tot_olf_m3 = 0  # downstream outflow

        self.tot_of_m3 = 0  # Overflow due to LF sat capacity reached

        self.tot_etp_m3 = 0
        self.tot_baseflow_m3 = 0
        self.tot_perc_z2_m3 = 0
        self.tot_runoff_m3 = 0

        # Stochastic / test parameters
        print("state:", m_state)
        # self.report(self.c1, 'c1')
        # self.report(self.c2, 'c2')
        # self.report(self.drain_coef, 'd1')
        # self.report(self.z0, 'z0')
        # self.report(self.s1, 's1')
        # if self.PEST:
        #     self.report(self.pestmass_z0, 'z0mini')
        #     self.report(self.delta_z0_ini, 'z0dCini')
        #     self.report(self.app1, 'z0app1')
        #     self.report(self.app1delta, 'z0app1dC')

    def dynamic(self):
        jd_sim = self.jd_start + self.jd_cum
        # Tells which row in crop table to select, based on the nominal.landuse value
        fields = timeinputscalar('landuse.tss', nominal(self.landuse))

        # SEE: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_timeinput....html?highlight=timeinputscalar
        # returns value of land-use field (i.e. n = 22), per time step. (Layon)
        # So, at dt = 1
        # fields's values: 98 98	98	98	98	98	98	98	98	98	98	98	13	98	15	99	99	99	99	99	99	98

        " Crop Parameters "
        # SEE: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_lookup.html?highlight=lookupscalar
        setglobaloption('matrixtable')  # allows lookupscalar to read more than 2 expressions.
        crop_type = lookupscalar('croptable.tbl', 1, fields)  # (table, col-value, row-value)
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

        p_tab = lookupscalar('croptable.tbl', 16,
                             fields)  # depletable theta before water stress (Allen1998, Table no.22)
        # p_tab (-) according to Allen 1998, Table 22 (now using FAO source)
        # Sugar beet = 0.55
        # Corn = 0.55
        # Grazing Pasture = 0.6
        # Spring Wheat = 0.55
        # Winter Wheat = 0.55
        # Apple trees = 0.5

        """ Soil physical parameters """
        # Saturated moisture capacity is equal for depth0 and depth1
        theta_sat_z0z1 = lookupscalar('croptable.tbl', 17, fields)  # saturated moisture of the first layer # [-]
        theta_fcap_z0z1 = lookupscalar('croptable.tbl', 18,
                                       fields)  # field capacity of 1st layer (equal for D0 and k=1)
        theta_sat_z2 = lookupscalar('croptable.tbl', 19, fields)  # saturated moisture of 2nd layer
        theta_fcap_z2 = lookupscalar('croptable.tbl', 20, fields)  # field capacity of the 2nd layer
        theta_wp = lookupscalar('croptable.tbl', 21, fields)  # wilting point moisture
        # k_sat_z0z1 = lookupscalar('croptable.tbl', 22, fields)  # saturated conductivity of the first layer

        k_sat_z2 = lookupscalar('croptable.tbl', 23, fields)  # saturated conductivity of the second layer
        k_sat_z2 *= self.s2
        CN2 = lookupscalar('croptable.tbl', 24, fields)  # curve number of moisture condition II

        k_sat_z0z1 = timeinputscalar('ksats.tss', nominal(self.landuse))
        k_sat_z0z1 *= self.s1

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
        all_stages = len_grow_stage_ini + len_dev_stage + len_mid_stage + len_end_stage

        # updating of sowing date by land use
        sow_yy = ifthenelse(jd_sim < jd_sow + all_stages, sow_yy,
                            ifthenelse(jd_sim < jd_sow + all_stages + 365, sow_yy + 1,
                                       ifthenelse(jd_sim < jd_sow + all_stages + 730, sow_yy + 2,
                                                  ifthenelse(jd_sim < jd_sow + all_stages + 1095, sow_yy + 3,
                                                             ifthenelse(jd_sim < jd_sow + all_stages + 1460,
                                                                        sow_yy + 4,
                                                                        scalar(0))))))

        # Update sowing date / plant date
        jd_plant = convertJulian(sow_yy, sow_mm, sow_dd)

        jd_dev = jd_plant + len_grow_stage_ini
        jd_mid = jd_dev + len_dev_stage
        jd_late = jd_mid + len_mid_stage
        jd_end = jd_late + len_end_stage
        LAIful = max_LAI + 0.5

        # calculation of crop height
        height = ifthenelse(crop_type > scalar(1), max_height,
                            ifthenelse(jd_sim < jd_plant, scalar(0),
                                       ifthenelse(jd_sim < jd_mid + 0.5 * len_mid_stage,
                                                  max_height * (jd_sim - jd_plant) / (
                                                      len_grow_stage_ini + len_dev_stage + 0.5 * len_mid_stage),
                                                  ifthenelse(jd_sim < jd_end, max_height,
                                                             0))))
        # calculation of root depth
        # TODO: Check first argument is correct, find documentation for root depth & height
        root_depth = ifthenelse(crop_type > scalar(1), max_root_depth,
                                ifthenelse(jd_sim < jd_plant, scalar(0),
                                           ifthenelse(jd_sim < jd_mid + len_mid_stage / 2,
                                                      max_root_depth * (jd_sim - jd_plant) / (
                                                          len_grow_stage_ini + len_dev_stage + len_mid_stage / 2),
                                                      ifthenelse(jd_sim < jd_end, max_root_depth,
                                                                 scalar(0)))))

        # root dispersal for each soil layer (z)
        root_depth_z0 = ifthenelse(root_depth > self.z0, self.z0, root_depth)
        root_depth_z1 = ifthenelse(root_depth < self.z1, scalar(0),
                                   ifthenelse(root_depth < self.z1 + self.z0, root_depth - self.z0, self.z1))
        root_depth_z2 = ifthenelse(root_depth <= self.z0 + self.z1, scalar(0),
                                   ifthenelse(root_depth < self.tot_depth, root_depth - self.z1 - self.z0, self.z2))

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

        # Not in use for water balance, but used to estimate surface temp due to bio-cover.
        frac_soil_cover = etp_dict["f"]

        bio_cover = getBiomassCover(self, frac_soil_cover)

        ######################################################################################
        # Mixing layer: depth z0
        # State functions
        #########################
        # Temperature, z0
        temp_dict_z0 = getLayerTemp(self, 0, bio_cover, temp_bare_soil)
        self.temp_surf_fin = temp_dict_z0["temp_surface"]
        self.temp_z0_fin = temp_dict_z0["temp_layer"]

        #########################
        # Moisture, z0
        z0_moisture = getLayerMoisture(self, 0,
                                       precip, theta_wp, CN2, crop_type,
                                       jd_sim, jd_dev, jd_mid, jd_end, len_dev_stage,
                                       root_depth, pot_evapor, pot_transpir, depletable_water,
                                       k_sat_z0z1, root_depth_z0,
                                       theta_fcap_z0z1, theta_sat_z0z1)
        percolation_z0 = z0_moisture["percolate"]
        sat_excess_z0 = z0_moisture["satex"]
        # TODO: verify amounts lost to sub-layer are correct
        tot_percolation_z0 = percolation_z0 + sat_excess_z0
        runoff_z0 = z0_moisture["runoff"]
        etp_z0 = z0_moisture["ETP"]
        # Report the discharge map
        # discharge = accuflux(self.ldd, runoff_z0)
        # self.report(discharge, "dt" + str(self.jd_cum) + "discharge.map")
        lat_netflow_z0 = z0_moisture["lat_flow"]
        lat_inflow_z0 = z0_moisture["upstream_lat_inflow"]
        lat_outflow_z0 = z0_moisture["cell_lat_outflow"]
        overflow_height_z0 = z0_moisture["overflow_height"]

        if self.PEST:
            #########################
            # Mass Transfer, z0

            # Applications
            mass_applied = self.zero_map
            light_applied = self.zero_map
            heavy_applied = self.zero_map
            light_volat = self.zero_map

            # z0_mass_volatilized = self.zero_map
            if self.currentTimeStep() in self.app_days:
                mass_applied = ifthenelse(self.currentTimeStep() == self.app_days[0],
                                          self.app1,
                                          ifthenelse(self.currentTimeStep() == self.app_days[1], self.app2,
                                                     ifthenelse(self.currentTimeStep() == self.app_days[2], self.app3,
                                                                scalar(0))))

                self.cum_appl_ug += mass_applied
                self.pestmass_z0 += mass_applied

                light_applied = ifthenelse(self.currentTimeStep() == self.app_days[0],
                                           mass_applied / (1 + self.r_standard * (self.app1delta / 1000 + 1)),
                                           ifthenelse(self.currentTimeStep() == self.app_days[1],
                                                      mass_applied / (
                                                          1 + self.r_standard * (self.app2delta / 1000 + 1)),
                                                      ifthenelse(self.currentTimeStep() == self.app_days[2],
                                                                 mass_applied / (
                                                                     1 + self.r_standard * (self.app3delta / 1000 + 1)),
                                                                 scalar(0))))
                heavy_applied = mass_applied - light_applied

                self.lightmass_z0 += light_applied
                self.heavymass_z0 += heavy_applied

                # Isotopes change due to applications
                # self.delta_z0 = ((self.heavymass_z0 / self.lightmass_z0 - self.r_standard) /
                #                      self.r_standard) * 1000  # [permille]

                if self.TRANSPORT:
                    # Mass & delta volatilized
                    light_volat = getVolatileMass(self, self.temp_air, theta_sat_z0z1, self.lightmass_z0, "LL",
                                                  rel_diff_model='option-1', sorption_model="linear",
                                                  gas=True, isotopes=True)
                    heavy_volat = getVolatileMass(self, self.temp_air, theta_sat_z0z1, self.heavymass_z0, "HH",
                                                  rel_diff_model='option-1', sorption_model="linear",
                                                  gas=True, isotopes=True)

                    light_volat = ifthenelse(mass_applied > 0, light_volat, scalar(0))
                    heavy_volat = ifthenelse(mass_applied > 0, heavy_volat, scalar(0))

                    self.lightmass_z0 -= light_volat
                    self.heavymass_z0 -= heavy_volat

                # Final delta after application day (minus volat)
                # delta_before = self.delta_z0
                self.delta_z0 = ((self.heavymass_z0 / self.lightmass_z0 - self.r_standard) /
                                 self.r_standard) * 1000  # [permille]
                # self.report(delta_before, 'at0dC')
                # self.report(self.delta_z0, 'az0dC')

            # self.report(mass_applied, 'aMapp')
            # self.report(light_applied, 'aMappL')
            # self.report(heavy_applied, 'aMappH')
            # self.report(self.lightmass_z0, 'aMz0LL')
            # self.report(light_volat, 'aMz0Vol')
            if self.TRANSPORT:
                # Mass & delta run-off (RO)
                light_runoff = getRunOffMass(self, theta_sat_z0z1, precip, runoff_z0,
                                             self.lightmass_z0,
                                             transfer_model="simple-mt", sorption_model="linear",
                                             gas=True)
                heavy_runoff = getRunOffMass(self, theta_sat_z0z1, precip, runoff_z0,
                                             self.heavymass_z0,
                                             transfer_model="simple-mt", sorption_model="linear",
                                             gas=True)

                self.lightmass_z0 -= light_runoff  # ug
                self.heavymass_z0 -= heavy_runoff

                # Mass & delta leached (Deep Percolation - DP)
                z0_light_leached = getLeachedMass(self, 0, theta_sat_z0z1, theta_fcap_z0z1,
                                                  precip,
                                                  z0_moisture["theta_after_percolate"],
                                                  self.lightmass_z0,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)
                z0_heavy_leached = getLeachedMass(self, 0, theta_sat_z0z1, theta_fcap_z0z1,
                                                  precip,
                                                  z0_moisture["theta_after_percolate"],
                                                  self.heavymass_z0,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)
                # self.report(z0_light_leached, 'aLEACH')
                # self.report(self.lightmass_z0, 'aLEACHz')

                self.lightmass_z0 -= z0_light_leached
                self.heavymass_z0 -= z0_heavy_leached

                # Mass (LF)
                z0_light_latflux = getLatMassFlux(self, 0, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.lightmass_z0,
                                                  sorption_model="linear", gas=True)
                z0_heavy_latflux = getLatMassFlux(self, 0, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.heavymass_z0,
                                                  sorption_model="linear", gas=True)

                self.lightmass_z0 += z0_light_latflux['net_mass_latflux']
                self.heavymass_z0 += z0_heavy_latflux['net_mass_latflux']

                # net_latflux = z0_light_latflux['net_mass_latflux']
                # self.report(net_latflux, 'aLF')
                # self.report(self.lightmass_z0, 'aLFz')

            # Degradation
            if self.DEGRADE:
                z0_deg_light = getMassDegradation(self, 0,
                                                  theta_sat_z0z1, theta_fcap_z0z1, theta_wp,
                                                  self.lightmass_z0, frac="L", sor_deg_factor=1)
                z0_deg_heavy = getMassDegradation(self, 0,
                                                  theta_sat_z0z1, theta_fcap_z0z1, theta_wp,
                                                  self.heavymass_z0, frac="H", sor_deg_factor=1)

                self.lightmass_z0 = z0_deg_light["mass_tot_new"]
                z0_light_deg = z0_deg_light.get("mass_deg_aq") + z0_deg_light.get("mass_deg_ads")
                self.heavymass_z0 = z0_deg_heavy["mass_tot_new"]

            # Change in mass storage after degradation - Pesticide Mass
            ch_storage_z0_light = self.lightmass_z0 - self.lightmass_z0_ini
            self.lightmass_z0_ini = self.lightmass_z0

            # Cumulative
            # tot_runoff = light_runoff + heavy_runoff
            # self.cum_runoff_ug += tot_runoff
            # self.cum_latflux_ug_z0 += z0_light_latflux["net_mass_latflux"]
            self.delta_z0 = ((self.heavymass_z0 / self.lightmass_z0 - self.r_standard) /
                             self.r_standard) * 1000  # [permille]
            # self.report(self.delta_z0, 'az0dC')
            # self.report(self.lightmass_z0, 'az0ML')

        # Update state variables
        # Change in storage - Moisture
        self.theta_z0 = z0_moisture["theta_final"]
        ch_storage_z0_m3 = (self.theta_z0 * self.z0 * cellarea() / 1000) - \
                           (self.theta_z0_ini * self.z0 * cellarea() / 1000)
        self.theta_z0_ini = self.theta_z0

        # self.theta_z0tss.sample(self.theta_z0)
        # self.water_balance_z0tss.sample(z0_moisture["balance"])

        # Get Transect concentrations
        # Observed conc. is betw 1 and 8 ug/g dry soil (on transect)
        # Observed conc. can reach 20 ug/g dry soil (on single plot on application)

        # Isotope mass balance
        cell_mass = self.lightmass_z0 + self.heavymass_z0
        cell_mass_delta = cell_mass * self.delta_z0

        north_tot_mass = areatotal(cell_mass, self.north_wk)
        north_dC13 = areatotal(cell_mass_delta, self.north_wk) / north_tot_mass
        valley_tot_mass = areatotal(cell_mass, self.valley_wk)
        valley_dC13 = areatotal(cell_mass_delta, self.valley_wk) / valley_tot_mass
        south_tot_mass = areatotal(cell_mass, self.south_wk)
        south_dC13 = areatotal(cell_mass_delta, self.south_wk) / south_tot_mass

        # Concentrations
        north_ave_mass = north_tot_mass / scalar(30)
        north_ave_conc = 10e6 * (north_ave_mass /
                                 (cellarea() * self.smp_depth)) * 1 / (self.p_b * 10e03)  # ug/g soil
        valley_ave_mass = valley_tot_mass / scalar(25)
        valley_ave_conc = 10e6 * (valley_ave_mass /
                                  (cellarea() * self.smp_depth)) * 1 / (self.p_b * 10e03)  # ug/g soil
        south_ave_mass = south_tot_mass / scalar(26)
        south_ave_conc = 10e6 * (south_ave_mass /
                                 (cellarea() * self.smp_depth)) * 1 / (self.p_b * 10e03)  # ug/g soil

        self.north_conc_tss.sample(north_ave_conc)  # 30 obs
        self.valley_conc_tss.sample(valley_ave_conc)  # 25 obs
        self.south_conc_tss.sample(south_ave_conc)  # 26 obs

        self.north_dC13_tss.sample(north_dC13)
        self.valley_dC13_tss.sample(valley_dC13)
        self.south_dC13_tss.sample(south_dC13)

        #######################################################################################
        # Layer z = 1  (Layer with artificial drainage)

        # State functions
        # Temperature
        temp_dict_z1 = getLayerTemp(self, 1, bio_cover, temp_bare_soil)
        self.temp_z1_fin = temp_dict_z1["temp_layer"]
        # Moisture
        z1_moisture = getLayerMoisture(self, 1,
                                       precip, theta_wp, CN2, crop_type,
                                       jd_sim, jd_dev, jd_mid, jd_end, len_dev_stage,
                                       root_depth, pot_evapor, pot_transpir, depletable_water,
                                       k_sat_z0z1, root_depth_z1,
                                       theta_fcap_z0z1, theta_sat_z0z1,
                                       percolate=percolation_z0, satex=sat_excess_z0,
                                       PERCOL_z2=True,
                                       ADLF=self.ADLF, c_adr=self.c_adr)
        percolation_z1 = z1_moisture["percolate"]
        drain_outflow_z1 = z1_moisture["drain_lat_outflow"]
        lat_netflow_z1 = z1_moisture["lat_flow"]
        lat_outflow_z1 = z1_moisture["cell_lat_outflow"]
        lat_inflow_z1 = z1_moisture["upstream_lat_inflow"]
        overflow_height_z1 = z1_moisture["overflow_height"]
        etp_z1 = z1_moisture["ETP"]

        if self.PEST:
            #########################
            # Mass Transfer, z1
            # Mass volatilized = not relevant @z1!
            # Mass runoff = not relevant @z1!
            # Mass & delta leached (Deep Percolation - DP, z1)
            if self.TRANSPORT:
                self.lightmass_z1 += z0_light_leached
                self.heavymass_z1 += z0_heavy_leached
                # z1_mass_leached = getLeachedMass(self, 1, theta_sat_z0z1,
                #                                  precip,
                #                                  percolation_z1,
                #                                  z1_moisture["theta_after_percolate"],
                #                                  sorption_model="linear")

                # Mass & delta leached (Deep Percolation - DP)
                z1_light_leached = getLeachedMass(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  percolation_z1,
                                                  z1_moisture["theta_after_percolate"],
                                                  self.lightmass_z1,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)
                z1_heavy_leached = getLeachedMass(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  percolation_z1,
                                                  z1_moisture["theta_after_percolate"],
                                                  self.heavymass_z1,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)

                self.lightmass_z1 -= z1_light_leached
                self.heavymass_z1 -= z1_heavy_leached
                # self.report(percolation_z1, 'aPercol')
                # self.report(z1_light_leached, 'aLEACH')
                # self.report(self.lightmass_z1, 'aLEACHz')

                # Mass drainage (z1)
                z1_light_drain = getDrainMassFlux(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.lightmass_z1,
                                                  sorption_model='linear', gas=True)
                z1_heavy_drain = getDrainMassFlux(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.heavymass_z1,
                                                  sorption_model='linear', gas=True)

                self.lightmass_z1 -= z1_light_drain
                self.heavymass_z1 -= z1_heavy_drain

                # Mass latflux (LF), z1
                z1_light_latflux = getLatMassFlux(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.lightmass_z1,
                                                  sorption_model="linear", gas=True)
                z1_heavy_latflux = getLatMassFlux(self, 1, theta_sat_z0z1, theta_fcap_z0z1,
                                                  self.heavymass_z1,
                                                  sorption_model="linear", gas=True)

                self.lightmass_z1 += z1_light_latflux['net_mass_latflux']
                self.heavymass_z1 += z1_heavy_latflux['net_mass_latflux']

            # Degradation
            if self.DEGRADE:
                z1_deg_light = getMassDegradation(self, 1,
                                                  theta_sat_z0z1, theta_fcap_z0z1, theta_wp,
                                                  self.lightmass_z1, frac="L", sor_deg_factor=1)
                z1_deg_heavy = getMassDegradation(self, 1,
                                                  theta_sat_z0z1, theta_fcap_z0z1, theta_wp,
                                                  self.heavymass_z1, frac="H", sor_deg_factor=1)

                self.lightmass_z1 = z1_deg_light["mass_tot_new"]
                z1_light_deg = z1_deg_light.get("mass_deg_aq") + z1_deg_light.get("mass_deg_ads")
                self.heavymass_z1 = z1_deg_heavy["mass_tot_new"]

            self.delta_z1 = ((self.heavymass_z1 / self.lightmass_z1 - self.r_standard) /
                             self.r_standard) * 1000  # [permille]
            # self.report(self.delta_z1, 'az1dC')
            # self.report(self.lightmass_z0, 'az1ML')

            # Change in storage - Pesticide Mass
            # self.conc_z1 = self.lightmass_z1 / (self.theta_z1 * self.z1)  # mg/mm
            ch_storage_z1_light = self.lightmass_z1 - self.lightmass_z1_ini
            self.lightmass_z1_ini = self.lightmass_z1

            #######################
            # Cumulative counters
            # self.cum_latflux_ug_z1 += z1_light_latflux['net_mass_latflux']

        # Update state variables
        # Change in storage - Moisture
        self.theta_z1 = z1_moisture["theta_final"]
        ch_storage_z1_m3 = (self.theta_z1 * self.z1 * cellarea() / 1000) - \
                           (self.theta_z1_ini * self.z1 * cellarea() / 1000)
        self.theta_z1_ini = self.theta_z1

        #######################################################################################
        # Layer z = 2
        # Temperature
        temp_dict_z2 = getLayerTemp(self, 2, bio_cover, temp_bare_soil)
        self.temp_z2_fin = temp_dict_z2["temp_layer"]
        # Moisture
        # PERCOL_z2 = True
        # if self.currentTimeStep() >= 180:
        #     PERCOL_z2 = True
        z2_moisture = getLayerMoisture(self, 2,
                                       precip, theta_wp, CN2, crop_type,
                                       jd_sim, jd_dev, jd_mid, jd_end, len_dev_stage,
                                       root_depth, pot_evapor, pot_transpir, depletable_water,
                                       k_sat_z2, root_depth_z2,
                                       theta_fcap_z2, theta_sat_z2,
                                       percolate=percolation_z1, PERCOL_z2=self.PERCOL)

        percolation_z2 = z2_moisture["percolate"]
        lat_netflow_z2 = z2_moisture["lat_flow"]
        lat_outflow_z2 = z2_moisture["cell_lat_outflow"]
        lat_inflow_z2 = z2_moisture["upstream_lat_inflow"]
        overflow_height_z2 = z2_moisture["overflow_height"]
        etp_z2 = z2_moisture["ETP"]

        if self.PEST:
            #########################
            # Mass Transfer, z2
            # Mass volatilized = not relevant @z2!
            # Mass runoff = not relevant @z2!
            # Mass & delta leached (Deep Percolation - DP, z2)
            if self.TRANSPORT:
                self.lightmass_z2 += z1_light_leached
                self.heavymass_z2 += z1_heavy_leached

                z2_light_leached = getLeachedMass(self, 2, theta_sat_z2, theta_fcap_z2,
                                                  percolation_z2,
                                                  z2_moisture["theta_after_percolate"],
                                                  self.lightmass_z2,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)
                z2_heavy_leached = getLeachedMass(self, 1, theta_sat_z2, theta_fcap_z2,
                                                  percolation_z2,
                                                  z2_moisture["theta_after_percolate"],
                                                  self.heavymass_z2,
                                                  sorption_model="linear",
                                                  leach_model="mcgrath", gas=True)

                self.lightmass_z2 -= z2_light_leached
                self.heavymass_z2 -= z2_heavy_leached

                # Mass & delta latflux (LF), z2
                z2_light_latflux = getLatMassFlux(self, 2, theta_sat_z2, theta_fcap_z2,
                                                  self.lightmass_z2,
                                                  sorption_model="linear", gas=True)
                z2_heavy_latflux = getLatMassFlux(self, 2, theta_sat_z2, theta_fcap_z2,
                                                  self.heavymass_z2,
                                                  sorption_model="linear", gas=True)

                self.lightmass_z2 += z2_light_latflux['net_mass_latflux']
                self.heavymass_z2 += z2_heavy_latflux['net_mass_latflux']

            # Degradation
            if self.DEGRADE:
                z2_deg_light = getMassDegradation(self, 2,
                                                  theta_sat_z2, theta_fcap_z2, theta_wp,
                                                  self.lightmass_z2, frac="L", sor_deg_factor=1)
                z2_deg_heavy = getMassDegradation(self, 2,
                                                  theta_sat_z2, theta_fcap_z2, theta_wp,
                                                  self.heavymass_z2, frac="H", sor_deg_factor=1)

                self.lightmass_z2 = z2_deg_light["mass_tot_new"]
                z2_light_deg = z2_deg_light.get("mass_deg_aq") + z2_deg_light.get("mass_deg_ads")
                self.heavymass_z2 = z2_deg_heavy["mass_tot_new"]

        # Update state variable
        self.theta_z2 = z2_moisture["theta_final"]

        # Baseflow
        self.baseflow = (self.theta_z2 * self.z2 * self.gw_factor) / self.k_g  # [mm/d]

        if self.PEST:
            if self.TRANSPORT:
                conc_aq_light_z2 = getConcAq(self, 2, self.theta_z2, self.lightmass_z2,
                                             sorption_model="linear", gas=True)  # ug/L
                conc_aq_heavy_z2 = getConcAq(self, 2, self.theta_z2, self.heavymass_z2,
                                             sorption_model="linear", gas=True)  # ug/L

                baseflow_light = self.baseflow * conc_aq_light_z2 * cellarea()  # ug
                baseflow_heavy = self.baseflow * conc_aq_heavy_z2 * cellarea()  # ug

                self.lightmass_z2 -= baseflow_light
                self.heavymass_z2 -= baseflow_heavy

            self.delta_z2 = ((self.heavymass_z2 / self.lightmass_z2 - self.r_standard) /
                             self.r_standard) * 1000  # [permille]
            # self.report(self.delta_z2, 'az2dC')
            # self.report(self.lightmass_z2, 'az2ML')

        # State Update
        self.theta_z2 -= self.baseflow / self.z2

        # Storage hydro
        ch_storage_z2_m3 = (self.theta_z2 * self.z2 * cellarea() / 1000) - \
                           (self.theta_z2_ini * self.z2 * cellarea() / 1000)
        self.theta_z2_ini = self.theta_z2

        # Storage pest
        if self.PEST:
            # Change in storage - Pesticide Mass
            # self.conc_z2 = self.pestmass_z2 / (self.theta_z2 * self.z2)  # mg/mm
            ch_storage_z2_light = self.lightmass_z2 - self.lightmass_z2_ini
            self.lightmass_z2_ini = self.lightmass_z2

        ###########################################################################
        # FINAL MASS BALANCE #
        ######################
        # 'Sample' the time-series associated to each component (e.g. runoff) at the outlet or due to accuflux()

        ######################
        # Water Balance
        ######################
        # Precipitation total
        rain_m3 = self.mask * precip * cellarea() / 1000  # m3
        tot_rain_m3 = areatotal(rain_m3, self.is_catchment)
        self.tot_rain_m3_tss.sample(tot_rain_m3)

        # Discharge due to runoff at the outlet
        runoff_m3 = runoff_z0 * cellarea() / 1000  # m3
        accu_runoff_m3 = accuflux(self.ldd_subs, runoff_m3)
        out_runoff_m3 = areatotal(accu_runoff_m3, self.outlet_multi)
        # self.report(out_runoff_m3, 'roffm3')
        self.out_runoff_m3_tss.sample(out_runoff_m3)  # save to outlet
        # self.obs_cum_runoff_m3_tss.sample(out_runoff_m3)  # save to sample locations

        # Discharge due to z1 drainage
        cell_drain_z1_m3 = drain_outflow_z1 * cellarea() / 1000  # m3
        accu_o_drain_z1_m3 = accuflux(self.ldd_subs, cell_drain_z1_m3)  # m3
        o_drain_z1_m3 = areatotal(accu_o_drain_z1_m3, self.outlet_multi)
        self.out_accu_o_drain_m3_tss.sample(o_drain_z1_m3)

        # Discharge due to lateral flow

        # Overflow
        overflow_cellvol_z0 = overflow_height_z0 * cellarea() / 1000  # m3
        overflow_cellvol_z1 = overflow_height_z1 * cellarea() / 1000  # m3
        overflow_cellvol_z2 = overflow_height_z2 * cellarea() / 1000  # m3
        column_overflow = overflow_cellvol_z2 + overflow_cellvol_z1 + overflow_cellvol_z0
        accu_of_latflow_m3 = accuflux(self.ldd_subs, column_overflow)
        of_latflow_m3 = areatotal(accu_of_latflow_m3, self.outlet_multi)
        self.sat_accu_overflow_m3_tss.sample(of_latflow_m3)

        # In/Outflow (Top Layer only for every cell!!)
        in_latflow_z0_m3 = lat_inflow_z0 * cellarea() / 1000  # m3
        # self.report(in_latflow_z0_m3, "inLF")
        in_latflow_z1_m3 = lat_inflow_z1 * cellarea() / 1000
        in_latflow_z2_m3 = lat_inflow_z2 * cellarea() / 1000
        cell_lat_inflow_m3 = in_latflow_z0_m3 + in_latflow_z1_m3 + in_latflow_z2_m3
        # self.report(cell_lat_inflow_m3, 'cellIn')
        outlet_cell_inflow_m3 = areatotal(cell_lat_inflow_m3, self.outlet_multi)  # Only outlet cells
        self.out_cell_i_latflow_m3_tss.sample(outlet_cell_inflow_m3)
        catch_lat_inflow_m3 = areatotal(cell_lat_inflow_m3, self.is_catchment)  # Sum catchment
        self.out_accu_i_latflow_m3_tss.sample(catch_lat_inflow_m3)

        out_latflow_z0_m3 = lat_outflow_z0 * cellarea() / 1000  # m3
        out_latflow_z1_m3 = lat_outflow_z1 * cellarea() / 1000
        out_latflow_z2_m3 = lat_outflow_z2 * cellarea() / 1000
        cell_lat_outflow_m3 = out_latflow_z0_m3 + out_latflow_z1_m3 + out_latflow_z2_m3
        # self.report(cell_lat_outflow_m3, 'cellOut')
        outlet_cell_outflow_m3 = areatotal(cell_lat_outflow_m3, self.outlet_multi)  # Only outlet cells
        self.out_cell_o_latflow_m3_tss.sample(outlet_cell_outflow_m3)
        catch_lat_outflow_m3 = areatotal(cell_lat_outflow_m3, self.is_catchment)  # Sum catchment
        self.out_accu_o_latflow_m3_tss.sample(catch_lat_outflow_m3)

        # Netflow
        catch_n_latflow_m3 = catch_lat_inflow_m3 - catch_lat_outflow_m3
        self.out_accu_n_latflow_m3_tss.sample(catch_n_latflow_m3)

        # Percolation (only interested in the bottom-most layer, where mass leaves the model)
        percol_z0_m3 = percolation_z0 * cellarea() / 1000  # m3
        percol_z1_m3 = percolation_z1 * cellarea() / 1000  # m3
        percol_z2_m3 = percolation_z2 * cellarea() / 1000  # m3
        # percol_m3 = percol_z2_m3 + percol_z1_m3 + percol_z0_m3
        accu_out_percol_m3 = accuflux(self.ldd_subs, percol_z2_m3)
        # accu_out_percol_m3 = accuflux(self.ldd_subs, percol_m3)
        out_percol_m3 = areatotal(accu_out_percol_m3, self.outlet_multi)
        self.out_percol_z2_m3_tss.sample(out_percol_m3)

        # Evapotranspiration
        etp_z0_m3 = etp_z0 * cellarea() / 1000  # m3
        etp_z1_m3 = etp_z1 * cellarea() / 1000
        etp_z2_m3 = etp_z2 * cellarea() / 1000
        etp_m3 = etp_z0_m3 + etp_z1_m3 + etp_z2_m3
        accu_out_etp_m3 = accuflux(self.ldd_subs, etp_m3)
        out_etp_m3 = areatotal(accu_out_etp_m3, self.outlet_multi)
        self.out_etp_m3_tss.sample(out_etp_m3)

        # Baseflow discharge (z2)
        accu_baseflow = accuflux(self.ldd_subs, self.baseflow * cellarea() / 1000)  # m3
        out_baseflow_m3 = areatotal(accu_baseflow, self.outlet_multi)
        self.out_baseflow_m3_tss.sample(out_baseflow_m3)

        # Change in storage
        ch_storage_m3 = ch_storage_z0_m3 + ch_storage_z1_m3 + ch_storage_z2_m3
        accu_out_ch_storage_m3 = accuflux(self.ldd_subs, ch_storage_m3)
        out_ch_storage_m3 = areatotal(accu_out_ch_storage_m3, self.outlet_multi)
        self.out_ch_storage_m3_tss.sample(out_ch_storage_m3)

        # GLOBAL Water
        global_mb_water = (tot_rain_m3 -
                           out_etp_m3 - out_runoff_m3 -
                           out_percol_m3 + catch_n_latflow_m3 - of_latflow_m3 - o_drain_z1_m3 -
                           out_ch_storage_m3 - out_baseflow_m3)

        self.global_mb_water_tss.sample(global_mb_water)

        cell_vol_z0 = self.theta_z0 * self.z0 * cellarea() / 1000
        cell_vol_z1 = self.theta_z1 * self.z1 * cellarea() / 1000
        cell_vol_z2 = self.theta_z2 * self.z2 * cellarea() / 1000
        cell_vol_tot_m3 = cell_vol_z0 + cell_vol_z1 + cell_vol_z2
        vol_tot_m3 = accuflux(self.ldd_subs, cell_vol_tot_m3)
        multi_vol_tot_m3 = areatotal(vol_tot_m3, self.outlet_multi)
        self.storage_m3_tss.sample(multi_vol_tot_m3)

        # Runoff + accu_drainage + of_latflow_m3
        # in_vol_disch_m3 = out_runoff_m3 + i_latflow_m3
        # out_vol_disch_m3 = out_runoff_m3 + o_latflow_m3
        net_vol_disch_m3 = out_runoff_m3 + o_drain_z1_m3 + of_latflow_m3 + outlet_cell_outflow_m3 + out_baseflow_m3

        # self.i_Q_m3_tss.sample(in_vol_disch_m3)
        # self.o_Q_m3_tss.sample(out_vol_disch_m3)
        self.n_Q_m3_tss.sample(net_vol_disch_m3)

        if self.PEST:
            ######################
            # Pesticide Balance
            ######################

            # Applied ug on catchment
            accu_app_light = accuflux(self.ldd_subs, light_applied)  # ug
            # out_app_light = areatotal(accu_app_light, self.outlet_multi)  # ug
            out_app_light = areatotal(light_applied, self.is_catchment)  # ug
            self.out_app_L_tss.sample(out_app_light)

            # cum_appl_catch_light = accuflux(self.ldd_subs, self.cum_appl_ug)
            # Todo: Check if below same result
            # cum_appl_catch_mg = upstream(self.ldd_subs, self.cum_appl_ug)
            if self.TRANSPORT:
                # Degradation
                light_deg = z0_light_deg + z1_light_deg + z2_light_deg
                out_deg_light = areatotal(light_deg, self.is_catchment)

                # Volatilized
                out_volat_light = areatotal(light_volat, self.is_catchment)
                self.out_volat_L_tss.sample(out_volat_light)

                # Mass loss to run-off
                out_runoff_light = areatotal(light_runoff, self.is_catchment)
                self.out_runoff_L_tss.sample(out_runoff_light)
                # cum_out_runoff_light = accuflux(self.ldd_subs, self.cum_runoff_ug)

                # Mass loss to leaching  (zero as long as no percolation)
                # accu_leach_light_z0 = accuflux(self.ldd_subs, z0_light_leached)
                # accu_leach_light_z1 = accuflux(self.ldd_subs, z1_light_leached)
                # Accumulated upstream cells
                # self.report(accu_leach_light_z0, 'az0LCH')
                # self.report(accu_leach_light_z1, 'az1LCH')
                # Cell level
                # self.report(z0_light_leached, 'az0LCH')
                # self.report(z1_light_leached, 'az1LCH')
                # out_leach_light = areatotal(accu_leach_light, self.outlet_multi)
                # out_leach_light_z0 = areatotal(z0_light_leached, self.is_catchment)
                # out_leach_light_z1 = areatotal(z1_light_leached, self.is_catchment)
                out_leach_light = areatotal(z2_light_leached, self.is_catchment)
                # out_leach_light = out_leach_light_z0 + out_leach_light_z1 + out_leach_light_z2
                self.out_leach_L_tss.sample(out_leach_light)
                # cum_out_leach_light = accuflux(self.ldd_subs, self.cum_leached_ug_z2)

                # Drained mass (layer z1)
                out_drain_light = areatotal(z1_light_drain, self.is_catchment)
                self.out_drain_L_tss.sample(out_drain_light)
                # out_drain_heavy = areatotal(z1_heavy_drain, self.is_catchment)

                # Net loss to lateral flux
                mass_latflux_light = (z0_light_latflux["net_mass_latflux"] +
                                      z1_light_latflux["net_mass_latflux"] +
                                      z2_light_latflux["net_mass_latflux"])

                mass_loss_light_z0 = z0_light_latflux['mass_loss']
                mass_loss_light_z1 = z1_light_latflux['mass_loss']
                mass_loss_light_z2 = z2_light_latflux['mass_loss']

                mass_gain_light_z0 = z0_light_latflux['mass_gain']
                mass_gain_light_z1 = z1_light_latflux['mass_gain']
                mass_gain_light_z2 = z2_light_latflux['mass_gain']

                net_latflux_light = areatotal(mass_latflux_light, self.is_catchment)
                self.out_latflux_L_tss.sample(net_latflux_light)
                # tot_latflux_mg = self.cum_latflux_ug_z0 + self.cum_latflux_ug_z1 + self.cum_latflux_ug_z2
                # cum_out_latflux_light = accuflux(self.ldd_subs, tot_latflux_mg)

                # Baseflow flux
                out_baseflow_light = areatotal(baseflow_light, self.is_catchment)
                self.out_baseflow_L_tss.sample(out_baseflow_light)

                # Change in mass storage
                # ... per time step
                ch_storage_light = ch_storage_z0_light + ch_storage_z1_light + ch_storage_z2_light
                # self.report(ch_storage_z0_light, 'az0CH')
                # self.report(ch_storage_z1_light, 'az1CH')
                # self.report(ch_storage_z2_light, 'az2CH')
                # self.report(ch_storage_light, 'atotCH')

                # accu_ch_storage_light = accuflux(self.ldd_subs, ch_storage_light)
                out_ch_storage_light = areatotal(ch_storage_light, self.is_catchment)

                # self.report(accu_ch_storage_light, 'accuCH')
                # out_ch_storage_light = areatotal(accu_ch_storage_light, self.outlet_multi)
                # self.report(out_ch_storage_light, 'areaCH')
                self.out_chstorage_L_tss.sample(out_ch_storage_light)

                # out_ch_storage_light = accuflux(accu_ch_storage_light, self.outlet_multi)

                # cum_ch_storage_light = ch_storage_light - self.light_ini_storage_ug
                # cum_tot_ch_storage_light = accuflux(self.ldd_subs, cum_ch_storage_light)

                # TODO: still need to add outlet cells to MB
                light_mb_pest = (out_app_light - out_deg_light -
                                 out_volat_light - out_runoff_light -
                                 out_leach_light -
                                 out_drain_light + net_latflux_light -
                                 out_baseflow_light -
                                 out_ch_storage_light)
                # self.report(light_mb_pest, 'areaMB')
                self.global_mb_pest_tss.sample(light_mb_pest)

        self.jd_cum += self.jd_dt  # updating JDcum, currently dt = 1 day

        # Analysis (NASH Discharge)
        q_obs = timeinputscalar('q_obs_m3day.tss', nominal("outlet_true"))
        rest_obs = tot_rain_m3 / q_obs
        self.rest_obs_tss.sample(rest_obs)

        self.q_obs_tot += ifthenelse(q_obs >= 0, q_obs, 0)
        self.q_sim_tot += ifthenelse(q_obs >= 0, net_vol_disch_m3, 0)
        self.q_diff += ifthenelse(q_obs >= 0, (net_vol_disch_m3 - q_obs) ** 2, 0)
        self.q_var += ifthenelse(q_obs >= 0, (q_obs - 260.07) ** 2, 0)  # Mean discharge of data range = 260.07 m3/day
        nash_q = 1 - (self.q_diff / self.q_var)
        self.nash_q_tss.sample(nash_q)
        self.q_obs_tot_tss.sample(self.q_obs_tot)
        self.q_sim_tot_tss.sample(self.q_sim_tot)
        self.rain_cum_m3 += tot_rain_m3
        self.rain_obs_tot_tss.sample(self.rain_cum_m3)
        # aguila 1\res_nash_q_m3.tss 2\res_nash_q_m3.tss

        # Mass balance components
        self.tot_runoff_m3 += ifthenelse(q_obs >= 0, out_runoff_m3, 0)
        self.tot_runoff_m3_tss.sample(self.tot_runoff_m3)
        self.tot_perc_z2_m3 += ifthenelse(q_obs >= 0, out_percol_m3, 0)
        self.tot_perc_z2_m3_tss.sample(self.tot_perc_z2_m3)
        self.tot_etp_m3 += ifthenelse(q_obs >= 0, out_etp_m3, 0)
        self.tot_etp_m3_tss.sample(self.tot_etp_m3)
        # self.tot_baseflow_m3 += ifthenelse(q_obs >= 0, out_baseflow_m3, 0)
        # self.tot_baseflow_m3_tss.sample(self.tot_baseflow_m3)

        self.tot_drain_m3 += ifthenelse(q_obs >= 0, o_drain_z1_m3, 0)
        self.tot_accu_drain_m3_tss.sample(self.tot_drain_m3)

        # self.tot_ilf_m3 += ifthenelse(q_obs >= 0, outlet_lat_inflow_m3, 0)
        # self.tot_accu_i_latflow_m3_tss.sample(self.tot_ilf_m3)
        self.tot_olf_m3 += ifthenelse(q_obs >= 0, outlet_cell_outflow_m3, 0)
        self.tot_accu_o_latflow_m3_tss.sample(self.tot_olf_m3)
        # self.tot_nlf_m3 += ifthenelse(q_obs >= 0, n_latflow_m3, 0)
        # self.tot_accu_n_latflow_m3_tss.sample(self.tot_nlf_m3)
        self.tot_of_m3 += ifthenelse(q_obs >= 0, of_latflow_m3, 0)
        self.tot_accu_of_latflow_m3_tss.sample(self.tot_of_m3)

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
# aguila --scenarios='{1,2}' --multi=1x2  --timesteps=[175,179,1] aLEACH aLEACHz aLF aLFz
# aguila --scenarios='{1}'  --timesteps=[100,280,1] az0dC az1dC az2dC
# aguila --scenarios='{1}'  --timesteps=[1,280,1] aZ1LCH

# Time series
# aguila 1\res_nash_q_m3.tss 6\res_nash_q_m3.tss
# aguila 1\resM_norCONC.tss 1\resM_valCONC.tss 1\resM_souCONC.tss

nrOfSamples = 2  # Samples are each a MonteCarlo realization
firstTimeStep = 1
nTimeSteps = 280
myAlteck16 = BeachModel("clone_nom.map")  # an instance of the model, which inherits from class: DynamicModel
dynamicModel = DynamicFramework(myAlteck16, lastTimeStep=nTimeSteps,
                                firstTimestep=firstTimeStep)  # an instance of the Dynamic Framework
mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)

t0 = datetime.now()
# dynamicModel.run()
mcModel.run()
t1 = datetime.now()

duration = t1 - t0
print("Total minutes: ", duration.total_seconds() / 60)
print("Minutes/Yr: ", (duration.total_seconds() / 60) / (nTimeSteps - firstTimeStep) * 365)
