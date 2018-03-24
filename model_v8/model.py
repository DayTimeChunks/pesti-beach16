# -*- coding: utf-8 -*-

# from time import *
import time
from datetime import datetime
from hydro import *
from pesti import *

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
        self.zero_map = self.dem - self.dem  # Zero map to generate scalar maps
        self.mask = self.dem / self.dem
        self.datum_depth = (self.dem - mapminimum(self.dem)) * scalar(10 ** 3)  # mm

        # self.ldd_surf = lddcreate(self.dem_route, 1e31, 1e31, 1e31, 1e31)  # To route runoff
        self.ldd_subs = lddcreate(self.dem, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI

        self.outlet = self.readmap("outlet_true")

        # TODO: temporary fix to landuse.map with holes!!!
        self.landuse = self.readmap("fields_cover")

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
        # Pesticide
        self.global_mb_pest_tss = TimeoutputTimeseries("res_global_mb_pest", self, nominal("outlet_true"),
                                                       noHeader=False)
        self.out_delta_tss = TimeoutputTimeseries("out_delta", self, nominal("outlet_true"), noHeader=False)
        self.out_mass_tss = TimeoutputTimeseries("out_mass", self, nominal("outlet_true"), noHeader=False)

        # Rain
        self.tot_rain_m3_tss = TimeoutputTimeseries("res_accuRain_m3", self, nominal("outlet_true"), noHeader=False)
        # Runoff
        self.out_runoff_m3_tss = TimeoutputTimeseries("res_accuRunoff_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_runoff_m3_tss = TimeoutputTimeseries("res_totRunoff_m3", self, nominal("outlet_true"), noHeader=False)
        # Percolation
        # self.out_percol_m3_tss = TimeoutputTimeseries("res_outPercol_m3", self, nominal("outlet_true"), noHeader=False)
        self.out_percol_z2_m3_tss = TimeoutputTimeseries("res_accuPercol_z2_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_perc_z2_m3_tss = TimeoutputTimeseries("res_totPercol_z2_m3", self, nominal("outlet_true"), noHeader=False)

        # ETP
        self.out_etp_m3_tss = TimeoutputTimeseries("res_accuEtp_m3", self, nominal("outlet_true"), noHeader=False)
        self.tot_etp_m3_tss = TimeoutputTimeseries("res_totEtp_m3", self, nominal("outlet_true"), noHeader=False)

        # LF Drainage
        self.out_accu_o_drain_m3_tss = TimeoutputTimeseries("res_o_accuDrain_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_drain_m3_tss = TimeoutputTimeseries("res_o_totDrain_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        # LF options
        self.sat_accu_overflow_m3_tss = TimeoutputTimeseries("res_of_accuLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)

        self.tot_accu_of_latflow_m3_tss = TimeoutputTimeseries("res_of_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)

        self.out_accu_i_latflow_m3_tss = TimeoutputTimeseries("res_i_accuLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_i_latflow_m3_tss = TimeoutputTimeseries("res_i_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.out_accu_o_latflow_m3_tss = TimeoutputTimeseries("res_o_accuLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_o_latflow_m3_tss = TimeoutputTimeseries("res_o_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.out_accu_n_latflow_m3_tss = TimeoutputTimeseries("res_n_accuLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)
        self.tot_accu_n_latflow_m3_tss = TimeoutputTimeseries("res_n_totLatflow_m3", self, nominal("outlet_true"),
                                                              noHeader=False)

        self.out_ch_storage_m3_tss = TimeoutputTimeseries("res_accuChStorage_m3", self, nominal("outlet_true"),
                                                          noHeader=False)
        self.global_mb_water_tss = TimeoutputTimeseries("res_global_waterMB", self, nominal("outlet_true"),
                                                        noHeader=False)
        self.storage_m3_tss = TimeoutputTimeseries("res_accuStorage_m3", self, nominal("outlet_true"), noHeader=False)

        # Analysis
        # This is 'q' as time series.
        self.i_Q_m3_tss = TimeoutputTimeseries("res_i_accuVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.o_Q_m3_tss = TimeoutputTimeseries("res_o_accuVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.n_Q_m3_tss = TimeoutputTimeseries("res_n_accuVol_m3", self, nominal("outlet_true"), noHeader=False)
        self.q_obs_tot_tss = TimeoutputTimeseries("res_q_obs_tot_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # Equivalent to net_Q
        self.rain_obs_tot_tss = TimeoutputTimeseries("res_rain_obs_tot_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # Equivalent to net_Q
        self.rest_obs_tss = TimeoutputTimeseries("res_restit_obs_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # = rain/q_obs
        self.q_sim_tot_tss = TimeoutputTimeseries("res_q_sim_tot_m3", self, nominal("outlet_true"),
                                                  noHeader=False)  # This is 'Nash_q' as time series.
        self.nash_q_tss = TimeoutputTimeseries("res_nash_q_m3", self, nominal("outlet_true"),
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
        self.s1 = scalar(self.ini_param.get("s1"))# mm/day
        self.s2 = scalar(self.ini_param.get("s2"))  # coefficient to calibrate Ksat2



        # First-sensitivity:
        """ Physical parameters """
        self.c1 = scalar(self.ini_param.get("c1"))  # subsurface flow coefficient
        self.c2 = scalar(self.ini_param.get("c2"))  # not used (second layer)
        self.drain_coef = scalar(self.ini_param.get("drain_coef"))  # drainage coefficient
        self.k = scalar(self.ini_param.get("k"))  # coefficient of declining LAI in end stage

        # First run will return state = 1
        m_state = get_state(state)
        self.PERCOL = False  # z2 deep percolation (DP)
        self.ADLF = True
        self.c_adr = 0.25
        z2_factor = 0.9
        if m_state == 1:
            pass  # c1, c2 = 0.15, d1 = 0.8
        elif m_state == 2:
            z2_factor = 0.8
        elif m_state == 3:
            z2_factor = 0.7


        self.z0 = self.zero_map + 10
        self.z1 = self.zero_map + 300  # mm
        self.z2 = (self.datum_depth + 610 - self.z0 - self.z1)*z2_factor  # mm (300mm at outlet) * z2_factor
        self.tot_depth = self.z0 + self.z1 + self.z2

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
        # https://www.gsi-net.com
        self.k_h = max(scalar(self.ini_param.get("k_h")), scalar(3.13e-08))  # Henry's constant @ 20 C (dimensionless, Metolachlor)

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

        """
        Isotopes
        """
        self.r_standard = scalar(self.ini_param.get("r_standard"))  # VPDB
        self.alpha_iso = scalar(self.ini_param.get("alpha_iso"))  # 1 is no fractionation

        """
        Hydro Maps
        """
        # Initial moisture (Final from model v1, Sept 30, 2016)
        self.theta_z0 = readmap('ini_theta_z0')  # map of initial soil moisture in top layer (-)
        self.theta_z1 = readmap('ini_theta_z1')
        self.theta_z2 = max(readmap("ini_theta_z2"), scalar(0.1))

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
            # in ug/m2 = conc. (ug/g soil) * density (g/cm3) * (10^6 cm3 / m3) * (1 m/10^3 mm) * depth_layer (mm)
            self.smback_z0 = (self.zero_map + 0.06) * self.p_b * scalar(
                10 ** 6 / 10 ** 3) * self.z0  # Based on detailed soils
            self.smback_z1 = (self.zero_map + 0.03) * self.p_b * scalar(
                10 ** 6 / 10 ** 3) * self.z1  # Based on detailed soils
            self.smback_z2 = (self.zero_map + 0.00001) * self.p_b * scalar(10 ** 6 / 10 ** 3) * self.z2  # Assumed

            # Carbon Delta (Background)
            # Assumed theoretical max @99% deg Streitwieser Semiclassical Limits
            self.delta_z0 = self.zero_map - 23.7
            self.delta_z0_ini = self.delta_z0
            self.delta_z1 = self.zero_map - 23.7
            self.delta_z1_ini = self.delta_z1
            self.delta_z2 = self.zero_map - 23.7
            self.delta_z2_ini = self.delta_z2

            # Applications Mass
            # Product concentration (active ing.)
            double = scalar(2.0)  # ~ Dosage for corn when growing beet
            d_gold = scalar(915) * 10 ** 6  # ug/L S-met
            m_gold = scalar(960) * 10 ** 6  # ug/L

            # Dosages # L/Ha * 1Ha/1000m2 = L/m2
            d_beet = None
            d_corn = scalar(2.1) * 1 / 10 ** 4  # L/Ha * 1 Ha / 10000 m2
            m_beet = scalar(0.6) * 1 / 10 ** 4 * double
            m_corn = scalar(2.0) * 1 / 10 ** 4
            m_beet_Friess = scalar(0.6) * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)
            m_beet_Mathis = scalar(0.6) * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)

            # Assign dosages based on Farmer-Crop combinations [ug/m2]
            fa_cr = readmap("farmer_crop")  # Contains codes to assign appropriate dosage
            app_conc = (  # [ug/m2]
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
            # Pesticide applied (ug/m2) on Julian day 177 (March 25, 2016).
            # March 26th, Friess and Mathis
            self.app1 = ifthenelse(fa_cr == 1111, 1 * app_conc,
                                   # 1111 (Friess, Beet), 1112 (Friess-Corn),
                                   ifthenelse(fa_cr == 1112, 1 * app_conc,
                                              ifthenelse(fa_cr == 1711, 1 * app_conc,  # 1711 (Mathis-Beet)
                                                         0 * app_conc)))
            # Pesticide applied (ug/m2) on Julian day 197 (April 14, 2016).
            # April 14, Kopp and Burger
            self.app2 = ifthenelse(fa_cr == 1511, 1 * app_conc,  # 1511 (Burger-Beet)
                                   ifthenelse(fa_cr == 1611, 1 * app_conc,  # 1611 (Kopp-Beet),
                                              0 * app_conc))

            # Pesticide applied (ug/m2) on Julian day 238 (May 25, 2016).
            # May 25, Schmidt and Speich, and (out of transect): Friess and Mahler
            # Note: Speich and Friess could be 1 week later.
            self.app3 = ifthenelse(fa_cr == 1112, 1 * app_conc,  # 1112 (Friess-Corn)
                                   ifthenelse(fa_cr == 1212, 1 * app_conc,  # 1212 (Speich-Corn),
                                              ifthenelse(fa_cr == 1412, 1 * app_conc,  # 1412 (Schmitt-Corn),
                                                         ifthenelse(fa_cr == 1312, 1 * app_conc,  # 1312 (Mahler-Corn)
                                                                    0 * app_conc))))

            # Applications delta
            # Use map algebra to produce a initial signature map,
            # ATT: Need to do mass balance on addition of new layer.
            # where app1 > 0, else background sig. (plots with no new mass will be 0)
            self.app1delta = ifthenelse(self.app1 > 0, scalar(-32.3), scalar(-23.7))
            self.app2delta = ifthenelse(self.app2 > 0, scalar(-32.3), scalar(-23.7))
            self.app3delta = ifthenelse(self.app3 > 0, scalar(-32.3), scalar(-23.7))

            # Convert mg/m2 -> mg
            self.pestmass_z0 = self.smback_z0 * cellarea()  # mg
            self.pestmass_z0_ini = self.pestmass_z0
            self.pestmass_z1 = self.smback_z1 * cellarea()  # mg
            self.pestmass_z1_ini = self.pestmass_z1
            self.pestmass_z2 = self.smback_z2 * cellarea()  # mg
            self.pestmass_z2_ini = self.pestmass_z2

            # Cumulative maps
            self.pest_ini_storage_mg = self.pestmass_z0_ini + self.pestmass_z1_ini + self.pestmass_z2_ini
            self.cum_appl_mg = self.zero_map
            self.cum_runoff_mg = self.zero_map
            self.cum_leached_mg_z2 = self.zero_map  # Only bottom-most layer needed
            self.cum_latflux_mg_z0 = self.zero_map
            self.cum_latflux_mg_z1 = self.zero_map
            self.cum_latflux_mg_z2 = self.zero_map

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
        self.tot_perc_z2_m3 = 0
        self.tot_runoff_m3 = 0


        # Stochastic / test parameters
        print("state:", m_state)
        # self.report(self.c1, 'c1')
        # self.report(self.c2, 'c2')
        # self.report(self.drain_coef, 'd1')
        # self.report(self.z0, 'z0')
        # self.report(self.s1, 's1')
        self.report(self.pestmass_z0, 'z0mini')
        self.report(self.delta_z0_ini, 'z0dCini')
        self.report(self.app1, 'z0app1')
        self.report(self.app1delta, 'z0app1dC')


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
        temp_air = timeinputscalar('airTemp.tss', nominal('clone_nom'))
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
        self.theta_z0_ini = self.theta_z0

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

        # self.obs_runoff_m3_tss.sample(runoff_z0*4/1000)

        if self.PEST:
            #########################
            # Mass Transfer, z0
            mass_loss_dt_z0 = scalar(0)
            mass_gain_dt_z0 = scalar(0)

            # Background
            self.pestmass_z0_ini = self.pestmass_z0  # mg

            # Applications
            mass_applied = self.zero_map
            # z0_mass_volatilized = self.zero_map
            if self.currentTimeStep() in self.app_days:
                mass_applied = ifthenelse(self.currentTimeStep() == self.app_days[0],
                                          self.app1 * cellarea(),
                                          ifthenelse(self.currentTimeStep() == self.app_days[1], self.app2 * cellarea(),
                                                     ifthenelse(self.currentTimeStep() == self.app_days[2], self.app3 * cellarea(),
                                                                scalar(0))))  # [mg]
                self.cum_appl_mg += mass_applied
                self.pestmass_z0 += mass_applied  # mg
                mass_loss_dt_z0 += scalar(0)
                mass_gain_dt_z0 += mass_applied

                # Isotopes change due to applications
                self.delta_z0_ini = self.delta_z0
                delta_applied = ifthenelse(self.currentTimeStep() == 177, self.app1delta,
                                           ifthenelse(self.currentTimeStep() == 197, self.app2delta,
                                                      ifthenelse(self.currentTimeStep() == 238, self.app3delta,
                                                                 scalar(0))))  # [delta permille]
                # isotope mass balance (due to application only)
                self.delta_z0 = scalar(1) / self.pestmass_z0 * (
                    self.delta_z0_ini * self.pestmass_z0_ini + delta_applied * mass_applied)

                self.report(self.delta_z0, 'z0dCAp')
                self.report(self.pestmass_z0, 'z0mssAp')

                # Mass & delta volatilized
                mass_before_transport = self.pestmass_z0
                z0_mass_volatilized = getVolatileMass(self, self.app_days,  # Application days
                                                      temp_air, theta_sat_z0z1,
                                                      rel_diff_model='option-1', sorption_model="linear",
                                                      gas=True)
                self.pestmass_z0 -= z0_mass_volatilized["mass_loss"]
                self.delta_z0 = update_layer_delta(self, 0, "volat", z0_mass_volatilized, mass_before_transport)

                self.report(self.delta_z0, 'z0dCVl')
                self.report(self.pestmass_z0, 'z0mssVl')

            # TODO: Mass loss due to plant uptake!

            # Mass & delta run-off (RO)
            mass_before_transport = self.pestmass_z0
            # transfer_model = "d-mlm"
            z0_mass_runoff = getRunOffMass(self, theta_sat_z0z1,
                                           precip, runoff_z0,
                                           transfer_model="d-mlm", sorption_model="linear")
            self.pestmass_z0 -= z0_mass_runoff["mass_runoff"]  # mg
            self.delta_z0 = update_layer_delta(self, 0, "runoff", z0_mass_runoff, mass_before_transport)
            #self.report(self.pestmass_z0, 'z0ro_M')
            #self.report(self.delta_z0, 'z0ro_dC')

            # Mass & delta leached (Deep Percolation - DP)
            z0_mass_leached = getLeachedMass(self, 0, theta_sat_z0z1,
                                             precip,
                                             tot_percolation_z0,
                                             z0_moisture["theta_after_percolate"],
                                             sorption_model="linear",
                                             leach_model="mcgrath")
            mass_before_transport = self.pestmass_z0
            self.pestmass_z0 -= z0_mass_leached["mass_leached"]  # mg
            self.delta_z0 = update_layer_delta(self, 0, "leach", z0_mass_leached, mass_before_transport)
            #self.report(self.pestmass_z0, 'z0lch_M')
            #self.report(self.delta_z0, 'z0lch_dC')

            # Mass & delta latflux (LF)
            mass_before_transport = self.pestmass_z0
            z0_mass_latflux = getLatMassDeltaFlux(self, 0, theta_sat_z0z1, theta_fcap_z0z1, mass_before_transport)

            self.pestmass_z0 += z0_mass_latflux["net_mass_latflux"]  # mg
            self.delta_z0 = z0_mass_latflux["delta_layer"]
            f1 = z0_mass_latflux["f1"]
            f2 = z0_mass_latflux["f2"]
            f3 = z0_mass_latflux["f3"]
            self.report(f1, 'f1')
            self.report(f2, 'f2')
            self.report(f3, 'f3')
            self.report(self.pestmass_z0, 'z0lf_M')
            self.report(self.delta_z0, 'z0lf_dC')

            # Degradation
            mass_before_degradation = self.pestmass_z0
            deg_z0_dict = degrade(self, 0,
                                  theta_sat_z0z1, theta_sat_z2,
                                  theta_fcap_z0z1, theta_wp,
                                  sor_deg_factor=1)
            self.pestmass_z0 = deg_z0_dict["mass_light_fin"] + deg_z0_dict["mass_heavy_fin"]
            self.delta_z0 = (deg_z0_dict["mass_heavy_fin"] / deg_z0_dict[ "mass_light_fin"] - self.r_standard) / self.r_standard
            self.report(self.pestmass_z0, 'z0deg_M')
            self.report(self.delta_z0, 'z0deg_dC')
            self.report(self.r_standard, 'rstd')

            # # Testing
            # if self.jd_cum in {2, 5, 8}:
            #     t0 = time.time()
            #     self.report(self.delta_z1, "dz0_" + str(self.jd_cum))
            #     t1 = time.time()
            #     print("Total:", t1 - t0)

            # Change in mass storage after degradation - Pesticide Mass
            ch_storage_z0_mg = self.pestmass_z0 - self.pestmass_z0_ini
            self.pestmass_z0_ini = self.pestmass_z0

            # Cumulative
            self.cum_runoff_mg += z0_mass_runoff["mass_runoff"]
            self.cum_latflux_mg_z0 += z0_mass_latflux["net_mass_latflux"]

        # Update state variables
        # Change in storage - Moisture
        self.theta_z0 = z0_moisture["theta_final"]
        ch_storage_z0_m3 = (self.theta_z0 * self.z0 * cellarea() / 1000) - \
                           (self.theta_z0_ini * self.z0 * cellarea() / 1000)
        self.theta_z0_ini = self.theta_z0

        # self.theta_z0tss.sample(self.theta_z0)
        # self.water_balance_z0tss.sample(z0_moisture["balance"])

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
            self.pestmass_z1 += z0_mass_leached["mass_leached"]
            z1_mass_leached = getLeachedMass(self, 1, theta_sat_z0z1,
                                             precip,
                                             percolation_z1,
                                             z1_moisture["theta_after_percolate"],
                                             sorption_model="linear")
            mass_before_transport = self.pestmass_z1
            self.pestmass_z1 -= z1_mass_leached["mass_leached"]
            self.delta_z1 = update_layer_delta(self, 1, "leach", z1_mass_leached, mass_before_transport)

            # Mass & delta latflux (LF), z1
            mass_before_transport = self.pestmass_z1
            z1_mass_latflux = getLatMassDeltaFlux(self, 1, theta_sat_z0z1, theta_fcap_z0z1, mass_before_transport)
            self.pestmass_z1 += z1_mass_latflux["net_mass_latflux"]
            self.delta_z1 = z1_mass_latflux["delta_layer"]

            # Degradation
            mass_before_degradation = self.pestmass_z1
            deg_z1_dict = degrade(self, 1,
                                  theta_sat_z0z1, theta_sat_z2,
                                  theta_fcap_z0z1, theta_wp,
                                  sor_deg_factor=1)
            self.pestmass_z1 = deg_z1_dict["mass_light_fin"] + deg_z1_dict["mass_heavy_fin"]
            self.delta_z1 = (deg_z1_dict["mass_heavy_fin"] / deg_z1_dict[
                "mass_light_fin"] - self.r_standard) / self.r_standard

            # Change in storage - Pesticide Mass
            self.conc_z1 = self.pestmass_z1 / (self.theta_z1 * self.z1)  # mg/mm
            ch_storage_z1_mg = self.pestmass_z1 - self.pestmass_z1_ini
            self.pestmass_z1_ini = self.pestmass_z1

            #######################
            # Cumulative counters
            self.cum_latflux_mg_z1 += z1_mass_latflux["net_mass_latflux"]

        # Update state variables
        # Change in storage - Moisture
        self.theta_z1 = z1_moisture["theta_final"]
        ch_storage_z1_m3 = (self.theta_z1 * self.z1 * cellarea() / 1000) - \
                           (self.theta_z1_ini * self.z1 * cellarea() / 1000)
        self.theta_z1_ini = self.theta_z1

        # SAVE
        # self.theta_z1tss.sample(self.theta_z1)
        # self.water_balance_z1tss.sample(z1_moisture["balance"])

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
            self.pestmass_z2 += z1_mass_leached["mass_leached"]
            z2_mass_leached = getLeachedMass(self, 2, theta_sat_z2,
                                             precip,
                                             percolation_z2,
                                             z2_moisture["theta_after_percolate"],
                                             sorption_model="linear")
            mass_before_transport = self.pestmass_z2
            self.pestmass_z2 -= z2_mass_leached["mass_leached"]
            self.delta_z2 = update_layer_delta(self, 2, "leach", z2_mass_leached, mass_before_transport)

            # Mass & delta latflux (LF), z2
            mass_before_transport = self.pestmass_z2
            z2_mass_latflux = getLatMassDeltaFlux(self, 0, theta_sat_z0z1, theta_fcap_z0z1, mass_before_transport)
            self.pestmass_z2 += z2_mass_latflux["net_mass_latflux"]
            self.delta_z2 = z2_mass_latflux["delta_layer"]

            # Degradation
            mass_before_degradation = self.pestmass_z2
            deg_z2_dict = degrade(self, 2,
                                  theta_sat_z0z1, theta_sat_z2,
                                  theta_fcap_z0z1, theta_wp,
                                  sor_deg_factor=1)
            self.pestmass_z2 = deg_z2_dict["mass_light_fin"] + deg_z2_dict["mass_heavy_fin"]
            self.delta_z2 = (deg_z2_dict["mass_heavy_fin"] / deg_z2_dict[
                "mass_light_fin"] - self.r_standard) / self.r_standard

            # self.theta_z2tss.sample(self.theta_z2)
            # self.water_balance_z2tss.sample(z2_moisture["balance"])

            # Change in storage - Pesticide Mass
            self.conc_z2 = self.pestmass_z2 / (self.theta_z2 * self.z2)  # mg/mm
            ch_storage_z2_mg = self.pestmass_z2 - self.pestmass_z2_ini
            self.pestmass_z2_ini = self.pestmass_z2

            #################
            # Cumulative counters
            self.cum_leached_mg_z2 += z2_mass_leached["mass_leached"]
            self.cum_latflux_mg_z2 += z2_mass_latflux["net_mass_latflux"]

        # Update state variables
        # Change in storage - Moisture
        self.theta_z2 = z2_moisture["theta_final"]
        ch_storage_z2_m3 = (self.theta_z2 * self.z2 * cellarea() / 1000) - \
                           (self.theta_z2_ini * self.z2 * cellarea() / 1000)
        self.theta_z2_ini = self.theta_z2

        ###########################################################################
        # FINAL MASS BALANCE #
        ######################
        # 'Sample' the time-series associated to each component (e.g. runoff) at the outlet or due to accuflux()

        ######################
        # Water Balance
        ######################
        # Precipitation total
        rain_m3 = precip * cellarea() / 1000  # m3
        tot_rain_m3 = accuflux(self.ldd_subs, rain_m3)
        self.tot_rain_m3_tss.sample(tot_rain_m3)

        # Discharge due to runoff at the outlet
        runoff_m3 = runoff_z0 * cellarea() / 1000  # m3
        out_runoff_m3 = accuflux(self.ldd_subs, runoff_m3)
        self.out_runoff_m3_tss.sample(out_runoff_m3)  # save to outlet
        # self.obs_cum_runoff_m3_tss.sample(out_runoff_m3)  # save to sample locations

        # Discharge due to z1 drainage
        cell_drain_z1_m3 = drain_outflow_z1 * cellarea() / 1000  # m3
        o_drain_z1_m3 = accuflux(self.ldd_subs, cell_drain_z1_m3)  # m3
        self.out_accu_o_drain_m3_tss.sample(o_drain_z1_m3)

        # Discharge due to lateral flow

        # Overflow
        overflow_cellvol_z0 = overflow_height_z0 * cellarea() / 1000  # m3
        overflow_cellvol_z1 = overflow_height_z1 * cellarea() / 1000  # m3
        overflow_cellvol_z2 = overflow_height_z2 * cellarea() / 1000  # m3
        column_overflow = overflow_cellvol_z2 + overflow_cellvol_z1 + overflow_cellvol_z0
        of_latflow_m3 = accuflux(self.ldd_subs, column_overflow)
        self.sat_accu_overflow_m3_tss.sample(of_latflow_m3)

        # In/Outflow
        in_latflow_z0_m3 = lat_inflow_z0 * cellarea() / 1000  # m3
        in_latflow_z1_m3 = lat_inflow_z1 * cellarea() / 1000
        in_latflow_z2_m3 = lat_inflow_z2 * cellarea() / 1000
        column_lat_inflow_m3 = in_latflow_z0_m3 + in_latflow_z1_m3 + in_latflow_z2_m3
        i_latflow_m3 = accuflux(self.ldd_subs, column_lat_inflow_m3)
        self.out_accu_i_latflow_m3_tss.sample(i_latflow_m3)

        out_latflow_z0_m3 = lat_outflow_z0 * cellarea() / 1000  # m3
        out_latflow_z1_m3 = lat_outflow_z1 * cellarea() / 1000
        out_latflow_z2_m3 = lat_outflow_z2 * cellarea() / 1000
        column_lat_outflow_m3 = out_latflow_z0_m3 + out_latflow_z1_m3 + out_latflow_z2_m3
        o_latflow_m3 = accuflux(self.ldd_subs, column_lat_outflow_m3)
        self.out_accu_o_latflow_m3_tss.sample(o_latflow_m3)

        # Net lateral flow
        net_latflow_z0_m3 = lat_netflow_z0 * cellarea() / 1000  # m3
        net_latflow_z1_m3 = lat_netflow_z1 * cellarea() / 1000
        net_latflow_z2_m3 = lat_netflow_z2 * cellarea() / 1000
        net_latflow_m3 = net_latflow_z0_m3 + net_latflow_z1_m3 + net_latflow_z2_m3
        n_latflow_m3 = accuflux(self.ldd_subs, net_latflow_m3)
        self.out_accu_n_latflow_m3_tss.sample(n_latflow_m3)

        # Percolation (only interested in the bottom-most layer, where mass leaves the model)
        percol_z2_m3 = percolation_z2 * cellarea() / 1000  # m3
        out_percol_m3 = accuflux(self.ldd_subs, percol_z2_m3)
        self.out_percol_z2_m3_tss.sample(out_percol_m3)

        # Evapotranspiration
        etp_z0_m3 = etp_z0 * cellarea() / 1000  # m3
        etp_z1_m3 = etp_z1 * cellarea() / 1000
        etp_z2_m3 = etp_z2 * cellarea() / 1000
        etp_m3 = etp_z0_m3 + etp_z1_m3 + etp_z2_m3
        out_etp_m3 = accuflux(self.ldd_subs, etp_m3)
        self.out_etp_m3_tss.sample(out_etp_m3)

        # Change in storage
        ch_storage_m3 = ch_storage_z0_m3 + ch_storage_z1_m3 + ch_storage_z2_m3
        out_ch_storage_m3 = accuflux(self.ldd_subs, ch_storage_m3)
        self.out_ch_storage_m3_tss.sample(out_ch_storage_m3)

        # Need to subtract: of_latflow_m3
        global_mb_water = tot_rain_m3 - out_runoff_m3 - out_percol_m3 - out_etp_m3 + n_latflow_m3 - o_drain_z1_m3 - out_ch_storage_m3
        self.global_mb_water_tss.sample(global_mb_water)

        cell_vol_z0 = self.theta_z0 * self.z0 * cellarea() / 1000
        cell_vol_z1 = self.theta_z1 * self.z1 * cellarea() / 1000
        cell_vol_z2 = self.theta_z2 * self.z2 * cellarea() / 1000
        cell_vol_tot_m3 = cell_vol_z0 + cell_vol_z1 + cell_vol_z2
        vol_tot_m3 = accuflux(self.ldd_subs, cell_vol_tot_m3)
        self.storage_m3_tss.sample(vol_tot_m3)

        # Runoff + accu_latflow + accu_drainage + of_latflow_m3
        in_vol_disch_m3 = out_runoff_m3 + i_latflow_m3
        out_vol_disch_m3 = out_runoff_m3 + o_latflow_m3
        net_vol_disch_m3 = out_runoff_m3 + n_latflow_m3 + o_drain_z1_m3 + of_latflow_m3

        self.i_Q_m3_tss.sample(in_vol_disch_m3)
        self.o_Q_m3_tss.sample(out_vol_disch_m3)
        self.n_Q_m3_tss.sample(net_vol_disch_m3)

        if self.PEST:
            ######################
            # Pesticide Balance
            ######################
            # Applied mg on catchment
            appl_catch_mg = accuflux(self.ldd_subs, mass_applied)
            cum_appl_catch_mg = accuflux(self.ldd_subs, self.cum_appl_mg)
            # Todo: Check if below same result
            # cum_appl_catch_mg = upstream(self.ldd_subs, self.cum_appl_mg)

            # Loss to run-off
            out_runoff_mg = accuflux(self.ldd_subs, z0_mass_runoff["mass_runoff"])
            cum_out_runoff_mg = accuflux(self.ldd_subs, self.cum_runoff_mg)

            # Loss to air/volatilized
            if self.currentTimeStep() in self.app_days:
                m_vol = z0_mass_volatilized.get("mass_loss", self.zero_map)
                out_volat_mg = accuflux(self.ldd_subs, m_vol)
                # cum_out_volat_mg = accuflux(self.ldd, ...)
            else:
                out_volat_mg = self.zero_map

            # Loss to leaching
            out_leach_mg = accuflux(self.ldd_subs, z2_mass_leached["mass_leached"])
            cum_out_leach_mg = accuflux(self.ldd_subs, self.cum_leached_mg_z2)

            # Loss to lateral flux
            # ... per time step
            mass_latflux_mg = (z0_mass_latflux["net_mass_latflux"] +
                               z1_mass_latflux["net_mass_latflux"] +
                               z2_mass_latflux["net_mass_latflux"])
            out_latflux_mg = accuflux(self.ldd_subs, mass_latflux_mg)
            # ... to date (cumulative)
            tot_latflux_mg = self.cum_latflux_mg_z0 + self.cum_latflux_mg_z1 + self.cum_latflux_mg_z2
            cum_out_latflux_mg = accuflux(self.ldd_subs, tot_latflux_mg)

            # Change in mass storage
            # ... per time step
            ch_storage_mg = ch_storage_z0_mg + ch_storage_z1_mg + ch_storage_z2_mg
            tot_ch_storage_mg = accuflux(self.ldd_subs, ch_storage_mg)
            # ... to date (cumulative)
            cum_ch_storage_mg = ch_storage_mg - self.pest_ini_storage_mg
            cum_tot_ch_storage_mg = accuflux(self.ldd_subs, cum_ch_storage_mg)

            global_mb_pest = appl_catch_mg - out_runoff_mg - out_leach_mg + out_latflux_mg - tot_ch_storage_mg - out_volat_mg
            global_mb_pest_cum = cum_appl_catch_mg - cum_out_runoff_mg - cum_out_leach_mg + cum_out_latflux_mg - cum_tot_ch_storage_mg  # - cum_out_volat_mg
            self.global_mb_pest_tss.sample(global_mb_pest)

        self.jd_cum += self.jd_dt  # updating JDcum, currently dt = 1 day

        # Analysis (NASH Discharge)
        q_obs = timeinputscalar('q_obs_m3day.tss', nominal("outlet_true"))
        rest_obs = tot_rain_m3/q_obs
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

        self.tot_drain_m3 += ifthenelse(q_obs >= 0, o_drain_z1_m3, 0)
        self.tot_accu_drain_m3_tss.sample(self.tot_drain_m3)

        self.tot_ilf_m3 += ifthenelse(q_obs >= 0, i_latflow_m3, 0)
        self.tot_accu_i_latflow_m3_tss.sample(self.tot_ilf_m3)
        self.tot_olf_m3 += ifthenelse(q_obs >= 0, o_latflow_m3, 0)
        self.tot_accu_o_latflow_m3_tss.sample(self.tot_olf_m3)
        self.tot_nlf_m3 += ifthenelse(q_obs >= 0, n_latflow_m3, 0)
        self.tot_accu_n_latflow_m3_tss.sample(self.tot_nlf_m3)
        self.tot_of_m3 += ifthenelse(q_obs >= 0, of_latflow_m3, 0)
        self.tot_accu_of_latflow_m3_tss.sample(self.tot_of_m3)


        # Daily Maps
        # self.report(vol_disch_m3, 'q')  # discharge map for all cells (accuflux)
        # self.report(vol_tot_m3, 'storage')  # Total catchment storage m3

        # Produce only one map for every realization
        if self.currentTimeStep() == 279:
            self.report(self.theta_z0, 'z0theta')
            self.report(self.theta_z1, 'z1theta')
            self.report(self.theta_z2, 'z2theta')

        if self.currentTimeStep() == self.app_days[0]+1:
            self.report(self.delta_z0, 'z0dC')
            self.report(self.pestmass_z0, 'z0mss')

        # Visualization
        # aguila --scenarios='{1,2,3,4,5}' --multi=1x5  --timesteps=[170,280,1] q outlet_true.map

        # Nash
        # aguila 1\res_nash_q_m3.tss 2\res_nash_q_m3.tss 3\res_nash_q_m3.tss
        # aguila 6\res_nash_q_m3.tss 7\res_nash_q_m3.tss 8\res_nash_q_m3.tss
        # aguila 1\res_nash_q_m3.tss 6\res_nash_q_m3.tss
        # Tot discharge & Daily discharge
        # aguila 4\res_q_obs_tot_m3.tss 4\res_q_sim_tot_m3.tss 8\res_accuVol_m3.tss q_obs_m3day.tss
        # aguila 1\res_accuLatflow_m3.tss 1\res_accuVol_m3.tss q_obs_m3day.tss


    def postmcloop(self):
        pass
        #names = ["q"]  # Discharge, Nash_Discharge
        #mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        # aguila --timesteps=[170,280,1] q-ave q-var outlet_true.map
        # percentiles = [0.25, 0.5, 0.75]
        # mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())
        # aguila --quantiles=[0.25,0.75,0.25] --timesteps=[170,280,1] q


nrOfSamples = 2  # Samples are each a MonteCarlo realization
firstTimeStep = 175
nTimeSteps = 180
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

# Current test, with Permeable bottom layer!

# Implementing first sensitivity analysis via mc loop
# Initial parameters are changed by 25%
# Discharge change is evaluated:
# a) Temporally via tss files (stored in folder realizations)
# b) As total at the end of the simulation. Total for sim, is only counted if, obs is not NA.
# Parameters tested are:
# Baseline
# 1) drain_coef = 0.8 & latflow_coef (c1) = 0.25
# Case A
# 2) drain_coef + 0.15
# 3) drain_coef - 0.15
# 4) latflow_coef (c1) - 0.15
# 5) latflow_coef (c1) + 0.15

# Case B  - basement
# 6) latflow_coef (c2) += 0.15
# 7) latflow_coef (c2) -= 0.15

# Case C
# 8) z0 * 1.15
# 9) z0 * 0.75

# Case D
# 10) *1.15 Ksat
# 11) *0.75


# Checked-success:
# Testing whether Nash calculations work
#   -> Yes

# Checked-Failed:
# Check if postmcloop() makes an average, if only a map exists for 1 time step across realizations?
#   -> No, results in error!!
# Check if Nash calculation works with NA values on "q_obs"
#   -> Failed, with new import (unsure about NA values), now removed col.names
