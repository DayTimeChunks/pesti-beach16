from pcraster._pcraster import *
from pcraster.framework import *
import os
import time

import hydro
import pesti

print(os.getcwd())


class BeachModel(DynamicModel):
    def setDebug(self):
        pass

    def __init__(self, cloneMap):
        DynamicModel.__init__(self)
        setclone(cloneMap)

        # dem = self.dem

    def initial(self):
        """ Physical parameters """
        self.c1 = 0.25  # subsurface flow coefficient
        self.c2 = 0.25
        self.drain_coef = 0.8063  # drainage coefficient
        self.s1 = 1  # coefficient to calibrate Ksat1
        self.s2 = 0.5  # coefficient to calibrate Ksat2
        self.k = scalar(0.03)  # coefficient of declining LAI in end stage

        """ Soil Properties """
        self.p_b = 1.4  # Soil bulk density (g/cm^3)
        self.f_oc = 0.021  # Organic Carbon in soil without grass (kg/kg)

        """
        Sorption parameters
        """
        # K_oc - S-metolachlor (K_oc in ml/g)
        # Marie's thesis: log(K_oc)= 1.8-2.6 [-] -> k_oc = 63 - 398 (Alletto et al., 2013).
        # Pesticide Properties Database: k_oc = 120 ml/g (range: 50-540 mL/g)
        self.k_oc = 120  # ml/g
        self.k_d = self.k_oc * self.f_oc  # Dissociation coefficient K_d (mL/g = L/Kg)

        # Pesticide Properties Database states :
        # K_d=0.67; but here
        # K_d=120*0.021 = 2.52;  (L/kg)
        # Difference will lead to higher retardation factor (more sorption)

        """
        Volatilization parameters
        """
        # Henry's constant @ 20 C (dimensionless, Metolachlor)
        # https://www.gsi-net.com
        self.k_h = 3.1326141504e-008

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
        self.dt_50_ref = 15  # S-met (days)

        self.temp_ref = 20  # Temp.  reference

        self.beta_temperature = 1  # Need to find a correct value for 'B', exponent in moisture dependency (Dairon)
        self.alpha_temperature = 54000 / float(8314)  # Need to confirm units Ea = 54000 KJ/mol; R = 8.314 J/mol/Kelvin

        """
        Isotopes
        """
        # VPDB
        self.r_standard = 0.0112372
        self.alpha_iso = 1  # 1 = no fractionation

        """
        Loading maps
        """

        """
        Landscape & Hydro Maps
        """
        self.dem = self.readmap("dem_slope")  # 192 - 231 m a.s.l
        self.dem_route = self.readmap("dem_ldd")  # To route surface run-off
        self.zero_map = self.dem - self.dem  # Zero map to generate scalar maps
        mask = self.dem / self.dem

        self.ldd_surf = lddcreate(self.dem_route, 1e31, 1e31, 1e31, 1e31)  # To route runoff
        self.ldd_subs = lddcreate(self.dem, 1e31, 1e31, 1e31, 1e31)  # To route lateral flow & build TWI

        self.tot_depth = (self.dem - mapminimum(self.dem)) * scalar(10 ** 3)  # mm
        self.z0 = self.zero_map + 10  # mm
        self.z1 = self.zero_map + 140  # mm
        self.z2 = self.tot_depth - self.z0 - self.z1  # mm

        # Initial moisture (arbitrary, Oct, 2015)
        self.theta_z0 = self.zero_map + 0.25  # map of initial soil moisture in top layer (-)
        self.theta_z1 = self.zero_map + 0.25
        self.theta_z2 = self.zero_map + 0.25

        # Need initial states to compute change in storage after each run
        self.theta_z0_ini = self.theta_z0
        self.theta_z1_ini = self.theta_z1
        self.theta_z2_ini = self.theta_z2

        self.outlet = self.readmap("outlet")
        self.landuse = self.readmap("landuse")

        # Topographical Wetness Index
        self.up_area = accuflux(self.ldd, self.cell_area)
        self.slope_rad = sin(atan(max(slope(self.dem), 0.001)))  # Slope in radians
        self.wetness = ln(self.up_area / tan(self.slope_rad))

        """
        Pesticides Maps
        """
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
        double = 2.0  # ~ Dosage for corn when growing beet
        d_gold = 915 * 10 ** 6  # ug/L S-met
        m_gold = 960 * 10 ** 6  # ug/L

        # Dosages # L/Ha * 1Ha/1000m2 = L/m2
        d_beet = None
        d_corn = 2.1 * 1 / 10 ** 4  # L/Ha * 1 Ha / 10000 m2
        m_beet = 0.6 * 1 / 10 ** 4 * double
        m_corn = 2.0 * 1 / 10 ** 4
        m_beet_Friess = 0.6 * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)
        m_beet_Mathis = 0.6 * 1 / 10 ** 4 * (double + 1)  # (Likely larger dosage, early in the season)

        # Assign dosages based on Farmer-Crop combinations [ug/m2]
        fa_cr = readmap("farmer_crop")  # Contains codes to assign appropriate dosage
        app_conc = (  # [ug/m2]
            ifthenelse(fa_cr == 1111,  # 1111 (Friess, Beet)
                       m_beet_Friess * m_gold * mask,
                       ifthenelse(fa_cr == 1122,  # 1112 (Friess-Corn),
                                  m_corn * m_gold * mask,
                                  ifthenelse(fa_cr == 1212,  # 1212 (Speich-Corn),
                                             m_corn * m_gold * mask,
                                             ifthenelse(fa_cr == 1312,  # 1312 (Mahler-Corn),
                                                        m_corn * m_gold * mask,
                                                        ifthenelse(fa_cr == 1412,  # 1412 (Schmitt-Corn)
                                                                   d_corn * d_gold * mask,
                                                                   ifthenelse(fa_cr == 1511,  # 1511 (Burger-Beet)
                                                                              m_beet * m_gold * mask,
                                                                              # 1711 (Mathis-Beet),
                                                                              ifthenelse(fa_cr == 1711,
                                                                                         m_beet_Mathis * m_gold * mask,
                                                                                         # 1611 (Kopp-Beet)
                                                                                         ifthenelse(
                                                                                             fa_cr == 1611,
                                                                                             m_beet * m_gold * mask,
                                                                                             0 * mask))))))))
        )
        # Pesticide applied (ug/m2) on Julian day 177 (March 25, 2016).
        # March 26th, Friess and Mathis
        self.app1 = ifthenelse(fa_cr == 1111, 1*app_conc,
                               # 1111 (Friess, Beet), 1112 (Friess-Corn),
                               ifthenelse(fa_cr == 1112, 1*app_conc,
                                          ifthenelse(fa_cr == 1711, 1*app_conc,  # 1711 (Mathis-Beet)
                                                     0*app_conc)))
        # Pesticide applied (ug/m2) on Julian day 197 (April 14, 2016).
        # April 14, Kopp and Burger
        self.app2 = ifthenelse(fa_cr == 1511, 1*app_conc,  # 1511 (Burger-Beet)
                               ifthenelse(fa_cr == 1611, 1*app_conc,  # 1611 (Kopp-Beet),
                                          0*app_conc))

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
        self.pestmass_z0 = self.smback * self.cell_area  # mg
        self.pestmass_z0_ini = self.pestmass_z0
        self.pestmass_z1 = self.smback * self.cell_area  # mg
        self.pestmass_z1_ini = self.pestmass_z1
        self.pestmass_z2 = self.smback * self.cell_area  # mg
        self.pestmass_z2_ini = self.pestmass_z2

        """
        Temperature maps and params
        """
        self.lag = 0.8  # lag coefficient (-), 0 < lag < 1; -> in SWAT, lag = 0.80
        # Generating initial surface temp map (15 deg is arbitrary)
        self.temp_z0_fin = self.zero_map + 15
        self.temp_z1_fin = self.zero_map + 15
        self.temp_z2_fin = self.zero_map + 15
        self.temp_surf_fin = self.zero_map + 15

        # Maximum damping depth (dd_max)
        # The damping depth (dd) is calculated daily and is a function of max. damping depth (dd_max), (mm):
        self.dd_max = (2500 * self.p_b) / (self.p_b + 686 * exp(-5.63 * self.p_b))

        # TODO
        # Average Annual air temperature (celcius - Layon!! Not Alteckendorf yet!!)
        self.temp_ave_air = 12.2  # 12.2 is for Layon

        """
        Output & Observations (tss and observation maps)
        """
        # Output time series (tss)

        # Outlet
        ###########
        # Pesticide
        self.out_mb_pest_tss = TimeoutputTimeseries("out_mb_pest", self, "outlet.map", noHeader=False)
        self.out_delta_tss = TimeoutputTimeseries("out_delta", self, "outlet.map", noHeader=False)
        self.out_mass_tss = TimeoutputTimeseries("out_mass", self, "outlet.map", noHeader=False)

        # Hydro
        self.out_vol_m3_tss = TimeoutputTimeseries("out_vol_m3", self, "outlet.map", noHeader=False)
        self.out_runoff_m3_tss = TimeoutputTimeseries("out_runoff_m3", self, "outlet.map", noHeader=False)
        self.out_latflow_m3_ss = TimeoutputTimeseries("out_latflow_m3", self, "outlet.map", noHeader=False)
        self.out_percol_m3_tss = TimeoutputTimeseries("out_percol_m3", self, "outlet.map", noHeader=False)
        self.out_etp_m3_tss = TimeoutputTimeseries("out_etp_m3", self, "outlet.map", noHeader=False)
        self.tot_ch_storage_m3_tss = TimeoutputTimeseries("oot_ch_storage_m3", self, "outlet.map", noHeader=False)

        # Transects and detailed soils
        ###########
        self.obs_trans = self.readmap("weekly_smp")
        self.obs_detail = self.readmap("detailed_smp")

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
        self.jd_cum = 0
        self.jd_dt = 1  # Time step size (days)

        def dynamic(self):
            jd_sim = self.jd_start + self.jd_cum
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
            len_grow_stage_ini = lookupscalar('croptable.tbl', 5, fields)  # old: Lini. length of initial crop growth stage
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
            k_sat_z0z1 = lookupscalar('croptable.tbl', 22, fields)  # saturated conductivity of the first layer
            k_sat_z2 = lookupscalar('croptable.tbl', 23, fields)  # saturated conductivity of the second layer
            CN2 = lookupscalar('croptable.tbl', 24, fields)  # curve number of moisture condition II

            # adjusting K_sat
            k_sat_z0z1 *= self.s1
            k_sat_z2 *= self.s2

            "Assign time-series data to spatial location, map is implicitly defined as the clonemap."
            precip = timeinputscalar('rain.tss', 1)  # daily precipitation data as time series (mm)
            temp_bare_soil = timeinputscalar('T_bare.tss', nominal('clone'))
            # T_bare, based on SWAT (Neitsch2009), p.43. Radiation term = (Rs*(1-albedo)-14)/20
            # -> see ET0.xls, SWAT doc. and Allen et al., 2006 (FAO65), Albedo = 0.05 (bare soil)

