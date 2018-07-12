# -*- coding: utf-8 -*-

from hydro_v2 import *
from nash import *
from pesti_v2 import *
from soil_samples import *
from test_suite import *
from balance import *

print(os.getcwd())

global state
state = -1

"""

This model inherits fro, Morris,

Changes:

- not a montecarlo
- parameters are drawn from the argument in the objective function
- TSS written needs to be only discharge

ATTENTION: Model still does not transport mass from Baseflow!!!

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
#     return str(obs_dict['data'][timestep-2][2])


def get_state(old_state):
    new_state = old_state + 1
    global state
    state = new_state
    return state


def start_jday():
    start_sim = 166  # 213 # 166
    return start_sim


class BeachModel(DynamicModel):
    def setDebug(self):
        pass

    def __init__(self, cloneMap, params):
        DynamicModel.__init__(self)
        setclone(cloneMap)

        self.params = params

    def initial(self):
        # This section includes all non-stochastic parameters.
        # Get initial parameters, make a dictionary of the raw file.
        import csv
        ini_path = 'initial.csv'
        self.ini_param = {}  # Dictionary to store the values
        with open(ini_path, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                self.ini_param[row[0].strip()] = float(row[1])

        self.num_layers = int(self.ini_param.get("layers"))
        # Hydrological scenarios
        self.bsmntIsPermeable = False  # basement percolation (DP)
        self.ADLF = True
        self.fixed_dt50 = False
        self.bioavail = True

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

        self.PEST = True
        self.TRANSPORT = True
        # Run fate processes
        self.ROM = True
        self.LCH = True
        self.ADRM = True
        self.LFM = True
        self.DEG = True

        self.TEST_LCH = False
        self.TEST_LFM = False
        self.TEST_DEG = False

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
        # self.report(self.mask, 'mask')
        self.fa_cr = readmap("farm_burn_v3")  # Contains codes to assign appropriate dosage
        self.aging = deepcopy(self.zero_map)  # Cumulative days after application on each pixel

        self.outlet_multi = self.readmap("out_multi_nom_v3")  # Multi-outlet with 0 or 2
        self.is_outlet = boolean(self.outlet_multi == 1)
        # self.report(self.is_outlet, 'OUT')

        importPlotMaps(self)

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

        self.q_m3day_mean = scalar(self.ini_param.get("ave_outlet_q_m3day"))

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
        self.out_baseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", self, nominal("outlet_v3"),
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
        self.nash_outlet_conc_tss = TimeoutputTimeseries("resNash_outConc_ugL", self, nominal("outlet_v3"),
                                                         noHeader=False)
        self.nash_outlet_iso_tss = TimeoutputTimeseries("resNash_outIso_delta", self, nominal("outlet_v3"),
                                                        noHeader=False)

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

        # Morris_error tests
        m_state = get_state(state)  # First run will return state = 0

        names = ['z3_factor',
                 'cZ0Z1', 'cZ',
                 'c_adr',
                 'k_g',
                 'gamma01', 'gammaZ',
                 'f_transp',
                 'FCz2', 'FCz',
                 'SATz2', 'SATz',
                 'WPz2', 'WPz'
                 ]

        z3_factor = self.mask * self.params[names.index('z3_factor')]

        """ Physical parameters for each layer """
        self.gamma = []  # coefficient to calibrate Ksat1
        self.c_lf = []
        for layer in range(self.num_layers):
            if layer < 2:
                self.gamma.append(self.mask * self.params[names.index('gamma01')])  # percolation coefficient
                self.c_lf.append(self.mask * self.params[names.index('cZ0Z1')])  # subsurface flow coefficient)
            else:
                self.gamma.append(self.mask * self.params[names.index('gammaZ')])  # percolation coefficient
                self.c_lf.append(self.mask * self.params[names.index('cZ')])

        self.c_adr = self.mask * self.params[names.index('c_adr')]
        self.k_g = self.mask * self.params[names.index('k_g')]  # [days]
        self.gw_factor = 1 - z3_factor  # Represents bottom-most portion of bottom layer
        self.f_transp = self.mask * self.params[names.index("f_transp")]  # [-]

        self.drainage_layers = [False, False, True, False, False]  # z2 turned on!!

        """
        Hydro Maps
        """

        # theta_sat_z2 = self.zero_map + readmap("thSATz2")  # 0.63  scalar(self.ini_param.get("sat_z2z3"))
        # # scalar(self.ini_param.get("sat_z2z3")) + mapnormal() * 0.04  # mean + 1SD(X)*0.04 = mean + (0.002)**0.5
        # theta_fcap_z2 = self.zero_map + readmap("thFCz2") # scalar(self.ini_param.get("fc_z2z3")) # => 0.38702
        # scalar(self.ini_param.get("fc_z2z3")) + mapnormal() * 0.04  # mean + 1SD(X)*0.04 = mean + (0.002)**0.5

        # Initial moisture (Final from model v1, Sept 30, 2016)
        self.theta = []
        self.theta_sat = []
        self.theta_fc = []
        self.theta_wp = []  # => 0.19
        for layer in range(self.num_layers):
            if layer < 2:
                self.theta_sat.append(deepcopy(self.zero_map))  # Dynamic TSS
                self.theta_fc.append(deepcopy(self.zero_map))
                self.theta_wp.append(scalar(self.ini_param.get("wp_z01")))
            elif layer == 2:
                self.theta_sat.append(self.mask * self.params[names.index('SATz2')])
                self.theta_fc.append(self.mask * self.params[names.index('FCz2')])
                self.theta_wp.append(self.mask * self.params[names.index('WPz2')])
            else:
                self.theta_sat.append(self.mask * self.params[names.index('SATz')])
                self.theta_fc.append(self.mask * self.params[names.index('FCz')])
                self.theta_wp.append(self.mask * self.params[names.index('WPz')])

            if start_jday() < 166:
                name = 'd14_theta_z' + str(layer)
                self.theta.append(readmap(name))
            else:
                name = 'd166_theta_z' + str(layer)
                self.theta.append(readmap(name))

            self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                           self.theta[layer])

        """ Soil Properties """
        self.p_bZ = scalar(self.ini_param.get("p_bZ"))  # Soil bulk density (g/cm^3)
        self.act_e = scalar(self.ini_param.get("activation_e"))  # Metolachlor Ea = 23.91 KJ/mol; @Jaikaew2017
        self.r_gas = scalar(self.ini_param.get("r_gas"))  # R = 8.314 J / mol Kelvin,

        """
        Layer depths
        """

        self.layer_depth = []
        self.tot_depth = deepcopy(self.zero_map)
        bottom_depth = deepcopy(self.zero_map)
        for layer in range(self.num_layers):
            if layer < self.num_layers - 2:  # 5 - 1 = 3 (i.e. z0,z1,z2)
                self.layer_depth.append(self.zero_map +
                                        scalar(self.ini_param.get('z' + str(layer))))
                self.tot_depth += self.layer_depth[layer]
                # self.report(self.layer_depth[layer], 'DepthZ' + str(layer))
            elif layer < self.num_layers - 1:  # 5 - 2 = 4 (i.e. z3)
                bottom_depth = (self.datum_depth +  # total height
                                scalar(self.ini_param.get('z' + str(layer))) + 100  # plus a min-depth
                                - self.tot_depth)
                self.layer_depth.append(bottom_depth * z3_factor)  # minus:(z0, z1, z2)*decreasing depth factor
                self.tot_depth += self.layer_depth[layer]
                # self.report(self.layer_depth[layer], 'DepthZ' + str(layer))
            else:  # Basement Layer = n5  (z4)
                self.layer_depth.append(bottom_depth * self.gw_factor)  # minus:(z0, z1, ...)*decreasing depth factor
                self.tot_depth += self.layer_depth[layer]
                # self.report(self.layer_depth[layer], 'DepthZ' + str(layer))

            if self.TEST_depth:
                checkLayerDepths(self, layer)

        self.smp_depth = self.layer_depth[0]

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
        self.dd_max = (scalar(2500) * self.p_bZ) / (self.p_bZ + 686 * exp(-5.63 * self.p_bZ))

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

        # Nash discharge
        self.days_cum = 0  # Track no. of days with data
        self.q_diff = 0
        self.q_var = 0
        self.q_obs_cum = 0
        self.q_sim_cum = 0  # net total disch
        self.q_sim_ave = 0

        # Nash concentration outlet
        self.out_conc_diff = 0
        self.out_conc_var = 0
        self.out_lnconc_diff = 0
        self.out_lnconc_var = 0

        # Nash isotopes outlet
        self.out_iso_diff = 0
        self.out_iso_var = 0

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

        # Need initial states to compute change in storage after each run
        self.theta_ini = deepcopy(self.theta)

    def dynamic(self):

        jd_sim = self.jd_start + self.jd_cum
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
        min_root_depth = scalar(150)  # Seeding depth [mm], Allen advices: 0.15 to 0.20 m

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
        """
        # Basement layers defined under initial()
        self.theta_sat[0] = timeinputscalar('thetaSat_agr.tss', nominal(self.landuse))  # saturated moisture # [-]
        self.theta_sat[1] = deepcopy(self.theta_sat[0])
        self.theta_fc[0] = timeinputscalar('thetaFC_agr.tss', nominal(self.landuse))  # * self.fc_adj  # field capacity
        self.theta_fc[1] = deepcopy(self.theta_fc[0])

        if self.TEST_thProp:
            checkMoistureProps(self, self.theta_sat, 'aSATz')
            checkMoistureProps(self, self.theta_fc, 'aFCz')

        self.p_bAgr = timeinputscalar('p_b_agr.tss', nominal(self.landuse))
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
        # root_depth_tot2 = timeinputscalar('height.tss', nominal(self.landuse)) * self.root_adj
        # root_depth_tot2 *= 10 ** 3  # Convert to mm

        root_depth_tot = ifthenelse(jd_sim > jd_mid, max_root_depth,  # Passed development stage
                                    ifthenelse(jd_sim > jd_plant,  # Planted, but before attaining max root depth
                                               min_root_depth + ((max_root_depth - min_root_depth) *
                                                                 ((jd_sim - jd_plant) / (jd_mid - jd_plant))),
                                               scalar(0)))  # Before planting
        # self.report(root_depth_tot, 'RDtot1')

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

        # Get potential evapotranspiration for all layers
        etp_dict = getPotET(self, sow_yy, sow_mm, sow_dd, root_depth_tot, min_root_depth,
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

        # self.report(pot_transpir, 'aPotTRA')
        # self.report(pot_evapor, 'aPotEVA')

        # Not in use for water balance, but used to estimate surface temp due to bio-cover.
        frac_soil_cover = etp_dict["f"]
        # self.report(frac_soil_cover, 'aFracCV')

        bio_cover = getBiomassCover(self, frac_soil_cover)
        # bcv should range 0 (bare soil) to 2 (complete cover)
        # self.report(bio_cover, 'aBCV')

        """
        Infiltration, runoff, & percolation (all layers)
        """
        infil_z0 = deepcopy(self.zero_map)
        runoff_z0 = deepcopy(self.zero_map)
        percolation = []

        permeable = True
        for layer in range(self.num_layers):
            if layer == 0:  # Layer 0
                z0_IRO = getTopLayerInfil(self, precip, CN2, crop_type,
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

                excess = max(self.theta[layer] - self.theta_sat[layer], scalar(0))
                if mapmaximum(excess) > 0:
                    val = float(mapmaximum(excess))
                    self.theta[layer] = ifthenelse(self.theta[layer] > self.theta_sat[layer], self.theta_sat[layer],
                                                   self.theta[layer])
                    if float(val) > float(1e-02):
                        print("Corrected Percolation(), SAT was exceeded, layer " + str(layer) + ' by ' + str(val))

                excessLj = max(self.theta[layer + 1] - self.theta_sat[layer + 1], scalar(0))
                if mapmaximum(excessLj) > 0:
                    val = float(mapmaximum(excessLj))
                    self.theta[layer + 1] = ifthenelse(self.theta[layer + 1] > self.theta_sat[layer + 1],
                                                       self.theta_sat[layer + 1],
                                                       self.theta[layer + 1])
                    if float(val) > float(1e-02):
                        print("Corrected Percolation(), SAT exceeded, layer " + str(layer + 1) + ' by ' + str(val))

                percolation.append(getPercolation(self, layer, k_sat[layer], isPermeable=permeable))  # [mm]
                water_flux = infil_z0 + infil_z1 + percolation[layer]

                SW1b = self.theta[layer] * self.layer_depth[layer] - percolation[layer]
                self.theta[layer] = SW1b / self.layer_depth[layer]

                # Discharge due to runoff at the outlet
                runoff_m3 = runoff_z0 * cellarea() / 1000  # m3
                accu_runoff_m3 = accuflux(self.ldd_subs, runoff_m3)
                out_runoff_m3 = areatotal(accu_runoff_m3, self.outlet_multi)

            else:  # Layers 2, 1, 3 & 4
                if layer == (self.num_layers - 1):
                    permeable = self.bsmntIsPermeable

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

                if mapmaximum(self.theta[layer]) > mapmaximum(self.theta_sat[layer]):
                    val = float(mapmaximum(self.theta[layer])) - float(mapmaximum(self.theta_sat[layer]))
                    if float(val) > float(1e-06):
                        print(
                                "Error at Percolation() layers: 2, 1, 3, SAT exceeded, layer " + str(
                            layer) + ' by ' + str(
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

                SW3 = self.theta[layer] * self.layer_depth[layer] - percolation[layer]
                self.theta[layer] = SW3 / self.layer_depth[layer]

            if mapminimum(self.theta[layer]) < 0:
                print('Error at Percolation, Layer ' + str(layer))

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

        SW4 = self.theta[adr_layer] * self.layer_depth[adr_layer] - cell_drainge_outflow
        self.theta[adr_layer] = SW4 / self.layer_depth[adr_layer]

        if mapminimum(self.theta[adr_layer]) < 0:
            print('Error at ADR, Layer ' + str(adr_layer))

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
                act_evaporation_layer = getActualEvap(self, layer, pot_evapor_layer, run=self.ETP)
            # Evaporation
            SW5 = self.theta[layer] * self.layer_depth[layer] - act_evaporation_layer
            self.theta[layer] = SW5 / self.layer_depth[layer]
            evap.append(act_evaporation_layer)
            evap_m3 += evap[layer] * cellarea() / 1000  # m3

            # Transpiration
            act_transpir_layer = getActualTransp(self, layer, root_depth_tot, root_depth[layer],
                                                 pot_transpir, depletable_water, run=self.ETP)
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
        for layer in range(self.num_layers):
            if layer < (self.num_layers - 1):
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

            else:  # Basement layer
                latflow_outlet_mm.append(deepcopy(self.zero_map))
                latflow_cell_mm.append(deepcopy(self.zero_map))

            cell_lat_outflow_m3 += latflow_outlet_mm[layer] * cellarea() / 1000  # m3, only out!

            if mapminimum(self.theta[layer]) < 0:
                print('Error at LF, Layer ' + str(layer))

        # Lateral flow (Outlet discharge)
        outlet_latflow_m3 = areatotal(cell_lat_outflow_m3, self.outlet_multi)  # Only outlet cells

        # Baseflow
        # Considering part of basement layer as completely saturated
        # baseflow_mm = (self.theta_sat[-2] * self.layer_depth[-2] * self.gw_factor) / self.k_g  # [mm/d]
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

        for layer in range(self.num_layers):
            # Temperature
            temp_dict = getLayerTemp(self, layer, bio_cover, temp_bare_soil)
            self.temp_surf_fin = temp_dict["temp_surface"]
            self.temp_fin[layer] = temp_dict["temp_layer"]

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

        if self.TEST_theta and self.currentTimeStep() % 2 == 0:
            checkMoisture(self, self.theta, 'athz')

        """
        Nash computations
        """
        # Total days with data (needed for mean calculations)
        self.days_cum += ifthenelse(q_obs >= 0, scalar(1), scalar(0))

        # Analysis (NASH Discharge)
        reportNashHydro(self, q_obs, tot_vol_disch_m3)


mini_test = True
if mini_test:
    params = [0.99,
              1, 1,
              1,
              3650,
              1, 1,
              1,
              0.41, 0.41,
              0.61, 0.53,
              0.19, 0.19
              ]
    firstTimeStep = start_jday()  # 166 -> 14/03/2016
    nTimeSteps = 200  # 360
    myAlteck16 = BeachModel("clone_nom.map", params)
    dynamicModel = DynamicFramework(myAlteck16, lastTimeStep=nTimeSteps,
                                    firstTimestep=firstTimeStep)  # an instance of the Dynamic Framework

    dynamicModel.run()
