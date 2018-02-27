
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
        self.zero_map = self.dem - self.dem  # Zero map to generate scalar maps

        self.slope_rad = sin(atan(max(slope(self.dem), 0.001)))  # Slope in radians
        self.dem_route = self.readmap("dem_ldd")  # To route surface run-off
        self.ldd = self.readmap("ldd")  # To route subsurface lateral flow & build TWI
        # TODO: Check best ldd ?
        # self.ldd = lddcreate(self.dem, 10, 1e31, 1e31, 1e31)  # There are differences btw. ldd's

        self.tot_depth = (self.dem - mapminimum(self.dem))*scalar(10**3)  # mm
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
        self.wetness = ln(self.up_area / tan(self.slope_rad))

        """
        Pesticides Maps
        """
        self.smback = self.readmap("sm_back")  # Map initial/background pesticide in soil (mg/m2)
        self.app0 = self.readmap("app0")  # Pesticide applied (mg/m2) on Julian day 177 (March 25, 2016).
        self.app1 = self.readmap("app1")  # Pesticide applied (mg/m2) on Julian day 197 (April 14, 2016).
        self.app2 = self.readmap("app3")  # Pesticide applied (mg/m2) on Julian day 238 (May 25, 2016).

        # Use map algebra to produce a initial signature map,
        # where app1 > 0, else background sig. (plots with no new mass will be 0)
        # ATT: Need to do mass balance on addition of new layer.
        self.app0delta = ifthenelse(self.app0 > 0, scalar(-32.3), scalar(-23.7))
        self.app1delta = ifthenelse(self.app2 > 0, scalar(-32.3), scalar(-23.7))
        self.app2delta = ifthenelse(self.app3 > 0, scalar(-32.3), scalar(-23.7))

        # Convert mg/m2 -> mg
        self.pestmass_z0 = self.smback * self.cell_area  # mg
        self.pestmass_z0_ini = self.pestmass_z0
        self.pestmass_z1 = self.smback * self.cell_area  # mg
        self.pestmass_z1_ini = self.pestmass_z1
        self.pestmass_z2 = self.smback * self.cell_area  # mg
        self.pestmass_z2_ini = self.pestmass_z2

        """
        Temperature Maps and params
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
        self.out_runoff_m3_tss = TimeoutputTimeseries("Out_runoff_m3", self, "outlet.map", noHeader=False)
        self.out_latflow_m3_ss = TimeoutputTimeseries("Out_latflow_m3", self, "outlet.map", noHeader=False)
        self.out_percol_m3_tss = TimeoutputTimeseries("Out_percol_m3", self, "outlet.map", noHeader=False)
        self.out_etp_m3_tss = TimeoutputTimeseries("Out_etp_m3", self, "outlet.map", noHeader=False)
        self.tot_ch_storage_m3_tss = TimeoutputTimeseries("Tot_ch_storage_m3", self, "outlet.map", noHeader=False)

        # Transects and detailed soils
        ###########
        self.obs_trans = self.readmap("observe_trans")
        self.obs_detail = self.readmap("observe_detail")

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