
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
            Landscape & Hydro
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
            Pesticides
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
        # Output & Observations
        """
        self.obs = self.readmap("observe")




