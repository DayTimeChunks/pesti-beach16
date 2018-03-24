# -*- coding: utf-8 -*-

from pcraster._pcraster import *
from pcraster.framework import *
# import os
# import time


def get_conc_aq(model, layer, theta_sat, sorption_model="linear", gas=True):
    theta_layer = temp_layer = depth = mass_layer = None
    if layer == 0:
        theta_layer = model.theta_z0
        depth = model.z0  # mm
        mass_layer = model.pestmass_z0  # mg
    elif layer == 1:
        theta_layer = model.theta_z1
        depth = model.z1  # mm
        mass_layer = model.pestmass_z1  # mg
    elif layer == 2:
        theta_layer = model.theta_z2
        depth = model.z2  # mm
        mass_layer = model.pestmass_z2  # mg
    else:
        print("incorrect number of layers, raising error")
        raise NotImplementedError

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        retard_layer = 1

    if gas:
        # Leistra et al., 2001
        theta_gas = theta_sat - theta_layer
        conc_layer_aq = mass_layer / ((cellarea() * depth) *
                                      (theta_gas * model.k_h + theta_layer * retard_layer))  # mg/L
    else:
        # Whelan, 1987 # No gas phase considered
        conc_layer_aq = (mass_layer / cellarea()) / (theta_layer * retard_layer * depth)  # mg/L
    return conc_layer_aq


def get_conc_from_mass(model, layer,
                       theta_sat_z0z1, theta_sat_z2,
                       return_phase=None, gas=True, isotopes=True, mass=None):
    if layer == 0:
        theta_aq = model.theta_z0
        theta_sat = theta_sat_z0z1
        depth = model.z0  # mm
        mass = mass if isotopes else model.pestmass_z0
    elif layer == 1:
        theta_aq = model.theta_z1
        theta_sat = theta_sat_z0z1
        depth = model.z1  # mm
        mass = mass if isotopes else model.pestmass_z0
    elif layer == 2:
        theta_aq = model.theta_z2
        theta_sat = theta_sat_z2
        depth = model.z2  # mm
        mass = mass if isotopes else model.pestmass_z0
    else:
        print("incorrect number of layers, raising error")
        raise NotImplementedError

    if gas:
        theta_gas = theta_sat - theta_aq
        if return_phase == "aq":  # [mg pest/L water]
            conc_phase = mass / (cellarea() * depth *
                                 (theta_gas * model.k_h + theta_aq + model.p_b * model.k_d))
        elif return_phase == "ads":  # [mg pest/Kg soil]
            conc_phase = (mass * model.k_d) / (cellarea() * depth *
                                               (theta_gas * model.k_h + theta_aq + model.p_b * model.k_d))
        else:
            print("No phase model chosen from: ('aq' or 'ads')")
            raise NotImplementedError
    else:
        print("No implementation without gas available")
        raise NotImplementedError

    return conc_phase


def get_mass_from_aq(model, layer,
                     theta_sat_z0z1, theta_sat_z2,
                     conc_aq, gas=True, sorption_model="linear"):
    """
    Used to update mass due to equilibrium conditions after degradation
    :param model:
    :param layer:
    :param conc_aq:
    :param gas:
    :param sorption_model:
    :return:
    """
    if layer == 0:
        theta_aq_layer = model.theta_z0
        theta_sat = theta_sat_z0z1
        depth = model.z0  # mm
    elif layer == 1:
        theta_aq_layer = model.theta_z1
        theta_sat = theta_sat_z0z1
        depth = model.z1  # mm
    elif layer == 2:
        theta_aq_layer = model.theta_z2
        theta_sat = theta_sat_z2
        depth = model.z2  # mm
    else:
        print("incorrect number of layers, raising error")
        raise NotImplementedError

    if gas:
        if sorption_model == "freundlich":
            print("Freundlich not implemented, running Linear Sorption")
            # Leistra et al., 2001
            theta_gas = theta_sat - theta_aq_layer
            mass_layer = cellarea() * depth * (theta_gas * model.k_h * conc_aq  # gas
                                                    + theta_aq_layer * conc_aq  # aqueous
                                                    + model.p_b * model.k_d * conc_aq)  # sorbed
        else:
            # Leistra et al., 2001
            theta_gas = theta_sat - theta_aq_layer
            mass_layer = cellarea() * depth * (theta_gas * model.k_h * conc_aq  # gas
                                                    + theta_aq_layer * conc_aq  # aqueous
                                                    + model.p_b * model.k_d * conc_aq)  # sorbed
    else:
        print("Running Linear Sorption, no Gas Phase")
        # Whelan, 1987 # No gas phase considered
        mass_layer = cellarea() * depth * (theta_aq_layer * conc_aq
                                                + model.p_b * model.k_d * conc_aq)  # mg
    return mass_layer  # mg


def get_mass_from_ads(model, layer,
                      theta_sat_z0z1, theta_sat_z2,
                      conc_ads, gas=True, sorption_model="linear"):
    if layer == 0:
        theta_aq_layer = model.theta_z0
        theta_sat = theta_sat_z0z1
        depth = model.z0  # mm
    elif layer == 1:
        theta_aq_layer = model.theta_z1
        theta_sat = theta_sat_z0z1
        depth = model.z1  # mm
    elif layer == 2:
        theta_aq_layer = model.theta_z2
        theta_sat = theta_sat_z2
        depth = model.z2  # mm
    else:
        print("incorrect number of layers, raising error")
        raise NotImplementedError
    if gas:
        theta_gas = theta_sat - theta_aq_layer
        if sorption_model == "freundlich":
            print("Freundlich not implemented, running Linear Sorption")
            # Leistra et al., 2001
            mass_layer = cellarea() * depth * (theta_gas * model.k_h / model.k_d * conc_ads  # gas
                                                    + theta_aq_layer * conc_ads / model.k_d  # aqueous
                                                    + model.p_b * conc_ads)  # sorbed
        else:
            # Leistra et al., 2001
            mass_layer = cellarea() * depth * (theta_gas * model.k_h / model.k_d * conc_ads  # gas
                                                    + theta_aq_layer * conc_ads / model.k_d  # aqueous
                                                    + model.p_b * conc_ads)  # sorbed
    else:
        print("Running Linear Sorption, no Gas Phase")
        # Whelan, 1987 # No gas phase considered
        mass_layer = cellarea() * depth * (theta_aq_layer * conc_ads / model.k_d
                                                + model.p_b * conc_ads)  # mg
    return mass_layer  # mg




def getRsample(model, layer):
    if layer == 0:
        delta = model.delta_z0
    elif layer == 1:
        delta = model.delta_z1
    elif layer == 2:
        delta = model.delta_z2
    else:
        print("Inaccurate layer index")
        raise NotImplementedError

    r_sample = model.r_standard * (delta + 1)
    return r_sample


def get_species(model, r_sample, conc_aq_tot):
    delta13c = (r_sample - model.r_standard) / model.r_standard
    conc_aq_light = conc_aq_tot / (model.r_standard * (delta13c + 1) + 1)
    conc_aq_heavy = conc_aq_tot - conc_aq_light
    return {"light": conc_aq_light, "heavy": conc_aq_heavy}


def degrade(model, layer,
            theta_sat_z0z1, theta_sat_z2, theta_fcap, theta_wp,
            sor_deg_factor=1,
            units='m'):
    """
    :type model: instance of "self"
    :param layer: layer integer (i.e., 0, 1, or 2)
    :param theta_fcap: field capacity for respective layer
    :param theta_wp: wilting point for respective layer
    :param mass_layer_t: current pesticide mass of layer
    :param sor_deg_factor: degradation factor adjusting rate in sorbed phase
    :param units: string points to necessary units conversion

    This function implements Dairon's temperature and moisture dependence
    on 1st order degradation kinetics.
    """
    # If bulk density is in Kg/m3,
    # need to convert layer's depth to [m]
    if units == 'm':
        unit_f = 10 ** 3  # -> depth in m
    else:
        unit_f = "specify correct unit"

    # Get mass applied if top-soil
    if layer == 0:
        theta_aq_layer = model.theta_z0
        temp_layer = model.temp_z0_fin
        depth = model.z0 / unit_f
        mass_layer = model.pestmass_z0

    elif layer == 1:
        theta_aq_layer = model.theta_z1
        temp_layer = model.temp_z1_fin
        depth = model.z1 / unit_f
        mass_layer = model.pestmass_z1

    elif layer == 2:
        theta_aq_layer = model.theta_z2
        temp_layer = model.temp_z2_fin
        depth = model.z2 / unit_f
        mass_layer = model.pestmass_z2

    else:
        print("incorrect number of layers, raising error")
        raise NotImplementedError

    # Moisture factor in biodegradation
    # F_Theta_1
    theta_factor = ifthenelse(theta_aq_layer <= 0.5 * theta_wp, 0,
                              ifthenelse(theta_aq_layer <= theta_fcap,
                                         (((theta_aq_layer - 0.5 * theta_wp) / (
                                             theta_fcap - theta_wp)) ** model.beta_moisture),
                                         1))
    # Temperature factor in biodegradation
    # F_T_1
    temp_factor = ifthenelse(temp_layer <= 0, 0,
                             ifthenelse(temp_layer < 5,
                                        (temp_layer / 5) * exp(model.alpha_temperature) * (5 - model.temp_ref),
                                        exp(model.alpha_temperature * (5 - model.temp_ref))
                                        )
                             )
    # Half-life as a function of temperature and moisture
    dt_50 = max(model.dt_50_ref * theta_factor * temp_factor, 0)

    # Convert to degradation constant
    # Deg in dissolved phase
    k_b = ifthenelse(dt_50 > 0, ln(2) / dt_50, 0)  # Constant of degradation (-) is dynamic, based on Theta and Temp.
    # Deg in sorbed phase (now assumed equal)
    k_bs = k_b * sor_deg_factor

    # Step 0 - Obtain species concentration (aqueous phase)
    conc_layer_aq = get_conc_from_mass(model, layer,
                                       theta_sat_z0z1, theta_sat_z2,
                                       return_phase="aq", gas=True, isotopes=False)
    r_sample = getRsample(model, layer)
    fractions = get_species(model, r_sample, conc_layer_aq)
    conc_aq_light_ini = fractions["light"]
    conc_aq_heavy_ini = fractions["heavy"]

    # Step 1 - Degrade
    # First order degradation kinetics
    conc_aq_light_new = conc_aq_light_ini * exp(-k_b * model.jd_dt)
    conc_aq_heavy_new = conc_aq_heavy_ini * exp(- model.alpha_iso * k_b * model.jd_dt)

    # Step 2 - Re-equilibrating sorbed pesticide, based on new total mass (i.e., after aqueous degradation)
    mass_light_new = get_mass_from_aq(model, layer,
                                      theta_sat_z0z1, theta_sat_z2,
                                      conc_aq_light_new)  # mg
    mass_heavy_new = get_mass_from_aq(model, layer,
                                      theta_sat_z0z1, theta_sat_z2,
                                      conc_aq_heavy_new)  # mg
    conc_ads_light = get_conc_from_mass(model, layer,
                                        theta_sat_z0z1, theta_sat_z2,
                                        return_phase="ads", gas=True, isotopes=True, mass=mass_light_new)
    conc_ads_heavy = get_conc_from_mass(model, layer,
                                        theta_sat_z0z1, theta_sat_z2,
                                        return_phase="ads", gas=True, isotopes=True, mass=mass_heavy_new)

    # Step 3 - Degrade sorbed fraction (currently same rate as in aqueous)
    conc_ads_light_new = conc_ads_light * exp(-k_bs * model.jd_dt)
    conc_ads_heavy_new = conc_ads_heavy * exp(- model.alpha_iso * k_bs * model.jd_dt)

    # Step 4 - Re-equilibrating aqueous pesticide, based on new total mass (i.e., after sorbed degradation)
    mass_light_fin = get_mass_from_ads(model, layer, theta_sat_z0z1, theta_sat_z2, conc_ads_light_new)
    mass_heavy_fin = get_mass_from_ads(model, layer, theta_sat_z0z1, theta_sat_z2, conc_ads_heavy_new)
    conc_aq_light_fin = get_conc_from_mass(model, layer,
                                           theta_sat_z0z1, theta_sat_z2,
                                           return_phase="aq", gas=True, isotopes=True, mass=mass_light_fin)
    conc_aq_heavy_fin = get_conc_from_mass(model, layer,
                                           theta_sat_z0z1, theta_sat_z2,
                                           return_phase="aq", gas=True, isotopes=True, mass=mass_heavy_fin)

    # Change in mass due to aqueous & adsorbed degradation
    mass_light_fin = get_mass_from_aq(model, layer,
                                      theta_sat_z0z1, theta_sat_z2,
                                      conc_aq_light_fin)
    mass_heavy_fin = get_mass_from_aq(model, layer,
                                      theta_sat_z0z1, theta_sat_z2,
                                      conc_aq_heavy_fin)

    mass_deg_light = (mass_light_fin -
                      get_mass_from_aq(model, layer,
                                       theta_sat_z0z1, theta_sat_z2,
                                       conc_aq_light_ini))

    mass_deg_heavy = (mass_heavy_fin -
                      get_mass_from_aq(model, layer,
                                       theta_sat_z0z1, theta_sat_z2,
                                       conc_aq_heavy_ini))

    return {'mass_light_fin': mass_light_fin,
            'mass_heavy_fin': mass_heavy_fin,
            'mass_deg_light': mass_deg_light,
            'mass_deg_heavy': mass_deg_heavy
            }


def update_layer_delta(model, layer, process, mass_process,
                       mass_before_transport):
    if layer == 0:
        delta_layer = model.delta_z0
        delta_layer_above = None
        mass_layer = model.pestmass_z0
    elif layer == 1:
        delta_layer = model.delta_z1
        delta_layer_above = model.delta_z0
        mass_layer = model.pestmass_z1
    elif layer == 2:
        delta_layer = model.delta_z2
        delta_layer_above = model.delta_z1
        mass_layer = model.pestmass_z2

    if process == "volat":
        mass_loss = mass_process["mass_loss"]
        mass_gain = 0
        delta_gain = 0
        delta_loss = delta_layer
    elif process == "runoff":
        mass_loss = mass_process["mass_runoff"]
        mass_gain = 0
        delta_gain = 0
        delta_loss = delta_layer
    elif process == "leach":
        mass_leached = mass_process["mass_leached"]  # mg
        mass_loss = mass_leached
        mass_gain = 0
        delta_gain = delta_layer_above
        delta_loss = delta_layer
    elif process == "latflux":
        print('Latflux implemented in getLatMassDeltaFlux()')
        raise NotImplementedError
        # mass_latflux = mass_process["net_mass_latflux"]  # mg
        # mass_loss = mass_process["cell_mass_loss_downstream"]
        # mass_gain = mass_process["upstream_mass_inflow"]
        # mass_tot = mass_before_transport + mass_gain - mass_loss
        # # 1st fraction term (f1), must: f1 > 1, with net gain
        # f1 = (mass_before_transport * delta_layer) / mass_tot
        # # 3rd fraction (f3)
        # cell_massDelta_loss = max(c * (depth * theta_layer - depth * theta_fcap), scalar(0)) * conc_layer_aq * delta_layer
        # f3 = cell_massDelta_loss / mass_tot
        # # 2nd fraction (f2)
        # cell_massdelta_gain = (model.wetness * accuflux(model.ldd_subs, cell_massDelta_loss)) / accuflux(model.ldd_subs, model.wetness)
        # f2 = cell_massdelta_gain / mass_tot
    else:
        raise NotImplementedError

    if process == "latflux":
        pass  # delta_int = f1 + f2 - f3
        # Not convinced of this:
        # delta2_f2 = accuflux(model.ldd_subs, mass_gain*delta_layer)/accuflux(model.ldd_subs, mass_tot)
    else:
        mass_tot = mass_before_transport + mass_gain - mass_loss
        delta_int = ((1/mass_tot) *
                     (delta_layer * mass_before_transport +  # initial
                      delta_gain * mass_gain -  # mass_in
                      delta_loss * mass_loss))  # mass_out

    # return {"delta_int": delta_int, "mass_layer": mass_layer}
    return delta_int


def getLatMassDeltaFlux(model, layer, theta_sat, theta_fcap, mass_before_transport):
    if layer == 0:
        depth = model.z0
        delta_layer = model.delta_z0
        theta_layer = model.theta_z0
        mass_layer = model.pestmass_z0
        c = model.c1
    elif layer == 1:
        depth = model.z1
        delta_layer = model.delta_z1
        theta_layer = model.theta_z1
        mass_layer = model.pestmass_z1
        c = model.c1
    elif layer == 2:
        store = False
        depth = model.z2
        delta_layer = model.delta_z2
        theta_layer = model.theta_z2
        mass_layer = model.pestmass_z2
        c = model.c2

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        retard_layer = 1

    if gas:
        # Leistra et al., 2001
        theta_gas = theta_sat - theta_layer
        conc_layer_aq = mass_layer / ((cellarea() * depth) *
                                      (theta_gas * model.k_h + theta_layer * retard_layer))  # mg/L
    else:
        # Whelan, 1987 # No gas phase considered
        conc_layer_aq = (mass_layer / cellarea()) / (theta_layer * retard_layer * depth)  # mg/L

    # W(j/i)
    rel_wetness = model.wetness/accuflux(model.ldd_subs, model.wetness)

    # Cell mass loss/gain (to update only mass)
    mass_loss = max(conc_layer_aq * (c * (depth * theta_layer - depth * theta_fcap)), scalar(0))
    mass_gain = rel_wetness * accuflux(model.ldd_subs, cell_mass_loss_downstream)
    net_mass_latflux = mass_gain - mass_loss

    # massDelta  loss/gain (to update only delta)
    massdC_loss = max(delta_layer * conc_layer_aq * (c * (depth * theta_layer - depth * theta_fcap)), scalar(0))
    massdC_gain = rel_wetness * accuflux(model.ldd_subs, massdC_loss)

    # Isotope mass balance update three terms:
    mass_tot = mass_before_transport + mass_gain - mass_loss
    # 1st fraction term (f1), must: f1 > 1, with net gain
    f1 = (mass_before_transport / mass_tot) * delta_layer
    f2 = massdC_gain / mass_tot
    f3 = massdC_loss / mass_tot
    delta_layer = f1 + f2 - f3

    return {
        'mass_loss': mass_loss,  # mg
        'mass_gain': mass_gain,  # mg
        'net_mass_latflux': net_mass_latflux,  # mg
        'massdC_loss': massdC_loss,  # mg*dC
        'massdC_gain': massdC_gain,  # mg*dC
        'f1': f1,  # could be > 1  (if net gain)
        'f2': f2,
        'f3': f3,
        'delta_layer': delta_layer
    }


def getVolatileMass(model, app_days, temp_air, theta_sat, rel_diff_model="option-1", sorption_model="linear", gas=True):
    # Volatilize only during peak volatilization time i.e., first 24 hrs, @Prueger2005.
    if model.jd_cum in app_days:
        theta_layer = model.theta_z0
        mass_layer = model.pestmass_z0  # ug
        depth_m = model.z0 * 1/10**3  # Convert to m (needed for final mass computation on cell basis)
        # Air boundary layer, assumed equivalent to top soil thickness
        thickness_a = depth_m  # m (Thickness air boundary layer)
        # D_ar (metolachlor) = 0.03609052694; Diffusion coefficient in air (cm^2/s); https://www.gsi-net.com
        diff_ar = 0.03609 * 86400 * 1/10**4  # m2/d (Diff. coeff in air at reference Temp., in Kelvin, D_a,r)
        # TODO: define actual temperature, the value here is temp_bare_soil, not temp_air
        diff_a = (temp_air / 293.15) ** 1.75 * diff_ar  # m2/d (Diff. coefficient adjusted to air Temp., D_a)

        if rel_diff_model == "option-1":
            # Millington and Quirk, 1960 (in Leistra, 2001, p.48)
            # a,b parameters: Jin and Jury, 1996 (in Leistra, 2001)
            diff_relative_gas = diff_a * (theta_sat - theta_layer) ** 2 / theta_sat ** (2 / 3)  # m2/d
        elif rel_diff_model == "option-2":
            # Currie 1960 (in Leistra, 2001)
            # a,b parameters: Baker, 1987 (in Leistra, 2001)
            diff_relative_gas = diff_a * 2.5 * (theta_sat - theta_layer) ** 3  # m2/d
        else:
            print("No appropriate relative diffusion parameter option chosen for getVolatileMass()")
            diff_relative_gas = diff_a  # m2/d

        r_a = thickness_a / diff_a  # d/m resistance for transport through boundary air layer
        r_s = (0.5 * depth_m) / diff_relative_gas  # d/m

        if sorption_model == "linear":
            # Retardation factor
            retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
        else:
            retard_layer = 1
        if gas:
            # Leistra et al., 2001
            theta_gas = theta_sat - theta_layer
            # Remember that: m2 * mm = L
            conc_layer_aq = mass_layer / ((cellarea() * model.z0) *
                                          (theta_gas * model.k_h + theta_layer * retard_layer))  # ug/L
        else:
            # Whelan, 1987 # No gas phase considered
            conc_layer_aq = (mass_layer / cellarea()) / (theta_layer * retard_layer * model.z0)  # ug/L

        # Convert ug/L to ug/m3, as will be multiplying by cell's area in m2
        conc_layer_aq *= 10 ** 3  # ug/m3
        conc_gas_layer = conc_layer_aq / model.k_h  # Henry's  ug/m3
        volat_flux = (conc_gas_layer / (r_a + r_s)) * cellarea()  # ug/day
    else:
        volat_flux = 0
    return {"mass_loss": volat_flux}  # ug / day


def getKfilm(model, runoffvelocity):
    """
    Note: Model uses run-off (mm) per day (i.e. timestep) as runoff velocity.
    Chemical parameter source:
    http://www.gsi-net.com/en/publications/gsi-chemical-database/single/377.html
    """
    # Dynamic viscosity of water (\mu) @25 Celsius = 8.9e-04 [Pa s]
    #   1 Pa = 1 N/(m s^2) = 1 Kg/(m s^2)
    #   Convert to g/(cm s): dyn_visc = 8.9e-03 [g/cm s]
    dyn_visc = 8.9e-03  # [g/cm s] @25 degrees, [@Shi2011]:\mu

    # Solute diffusivity in water (D_w)
    # Metolachlor = 5.0967719112e-006 (cm2 / s)
    diff_solute = 5.0967719112e-006  # [cm2 / s], [@Shi2011]:D_w
    Sc = dyn_visc / (model.p_b * diff_solute)  # (-) Schmidt number, [@Shi2011]:S_c

    # Reynolds number (dimensionless), 86400s = 1 day
    cell_length = 2 * 10 ** 3  # mm [@Shi2011]:L
    re = (model.p_b * 1 / 10 ** 2 * runoffvelocity * cell_length) / (dyn_visc * 86400)  # [-]  [@Shi2011]:Re
    kl = 0.664 * ((diff_solute * 86400 * 10 ** 2) / cell_length) * re ** (float(1) / 2) * Sc ** (float(1) / 3)  # mm/day
    return kl  # mm/day
    # NOT USED:
    # Kinematic water viscosity at 20 deg. Celsius = 1.004 x 10-6 (m2/s)
    # kin_visc = 1.004e-6 * 86400 * 10**6  # mm2/day


def getRunOffMass(model, theta_sat, precip, runoff_mm,
                  transfer_model="simple-mt", sorption_model="linear",
                  gas=True):
    depth, theta_layer, delta_layer = model.z0, model.theta_z0, model.delta_z0

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        retard_layer = 1

    # Aqueous concentration calculation
    if gas:
        # Leistra et al., 2001
        theta_gas = theta_sat - model.theta_z0
        # TODO: Check that all theta_gas >= 0
        conc_layer_aq = model.pestmass_z0 / ((cellarea() * depth) *
                                             (theta_gas * model.k_h + theta_layer * retard_layer))  # mg/L
    else:
        # Whelan, 1987 # No gas phase considered
        conc_layer_aq = (model.pestmass_z0 / cellarea()) / (theta_layer * retard_layer * depth)  # mg/L

    if transfer_model == "simple-mt":
        mass_ro = conc_layer_aq * runoff_mm * cellarea()  # mg
        deltaMass_ro = mass_ro * delta_layer
    elif transfer_model == "nu-mlm-ro":
        # non-uniform-mixing-layer-model-runoff (nu-mlm-ro)
        # Considers a decrease in effective transfer as mixing layer depth increases
        # Adapted from Ahuja and Lehman, 1983 in @Shi2011,
        # Adaptation replaces Precip by Runoff amount.
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (runoff_mm * cellarea()) * exp(-b * model.z0) * conc_layer_aq  # mg
        deltaMass_ro = mass_ro * delta_layer
    elif transfer_model == "nu-mlm":
        # non-uniform-mixing-layer-model (nu-mlm)
        # Original from Ahuja and Lehman, 1983 in @Shi2011
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (precip * cellarea()) * exp(-b * model.z0) * conc_layer_aq  # mg
        deltaMass_ro = mass_ro * delta_layer
    elif transfer_model == "d-mlm":
        # distributed mixing-layer-model (d-mlm)
        # Adapted from Havis et al., 1992, and
        # taking the K_L definition for laminar flow from Bennett and Myers, 1982.
        mass_ro = getKfilm(model, runoff_mm) * cellarea() * conc_layer_aq  # mg
        deltaMass_ro = mass_ro * delta_layer
    else:
        print("Run-off transfer model not stated")
        return None
    return {"mass_runoff": mass_ro, "deltaMass_runoff": deltaMass_ro}


def getLeachedMass(model, layer, theta_sat,
                   precip,
                   tot_percolation,
                   theta_after_percolate,
                   sorption_model=None,
                   leach_model=None,
                   gas=True):
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
        mass_layer = model.pestmass_z0  # mg
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
        mass_layer = model.pestmass_z1
    elif layer == 2:
        store = True
        depth = model.z2
        theta_layer = model.theta_z2
        mass_layer = model.pestmass_z2

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        print("No sorption assumed, Ret. factor = 1")
        retard_layer = 1  # No retardation.

    if gas:
        # Leistra et al., 2001
        theta_gas = theta_sat - theta_layer
        conc_layer_aq = mass_layer / ((cellarea() * depth) *
                                      (theta_gas * model.k_h + theta_layer * retard_layer))  # mg/L
    else:
        # Whelan, 1987 # No gas phase considered
        conc_layer_aq = (mass_layer / cellarea()) / (theta_layer * retard_layer * depth)  # mg/L

    if leach_model == "mcgrath":
        if layer == 0:
            conc_layer_new_aq = conc_layer_aq * exp(-precip / theta_layer * retard_layer * depth)
            mass_layer = conc_layer_aq * (theta_layer * depth) * cellarea()  # mg
            mass_layer_new = conc_layer_new_aq * (theta_after_percolate * depth) * cellarea()  # mg
            mass_leached = mass_layer - mass_layer_new  # mg
        else:
            # McGrath not used in lower layers,
            # as formulation accounts for rainfall impact
            mass_leached = conc_layer_aq * tot_percolation * cellarea()  # mg
            mass_layer_new = mass_layer - mass_leached
    else:
        mass_leached = conc_layer_aq * tot_percolation * cellarea()  # mg
        mass_layer_new = mass_layer - mass_leached

    return {"mass_leached": mass_leached}  # mg


def getLateralMassFlux(model, layer, theta_sat, theta_fcap, sorption_model="linear", gas=True):
    """
    This function uses initial theta to compute mass flux, therefore model.theta_zX
    should not be updated before running this function.
    I.e., update layer's moisture and mass only after running
    both getLayerMoisture(...) and getLateralMassFlux(...)
    Alternatively, call "theta_ini" from the result of getLayerMoisture(...)
    Note: Both of these functions could be integrated to avoid double calculation?
    :param model:
    :param layer:
    :param theta_fcap:
    :param conc_layer:
    :return:
    """
    if layer == 0:
        depth = model.z0
        delta_layer = model.delta_z0
        theta_layer = model.theta_z0
        mass_layer = model.pestmass_z0
        c = model.c1
    elif layer == 1:
        depth = model.z1
        delta_layer = model.delta_z1
        theta_layer = model.theta_z1
        mass_layer = model.pestmass_z1
        c = model.c1
    elif layer == 2:
        store = False
        depth = model.z2
        delta_layer = model.delta_z2
        theta_layer = model.theta_z2
        mass_layer = model.pestmass_z2
        c = model.c2
    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        retard_layer = 1

    if gas:
        # Leistra et al., 2001
        theta_gas = theta_sat - theta_layer
        conc_layer_aq = mass_layer / ((cellarea() * depth) *
                                      (theta_gas * model.k_h + theta_layer * retard_layer))  # mg/L
    else:
        # Whelan, 1987 # No gas phase considered
        conc_layer_aq = (mass_layer / cellarea()) / (theta_layer * retard_layer * depth)  # mg/L

    # Cell mass loss
    cell_mass_loss_downstream = max(c * (depth * theta_layer - depth * theta_fcap), scalar(0)) * conc_layer_aq
    upstream_mass_inflow = (model.wetness * accuflux(model.ldd_subs, cell_mass_loss_downstream)) / accuflux(model.ldd_subs,
                                                                                                       model.wetness)
    net_mass_latflux = upstream_mass_inflow - cell_mass_loss_downstream
    deltaMass_latflux = delta_layer*net_mass_latflux
    return {"net_mass_latflux": net_mass_latflux,
            "deltaMass_latflux": deltaMass_latflux,
            "upstream_mass_inflow": upstream_mass_inflow,
            "cell_mass_loss_downstream": cell_mass_loss_downstream}  # mg

