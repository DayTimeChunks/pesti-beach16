# -*- coding: utf-8 -*-

from pcraster._pcraster import *
from pcraster.framework import *

# import os
# import time

DEBUG_pest = False
global DEBUG_pest


def getConcAq(model, layer, theta_sat, mass,
              sorption_model="linear", gas=True):
    # Note that p_b (g/cm3) x k_d (L/Kg) -> unit-less
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
    elif layer == 2:
        depth = model.z2
        theta_layer = model.theta_z2

    if sorption_model == "linear":
        # Retardation factor (dimensionless)
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        print("No sorption assumed, Ret. factor = 1")
        retard_layer = 1  # No retardation.

    if gas:  # Leistra et al., 2001
        theta_gas = max(theta_sat - theta_layer, scalar(0))
        conc_aq = mass / ((cellarea() * depth) *  # m2 * mm = L
                          (theta_gas / model.k_h +
                           theta_layer * retard_layer))  # ug/L cell volume
    else:  # No gas phase
        # Whelan, 1987
        conc_aq = (mass / (cellarea() * depth * theta_layer * retard_layer))

    return conc_aq


def getConcAds(model, layer, theta_sat, mass, gas=True):
    # mass / Kg soil
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
    elif layer == 2:
        depth = model.z2
        theta_layer = model.theta_z2

    if gas:
        theta_gas = max(theta_sat - theta_layer, scalar(0))
        # [mass pest/Kg soil]
        conc_ads = mass / ((cellarea() * depth) *
                           (theta_gas / (model.k_h * model.k_d) +
                            theta_layer / model.k_d +
                            model.p_b))
    else:
        print("No implementation without gas available")
        raise NotImplementedError

    return conc_ads


# def getMassAq(model, layer, theta_sat,
#               conc_aq, gas=True, sorption_model="linear"):
#     """
#     Used to update mass due to equilibrium conditions after degradation
#     :param model:
#     :param layer:
#     :param conc_aq:
#     :param gas:
#     :param sorption_model:
#     :return:
#     """
#     if layer == 0:
#         theta_aq_layer = model.theta_z0
#         # theta_sat = theta_sat_z0z1
#         depth = model.z0  # mm
#     elif layer == 1:
#         theta_aq_layer = model.theta_z1
#         # theta_sat = theta_sat_z0z1
#         depth = model.z1  # mm
#     elif layer == 2:
#         theta_aq_layer = model.theta_z2
#         # theta_sat = theta_sat_z2
#         depth = model.z2  # mm
#     else:
#         print("incorrect number of layers, raising error")
#         raise NotImplementedError
#
#     if gas:
#         if sorption_model == "freundlich":
#             print("Freundlich not implemented, running Linear Sorption")
#             # Leistra et al., 2001
#             theta_gas = max(theta_sat - theta_aq_layer, scalar(0))
#             mass_layer = cellarea() * depth * ((theta_gas / model.k_h) * conc_aq  # gas
#                                                + theta_aq_layer * conc_aq  # aqueous
#                                                + model.p_b * model.k_d * conc_aq)  # sorbed
#         else:
#             # Leistra et al., 2001
#             theta_gas = max(theta_sat - theta_aq_layer, scalar(0))
#             mass_layer = cellarea() * depth * ((theta_gas / model.k_h) * conc_aq  # gas
#                                                + theta_aq_layer * conc_aq  # aqueous
#                                                + model.p_b * model.k_d * conc_aq)  # sorbed
#     else:
#         print("Running Linear Sorption, no Gas Phase")
#         # Whelan, 1987 # No gas phase considered
#         mass_layer = cellarea() * depth * (theta_aq_layer * conc_aq
#                                            + model.p_b * model.k_d * conc_aq)
#     return mass_layer


def getVolatileMass(model, temp_air, theta_sat, mass, frac,
                    rel_diff_model="option-1", sorption_model="linear",
                    gas=True, isotopes=True, ):
    # Volatilize only during peak volatilization time i.e., first 24 hrs, @Prueger2005.
    theta_layer = model.theta_z0
    theta_gas = max(theta_sat - theta_layer, scalar(0))
    # Convert to m (needed for final mass computation on cell basis)
    depth_m = model.z0 * 1 / 10 ** 3
    # Air boundary layer, assumed as 2m high
    thickness_a = scalar(1.0)  # m
    # Diffusion coefficient in air (cm^2/s); https://www.gsi-net.com
    #  D_ar (metolachlor) = 0.03609052694,  at reference Temp., in Kelvin, D_a,r)
    diff_ar = 0.03609052694 * 86400.0 * 1.0 / 10 ** 4  # m2/d
    # Diffusion coefficient adjusted to air Temp. in Kelvin, D_a
    diff_a = ((temp_air + 273.15) / 293.15) ** 1.75 * diff_ar  # m2/d

    if rel_diff_model == "option-1":
        # Millington and Quirk, 1960 (in Leistra, 2001, p.48)
        # a,b parameters: Jin and Jury, 1996 (in Leistra, 2001)
        diff_relative_gas = (diff_a * theta_gas ** 2 /
                             theta_sat ** (2 / 3))  # m2/d
    elif rel_diff_model == "option-2":
        # Currie 1960 (in Leistra, 2001)
        # a,b parameters: Baker, 1987 (in Leistra, 2001)
        diff_relative_gas = diff_a * 2.5 * theta_gas ** 3  # m2/d
    else:
        print("No appropriate relative diffusion parameter chosen")
        diff_relative_gas = diff_a  # m2/d
    # Transport resistance through air (r_a) and soil (r_s) layer
    r_a = thickness_a / diff_a  # d/m
    r_s = (0.5 * depth_m) / diff_relative_gas  # d/m

    conc_aq = getConcAq(model, 0, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)
    # Convert ug/L to ug/m3, as will be multiplying by cell's area in m2
    conc_aq *= 10 ** 3  # ug/L * 10^3 L/m3
    conc_gas = conc_aq / model.k_h  # ug/L air
    volat_flux = (conc_gas / (r_a + r_s)) * cellarea()  # ug/day
    return volat_flux


# Runoff
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
    cell_length = 2 * 10 ** 3  # mm
    # Reynolds (Re), [-] (Shi et al., 2011)
    re = (model.p_b * 1 / 10 ** 2 * runoffvelocity * cell_length) / (dyn_visc * 86400)
    kl = (0.664 * ((diff_solute * 86400 * 10 ** 2) / cell_length) *
          re ** (float(1) / 2) * Sc ** (float(1) / 3))
    return kl  # mm/day


def getRunOffMass(model, theta_sat, precip, runoff_mm,
                  mass,
                  transfer_model="simple-mt", sorption_model="linear",
                  gas=True, debug=False):
    # Aqueous concentration
    conc_aq = getConcAq(model, 0, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)

    if transfer_model == "simple-mt":
        mass_ro = conc_aq * runoff_mm * cellarea()
    elif transfer_model == "nu-mlm-ro":
        # non-uniform-mixing-layer-model-runoff (nu-mlm-ro)
        # Considers a decrease in effective transfer as mixing layer depth increases
        # Adapted from Ahuja and Lehman, 1983 in @Shi2011,
        # Adaptation replaces Precip by Runoff amount.
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (runoff_mm * cellarea()) * exp(-b * model.z0) * conc_aq
    elif transfer_model == "nu-mlm":
        # non-uniform-mixing-layer-model (nu-mlm)
        # Original from Ahuja and Lehman, 1983 in @Shi2011
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (precip * cellarea()) * exp(-b * model.z0) * conc_aq
    elif transfer_model == "d-mlm":
        # distributed mixing-layer-model (d-mlm)
        # Adapted from Havis et al., 1992, and
        # taking the K_L definition for laminar flow from Bennett and Myers, 1982.
        mass_ro = getKfilm(model, runoff_mm) * cellarea() * conc_aq
    else:
        print("Run-off transfer model not stated")
        return None

    if debug:
        pass
        # model.report(conc_aq, 'aCo')
        # model.report(mass_ro, 'aMROa')
        # model.report(runoff_mm, 'aRO')
    return mass_ro


def getLeachedMass(model, layer, theta_sat, theta_fc,
                   water_flux,
                   theta_after_percolate,
                   mass,
                   sorption_model=None,
                   leach_model=None, gas=True, debug=False):
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
    elif layer == 2:
        depth = model.z2
        theta_layer = model.theta_z2

    # Aqueous concentration
    conc_aq = getConcAq(model, layer, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = scalar(1) + (model.p_b * model.k_d) / theta_layer
    else:
        print("No sorption assumed, Ret. factor = 1")
        retard_layer = scalar(1)  # No retardation.

    if leach_model == "mcgrath":
        if layer == 0:
            conc_aq_new = conc_aq * exp(-water_flux / (theta_layer * retard_layer * depth))
            mass_aq = conc_aq * (theta_layer * depth * cellarea())
            mass_aq_new = conc_aq_new * (theta_after_percolate * depth * cellarea())
            mass_leached = mass_aq - mass_aq_new
            if mapminimum(mass_aq_new) < 0:
                print("Error in Leached Model, layer: ", str(layer))
                model.report(mass_leached, 'aZ' + str(layer) + 'LCH')
        else:
            # McGrath not used in lower layers,
            # as formulation accounts for rainfall impact
            max_flux = max(min(water_flux, (theta_layer - theta_fc) * depth), scalar(0))
            mass_leached = conc_aq * max_flux * cellarea()
            mass_aq = conc_aq * (theta_layer * depth * cellarea())
            mass_aq_new = mass_aq - mass_leached
            if mapminimum(mass_aq_new) < -1e10-6:
                print("Error in Leached Model, layer: ", str(layer))
                model.report(mass_leached, 'aZ' + str(layer) + 'LCH')
            if mapminimum(mass_aq_new) < 0:
                mass_leached = max(mass_leached, scalar(0))
                # mass_aq_new = mass_aq - mass_leached
    else:
        mass_leached = conc_aq * water_flux * cellarea()
        mass_aq = conc_aq * (theta_layer * depth) * cellarea()
        mass_aq_new = mass_aq - mass_leached
        if mapminimum(mass_aq_new) < 0:
            print("Error in Leached Model")

    return mass_leached


def getLatMassFlux(model, layer, theta_sat, theta_fcap,
                   mass, sorption_model='linear', gas=True, debug=False):
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
        c = model.c1
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
        c = model.c1
    elif layer == 2:
        depth = model.z2
        theta_layer = model.theta_z2
        c = model.c2

    # Aqueous concentration
    conc_aq = getConcAq(model, layer, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)

    # W(j/i)
    rel_wetness = model.wetness / accuflux(model.ldd_subs, model.wetness)

    # Cell mass loss/gain (to update only mass)
    mass_loss = max(conc_aq * (c * (depth * theta_layer - depth * theta_fcap)), scalar(0))
    mass_gain = rel_wetness * accuflux(model.ldd_subs, mass_loss)
    net_mass_latflux = mass_gain - mass_loss

    return {
        'mass_loss': mass_loss,
        'mass_gain': mass_gain,
        'net_mass_latflux': net_mass_latflux
    }


def getDrainMassFlux(model, layer, theta_sat, theta_fcap, mass,
                     sorption_model='linear', gas=True, debug=False):
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
        c = model.c_adr
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
        c = model.c_adr
    elif layer == 2:
        depth = model.z2
        theta_layer = model.theta_z2
        c = model.c_adr

    # TODO: Conc. and mass loss are being computed with initial theta
    # but water loss potential accounts for intermediate updates to water content.
    # Aqueous concentration
    conc_aq = getConcAq(model, layer, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)
    # TODO: water loss here is not equivalent to water loss in hydro.py module for lateral flow
    # due to intermediate updates to water content
    # -> Thus likely leads to overestimation of mass transport.
    mass_loss = max(conc_aq * (c * (depth * theta_layer - depth * theta_fcap)), scalar(0))
    return mass_loss


def getMassDegradation(model, layer,
                       theta_sat, theta_fcap, theta_wp,
                       mass, frac="L", sor_deg_factor=1,
                       sorption_model="linear", gas=True):
    # Get mass applied if top-soil
    if layer == 0:
        theta_aq_layer = model.theta_z0
        temp_layer = model.temp_z0_fin
        depth = model.z0
    elif layer == 1:
        theta_aq_layer = model.theta_z1
        temp_layer = model.temp_z1_fin
        depth = model.z1
    elif layer == 2:
        theta_aq_layer = model.theta_z2
        temp_layer = model.temp_z2_fin
        depth = model.z2

    # F_Theta_1
    theta_factor = ifthenelse(theta_aq_layer <= 0.5 * theta_wp, scalar(0),
                              ifthenelse(theta_aq_layer <= theta_fcap,
                                         (((theta_aq_layer - 0.5 * theta_wp) / (
                                             theta_fcap - theta_wp)) ** scalar(model.beta_moisture)),
                                         scalar(1)))
    # Temperature factor in biodegradation
    # F_T_1
    # temp_factor = max(ifthenelse(temp_layer <= scalar(0), scalar(0),
    #                              ifthenelse(temp_layer < 5,
    #                                         (temp_layer / 5) * exp(model.alpha_temperature) * (5 - model.temp_ref),
    #                                         exp(model.alpha_temperature * (5 - model.temp_ref))
    #                                         )
    #                              ), scalar(0))
    temp_factor = exp((model.act_e / model.r_gas) * (1 / (model.temp_ref + 273.15) - 1 / (model.temp_air + 273.15)))

    # Half-life as a function of temperature and moisture
    dt_50 = max(model.dt_50_ref * theta_factor * temp_factor, scalar(0))
    # dt_50 = max(model.dt_50_ref)

    # Convert to degradation constant
    # Deg in dissolved phase
    k_b = ifthenelse(dt_50 > 0, ln(2) / dt_50,
                     scalar(0))  # Constant of degradation (-) is dynamic, based on Theta and Temp.
    # Deg in sorbed phase (now assumed equal)
    k_bs = k_b * sor_deg_factor

    # Step 0 - Obtain species concentration (all phases)
    conc_aq = getConcAq(model, layer, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas)  # mass/L
    conc_ads = getConcAds(model, layer, theta_sat, mass, gas=gas)  # mass/g soil
    conc_gas = conc_aq / model.k_h

    mass_aq = conc_aq * (theta_aq_layer * depth * cellarea())
    mass_ads = conc_ads * (depth * cellarea() * model.p_b)  # pb = g/cm3

    # Step 1 - Degrade phase fractions
    # First order degradation kinetics
    theta_gas = max(theta_sat - theta_aq_layer, scalar(0))
    if frac == "H":
        conc_aq_new = conc_aq * exp(-1 * model.alpha_iso * k_b * scalar(model.jd_dt))
        conc_ads_new = conc_ads * exp(-1 * model.alpha_iso * k_bs * scalar(model.jd_dt))
    else:  # Same for total conc as for "L"
        conc_aq_new = conc_aq * exp(-1 * k_b * scalar(model.jd_dt))
        conc_ads_new = conc_ads * exp(-1 * k_bs * scalar(model.jd_dt))

    # Step 2 - Convert to mass (i.e., after degradation in each phase)
    mass_aq_new = conc_aq_new * (theta_aq_layer * depth * cellarea())
    mass_ads_new = conc_ads_new * (model.p_b * depth * cellarea())  # pb = g/cm3
    mass_gas = conc_gas * (theta_gas * depth * cellarea())
    mass_tot_new = mass_aq_new + mass_ads_new + mass_gas
    mass_deg_aq = mass_aq - mass_aq_new
    mass_deg_ads = mass_ads - mass_ads_new
    # if frac == "L" and layer == 0:
    #     model.report(mass, 'adt0M')
    #     model.report(mass_tot_new, 'adt1M')
    #   # model.report(mass_aq_new, 'adt1M')
    return {"mass_tot_new": mass_tot_new,
            "mass_deg_aq": mass_deg_aq,
            "mass_deg_ads": mass_deg_ads}
