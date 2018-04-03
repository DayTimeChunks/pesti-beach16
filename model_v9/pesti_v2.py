# -*- coding: utf-8 -*-

from pcraster._pcraster import *
from pcraster.framework import *

# import os
# import time

DEBUG_pest = False
global DEBUG_pest


def getConcAq(model, layer, theta_sat, mass,
              sorption_model="linear", gas=True, isotopes=True):
    if layer == 0:
        depth = model.z0
        theta_layer = model.theta_z0
    elif layer == 1:
        depth = model.z1
        theta_layer = model.theta_z1
    elif layer == 2:
        store = False
        depth = model.z2
        theta_layer = model.theta_z2

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
    else:
        print("No sorption assumed, Ret. factor = 1")
        retard_layer = 1  # No retardation.

    if gas:  # Leistra et al., 2001
        theta_gas = max(theta_sat - theta_layer, scalar(0))
        conc_aq = mass / ((cellarea() * depth) *  # m2 * mm = L
                          (theta_gas / model.k_h +
                           theta_layer * retard_layer))  # ug/L
    else:  # No gas phase
        # Whelan, 1987
        conc_aq = (mass / (cellarea() * depth * theta_layer * retard_layer))  # ug/L

    return conc_aq


def getVolatileMass(model, temp_air, theta_sat, mass, frac,
                    rel_diff_model="option-1", sorption_model="linear",
                    gas=True, isotopes=True, ):
    # Volatilize only during peak volatilization time i.e., first 24 hrs, @Prueger2005.
    theta_layer = model.theta_z0
    theta_gas = max(theta_sat - theta_layer, scalar(0))
    # Convert to m (needed for final mass computation on cell basis)
    depth_m = model.z0 * 1 / 10 ** 3
    # Air boundary layer, assumed as 2m high
    thickness_a = scalar(2.0)  # m
    # Diffusion coefficient in air (cm^2/s); https://www.gsi-net.com
    #  D_ar (metolachlor) = 0.03609052694,  at reference Temp., in Kelvin, D_a,r)
    diff_ar = 0.03609 * 86400 * 1 / 10 ** 4  # m2/d
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
                        sorption_model=sorption_model, gas=gas, isotopes=isotopes)
    # Convert ug/L to ug/m3, as will be multiplying by cell's area in m2
    conc_aq *= 10 ** 3  # ug/L * 10^3 L/m3
    volat_flux = (conc_aq / (r_a + r_s)) * cellarea()  # ug/day
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
    cell_length = 2 * 10 ** 3  # mm [@Shi2011]:L
    re = (model.p_b * 1 / 10 ** 2 * runoffvelocity * cell_length) / (dyn_visc * 86400)  # [-]  [@Shi2011]:Re
    kl = 0.664 * ((diff_solute * 86400 * 10 ** 2) / cell_length) * re ** (float(1) / 2) * Sc ** (float(1) / 3)  # mm/day
    return kl  # mm/day


def getRunOffMass(model, theta_sat, precip, runoff_mm,
                  mass, frac,
                  transfer_model="simple-mt", sorption_model="linear",
                  gas=True, isotopes=True):
    # Aqueous concentration
    conc_aq = getConcAq(model, 0, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas, isotopes=isotopes)

    if transfer_model == "simple-mt":
        mass_ro = conc_aq * runoff_mm * cellarea()  # ug
    elif transfer_model == "nu-mlm-ro":
        # non-uniform-mixing-layer-model-runoff (nu-mlm-ro)
        # Considers a decrease in effective transfer as mixing layer depth increases
        # Adapted from Ahuja and Lehman, 1983 in @Shi2011,
        # Adaptation replaces Precip by Runoff amount.
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (runoff_mm * cellarea()) * exp(-b * model.z0) * conc_aq  # ug
    elif transfer_model == "nu-mlm":
        # non-uniform-mixing-layer-model (nu-mlm)
        # Original from Ahuja and Lehman, 1983 in @Shi2011
        b = 1  # [mm] Calibration constant, 1 >= b > 0 (b-ranges appear reasonable).
        # As b decreases, mass transfer increases, model.z0 in mm
        mass_ro = (precip * cellarea()) * exp(-b * model.z0) * conc_aq  # ug
    elif transfer_model == "d-mlm":
        # distributed mixing-layer-model (d-mlm)
        # Adapted from Havis et al., 1992, and
        # taking the K_L definition for laminar flow from Bennett and Myers, 1982.
        mass_ro = getKfilm(model, runoff_mm) * cellarea() * conc_aq  # ug
    else:
        print("Run-off transfer model not stated")
        return None

    if frac == 'LL':
        model.report(conc_aq, 'aCo')
        # model.report(mass_ro, 'aMROa')
        # model.report(runoff_mm, 'aRO')
    return mass_ro


def getLeachedMass(model, layer, theta_sat,
                   water_flux,
                   theta_after_percolate,
                   mass,
                   sorption_model=None,
                   leach_model=None, gas=True, isotopes=True):
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
    conc_aq = getConcAq(model, 0, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas, isotopes=isotopes)

    if sorption_model == "linear":
        # Retardation factor
        retard_layer = scalar(1) + (model.p_b * model.k_d) / theta_layer
    else:
        print("No sorption assumed, Ret. factor = 1")
        retard_layer = scalar(1)  # No retardation.

    if leach_model == "mcgrath":
        if layer == 0:
            conc_aq_new = conc_aq * exp(-water_flux / theta_layer * retard_layer * depth)
            mass_layer = conc_aq * (theta_layer * depth) * cellarea()  # mg
            mass_layer_new = conc_aq_new * (theta_after_percolate * depth) * cellarea()  # mg
            mass_leached = max(mass_layer - mass_layer_new, scalar(0))  # mg
        else:
            # McGrath not used in lower layers,
            # as formulation accounts for rainfall impact
            mass_leached = conc_aq * water_flux * cellarea()  # mg
            mass_new = mass - mass_leached
    else:
        mass_leached = conc_aq * water_flux * cellarea()  # mg
        mass_new = mass - mass_leached

    return mass_leached


def getLatMassFlux(model, layer, theta_sat, theta_fcap,
                   mass, sorption_model='linear', gas=True, isotopes=True):
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
    conc_aq = getConcAq(model, 0, theta_sat, mass,
                        sorption_model=sorption_model, gas=gas, isotopes=isotopes)

    # W(j/i)
    rel_wetness = model.wetness / accuflux(model.ldd_subs, model.wetness)

    # Cell mass loss/gain (to update only mass)
    mass_loss = max(conc_aq * (c * (depth * theta_layer - depth * theta_fcap)), scalar(0))
    mass_gain = rel_wetness * accuflux(model.ldd_subs, mass_loss)
    net_mass_latflux = mass_gain - mass_loss

    return {
        'mass_loss': mass_loss,  # ug
        'mass_gain': mass_gain,  # ug
        'net_mass_latflux': net_mass_latflux  # ug
    }

