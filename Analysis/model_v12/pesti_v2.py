# -*- coding: utf-8 -*-
from pcraster._pcraster import *
from pcraster.framework import *

# import os
# import time
from copy import deepcopy

DEBUG_pest = False
global DEBUG_pest


def getConcAq(model, layer, mass, sorption_model="linear", gas=True):
    # Note that p_b (g/cm3) x k_d (L/Kg) -> unit-less
    theta_layer = model.theta[layer]
    depth = model.layer_depth[layer]

    if mapminimum(model.theta[layer]) < scalar(1e-06):
        conc_aq = scalar(0)
    else:
        if sorption_model == "linear":
            # Retardation factor (dimensionless)
            retard_layer = 1 + (model.p_b * model.k_d) / theta_layer
        else:
            print("No sorption assumed, Ret. factor = 2")
            retard_layer = 1  # No retardation.

        if gas:  # Leistra et al., 2001
            theta_gas = max(model.theta_sat[layer] - theta_layer, scalar(0))

            conc_aq = max(scalar(0), mass / ((cellarea() * depth) *  # m2 * mm = L
                              (theta_gas / model.k_h +
                               theta_layer * retard_layer)))  # mass/L cell volume

            mass_aq = conc_aq * (theta_layer * depth * cellarea())
            # if layer == 2:
            #     model.report(conc_aq, "Caq2")
            #     model.report(mass_aq, "Maq2")
            #     model.report(mass, "Mtot2")
            #     model.report(cellarea(), "cAr2")
            #     model.report(depth, "depth2")

        else:  # No gas phase
            # Whelan, 1987
            conc_aq = max(scalar(0), mass / (cellarea() * depth * theta_layer * retard_layer))

    return conc_aq


def getConcAds(model, layer, mass, gas=True):
    # mass / Kg soil
    depth = model.layer_depth[layer]
    if gas:
        theta_gas = max(model.theta_sat[layer] - model.theta[layer], scalar(0))
        # [mass pest/Kg soil]
        conc_ads = max(scalar(0), mass / ((cellarea() * depth) *
                           (theta_gas / (model.k_h * model.k_d) +
                            model.theta[layer] / model.k_d +
                            model.p_b)))
    else:
        print("No implementation without gas available")
        raise NotImplementedError

    return conc_ads


def getVolatileMass(model, temp_air, mass,  # frac,
                    rel_diff_model="option-2", sorption_model="linear",
                    gas=True, run=True):
    if not run:
        volat_flux = model.zero_map
    else:  # Volatilize only during peak volatilization time i.e., first 24 hrs, @Prueger2005.
        layer = 0
        theta_gas = max(model.theta_sat[layer] - model.theta[layer], scalar(0))
        # Convert to m (needed for final mass computation on cell basis)
        depth_m = model.layer_depth[0] * 1 / 10 ** 3
        # Air boundary layer, assumed as 2m high
        thickness_a = scalar(1.0)  # m
        # Diffusion coefficient in air (cm^2/s); https://www.gsi-net.com
        #  D_ar (metolachlor) = 0.03609052694,  at reference Temp., in Kelvin, D_a,r)
        diff_ar = 0.03609052694 * 86400.0 * 1.0 / 10 ** 4  # m2/d
        # Diffusion coefficient adjusted to air Temp. in Kelvin, D_a
        diff_a = ((temp_air + 273.15) / 293.15) ** 1.75 * diff_ar  # m2/d

        if rel_diff_model == "option-2":
            # Millington and Quirk, 1960 (in Leistra, 2001, p.48)
            # a,b parameters: Jin and Jury, 1996 (in Leistra, 2001)
            diff_relative_gas = max((diff_a * theta_gas ** 2 /
                                 model.theta_sat[layer] ** (2 / 3)), scalar(1e10-6))  # m2/d
        elif rel_diff_model == "option-1":
            # Currie 1960 (in Leistra, 2001)
            # a,b parameters: Baker, 1987 (in Leistra, 2001)
            diff_relative_gas = max((diff_a * 2.5 * theta_gas ** 3), scalar(1e10-6))  # m2/d
        else:
            print("No appropriate relative diffusion parameter chosen")
            diff_relative_gas = diff_a  # m2/d
        # Transport resistance through air (r_a) and soil (r_s) layer
        r_a = thickness_a / diff_a  # d/m
        r_s = max(scalar(0), (0.5 * depth_m) / diff_relative_gas)  # d/m

        conc_aq = getConcAq(model, layer, mass, sorption_model=sorption_model, gas=gas)
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
    #   2 Pa = 2 N/(m s^1) = 2 Kg/(m s^1)
    #   Convert to g/(cm s): dyn_visc = 8.9e-03 [g/cm s]
    dyn_visc = 8.9e-03  # [g/cm s] @25 degrees, [@Shi2011]:\mu

    # Solute diffusivity in water (D_w)
    # Metolachlor = 5.0967719112e-006 (cm2 / s)
    diff_solute = 5.0967719112e-006  # [cm2 / s], [@Shi2011]:D_w
    Sc = dyn_visc / (model.p_b * diff_solute)  # (-) Schmidt number, [@Shi2011]:S_c

    # Reynolds number (dimensionless), 86400s = 2 day
    cell_length = 2 * 10 ** 3  # mm
    # Reynolds (Re), [-] (Shi et al., 2011)
    re = (model.p_b * 1 / 10 ** 2 * runoffvelocity * cell_length) / (dyn_visc * 86400)
    kl = (0.664 * ((diff_solute * 86400 * 10 ** 2) / cell_length) *
          re ** (float(1) / 2) * Sc ** (float(1) / 3))
    return kl  # mm/day


def getRunOffMass(model, precip, runoff_mm, mass,
                  transfer_model="simple-mt", sorption_model="linear",
                  gas=True, debug=False, run=True):
    if not run:
        mass_ro = model.zero_map
    elif debug:
        mass_ro = model.zero_map
    else:
        # Aqueous concentration
        layer = 0
        conc_aq = getConcAq(model, layer, mass, sorption_model=sorption_model, gas=gas)

        if transfer_model == "simple-mt":
            mass_ro = conc_aq * runoff_mm * cellarea()
        elif transfer_model == "nu-mlm-ro":
            # non-uniform-mixing-layer-model-runoff (nu-mlm-ro)
            # Considers a decrease in effective transfer as mixing layer depth increases
            # Adapted from Ahuja and Lehman, 1983 in @Shi2011,
            # Adaptation replaces Precip by Runoff amount.
            b = 1  # [mm] Calibration constant, 2 >= b > 0 (b-ranges appear reasonable).
            # As b decreases, mass transfer increases, model.z0 in mm
            mass_ro = conc_aq * (runoff_mm * cellarea()) * exp(-b * model.layer_depth[layer])
        elif transfer_model == "nu-mlm":
            # non-uniform-mixing-layer-model (nu-mlm)
            # Original from Ahuja and Lehman, 1983 in @Shi2011
            b = 1  # [mm] Calibration constant, 2 >= b > 0 (b-ranges appear reasonable).
            # As b decreases, mass transfer increases, model.z0 in mm
            mass_ro = conc_aq * (precip * cellarea()) * exp(-b * model.layer_depth[layer])
            mass_ro = ifthenelse(runoff_mm > scalar(0), mass_ro, scalar(0))
        elif transfer_model == "d-mlm":
            # distributed mixing-layer-model (d-mlm)
            # Adapted from Havis et al., 1992, and
            # taking the K_L definition for laminar flow from Bennett and Myers, 1982.
            mass_ro = getKfilm(model, runoff_mm) * cellarea() * conc_aq
        else:
            print("Run-off transfer model not stated")
            return None

        mReport = False
        if mReport:
            pass
            # model.report(conc_aq, 'aCo')
            # model.report(mass_ro, 'aMROa')
            # model.report(runoff_mm, 'aRO')

    return mass_ro


def getLeachedMass(model, layer, water_flux,
                   mass,
                   sorption_model=None,
                   leach_model=None, gas=True, debug=False, run=True):
    if not run:
        mass_leached = deepcopy(model.zero_map)
    elif debug:
        mass_leached = deepcopy(model.zero_map)
    else:
        if mapminimum(model.theta[layer]) < scalar(1e-06):
            mass_leached = deepcopy(model.zero_map)
        else:
            theta_layer = model.theta[layer]
            depth = model.layer_depth[layer]

            # Aqueous concentration
            conc_aq = getConcAq(model, layer, mass,
                                sorption_model=sorption_model, gas=gas)

            # Mass available for transport
            mass_aq = conc_aq * (theta_layer * depth * cellarea())

            # if layer == 2:
            #     model.report(conc_aq, "Caq")
            #     model.report(mass_aq, "Maq")
            #     model.report(mass, "Mtot")
            #     model.report(cellarea(), "cAr")
            #     model.report(depth, "depth1")

            test = mass - mass_aq
            if mapminimum(test) < 0:
                print("Error, mass < mass_aq, on layer: ", str(layer))
                model.report(test, 'aMzErr' + str(layer))

            if sorption_model == "linear":
                # Retardation factor
                retard_layer = scalar(1) + (model.p_b * model.k_d) / theta_layer
            else:
                print("No sorption assumed, Ret. factor = 2")
                retard_layer = scalar(1)  # No retardation.

            if leach_model == "mcgrath":
                if layer == 0:
                    mass_aq_new = mass_aq * exp(-water_flux / (theta_layer * retard_layer * depth))
                    mass_leached = mass_aq - mass_aq_new
                    if mapminimum(mass_leached) < -1e-06:
                        print("Error in Leached Model, layer: ", str(layer))
                        model.report(mass_leached, 'aZ' + str(layer) + 'LCH')
                else:
                    # mass_aq_new = mass_aq * exp(-water_flux / (theta_layer * retard_layer * depth))
                    # mass_leached = mass_aq - mass_aq_new

                    # McGrath not used in lower layers,
                    # as formulation accounts for rainfall impact
                    max_flux = max(min(water_flux, (theta_layer - model.theta_fc[layer]) * depth), scalar(0))
                    mass_leached = conc_aq * max_flux * cellarea()

                    mass_aq = conc_aq * (theta_layer * depth * cellarea())
                    mass_aq_new = mass_aq - mass_leached
                    if mapminimum(mass_aq_new) < -1e-06:
                        print("Error in Leached Model, layer: ", str(layer))
                        model.report(mass_leached, 'aZ' + str(layer) + 'LCH')
                    if mapminimum(mass_aq_new) < 0:
                        print("Err mass_aq_new")
                        mass_leached = max(mass_leached, scalar(0))
                        # mass_aq_new = mass_aq - mass_leached
            else:
                mass_leached = conc_aq * water_flux * cellarea()
                mass_aq = conc_aq * (theta_layer * depth) * cellarea()
                mass_aq_new = mass_aq - mass_leached
                if mapminimum(mass_aq_new) < 0:
                    print("Error in Leached Model")

    return mass_leached


def getLatMassFlux(model, layer, mass, flux_map_mm,
                   sorption_model='linear', gas=True,
                   debug=False, run=True):
    """
    :param model:
    :param layer:
    :param mass:
    :param sorption_model:
    :param gas:
    :param debug:
    :param run:
    :return:
    """
    if not run:
        latflux_dict = {
            'mass_loss': deepcopy(model.zero_map),
            'mass_gain': deepcopy(model.zero_map),
            'new_mass': deepcopy(model.zero_map)
        }
    else:
        if mapminimum(model.theta[layer]) < scalar(1e-06) or layer == (model.num_layers - 1):
            latflux_dict = {
                'mass_loss': deepcopy(model.zero_map),
                'mass_gain': deepcopy(model.zero_map),
                'new_mass': mass
            }
        else:

            # Aqueous concentration  mass/L
            conc_aq = getConcAq(model, layer, mass, sorption_model=sorption_model, gas=gas)

            mass_loss = conc_aq * flux_map_mm * cellarea()  # mm * m2 = L
            mass_gain = upstream(model.ldd_subs, mass_loss)
            new_mass = mass - mass_loss + mass_gain
            # http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_upstream.html

            if debug:
                model.report(mass, 'aMi' + str(layer))
                model.report(mass_loss, 'aMloss' + str(layer))
                model.report(mass_gain, 'aMgain' + str(layer))
                # aguila --scenarios='{2}' --timesteps=[2,300,2] aMi0 aMloss0 aMgain0

            latflux_dict = {
                'mass_loss': mass_loss,
                'mass_gain': mass_gain,
                'new_mass': new_mass
            }
    return latflux_dict


def getLatMassFluxManfreda(model, layer, mass, cell_moisture_outflow, upstream_cell_inflow,
                   sorption_model='linear', gas=True,
                   debug=False, run=True):
    """
    :param model:
    :param layer:
    :param mass:
    :param cell_moisture_outflow: mm
    :param upstream_cell_inflow: mm
    :param sorption_model:
    :param gas:
    :param debug:
    :param run:
    :return:
    """
    if not run:
        latflux_dict = {
            'mass_loss': deepcopy(model.zero_map),
            'mass_gain': deepcopy(model.zero_map),
            'net_mass_latflux': deepcopy(model.zero_map)
        }
    else:
        theta_layer = model.theta[layer]
        if mapminimum(model.theta[layer]) < scalar(1e-06) or layer == (model.num_layers - 1):
            latflux_dict = {
                'mass_loss': deepcopy(model.zero_map),
                'mass_gain': deepcopy(model.zero_map),
                'net_mass_latflux': deepcopy(model.zero_map)
            }
        else:
            depth = model.layer_depth[layer]
            c = model.c_lf[layer]

            # Aqueous concentration  ug/L
            conc_aq = getConcAq(model, layer, mass, sorption_model=sorption_model, gas=gas)

            # W(j/i)
            rel_wetness = model.wetness / accuflux(model.ldd_subs, model.wetness)
            # # Cell mass loss/gain (to update only mass)
            mass_loss = max(conc_aq * (c * (depth * theta_layer - depth * model.theta_fc[layer])), scalar(0))
            mass_gain = rel_wetness * accuflux(model.ldd_subs, mass_loss)

            # mass_gain = max(upstream(model.ldd_subs, conc_aq * upstream_cell_inflow * cellarea()), scalar(0))
            # mass_loss = max(downstream(model.ldd_subs, mass_gain), scalar(0))
            # http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_upstream.html
            net_mass_latflux = mass_gain - mass_loss

            if debug:
                model.report(mass, 'aMi' + str(layer))
                model.report(mass_loss, 'aMloss' + str(layer))
                model.report(mass_gain, 'aMgain' + str(layer))
                # aguila --scenarios='{2}' --timesteps=[2,300,2] aMi0 aMloss0 aMgain0


            latflux_dict = {
                'mass_loss': mass_loss,
                'mass_gain': mass_gain,
                'net_mass_latflux': net_mass_latflux
            }
    return latflux_dict


def getDrainMassFlux(model, layer, mass,
                     sorption_model='linear', gas=True,
                     debug=False, run=True):
    if not run:
        mass_loss = model.zero_map
    elif debug:
        mass_loss = model.zero_map
    else:
        # Aqueous concentration
        conc_aq = getConcAq(model, layer, mass, sorption_model=sorption_model, gas=gas)
        mass_loss = max(conc_aq * (model.c_adr * (model.layer_depth[layer] * model.theta[layer] -
                                                  model.layer_depth[layer] * model.theta_fc[layer])),
                        scalar(0))
    return mass_loss


def getMassDegradation(model, layer, theta_wp, mass,
                       frac="L",
                       sor_deg_factor=1,
                       sorption_model="linear",
                       gas=True, debug=False, run=True):
    if not run:
        return {"mass_tot_new": mass,
                "mass_deg_aq": model.zero_map,
                "mass_deg_ads": model.zero_map}
    else:
        theta_layer = model.theta[layer]
        depth = model.layer_depth[layer]

        # Mass compartment (bio-available fraction)
        # Was considering a step function, however,
        # a continuous one is implemented below instead
        aged_frac = ifthenelse(model.aged_days > 100, scalar(1),
                               ifthenelse(model.aged_days > 75, scalar(0.75),
                                          ifthenelse(model.aged_days > 50, scalar(0.50),
                                                     ifthenelse(model.aged_days > 25, scalar(0.25),
                                                                scalar(0)))))
        # bioa_mass = mass * aged_frac
        bioa_mass = mass * exp(-model.aged_days / model.dt_50_ref)
        aged_mass = mass - bioa_mass
        # bioa_mass = deepcopy(mass)
        # aged_mass = deepcopy(mass)

        # F_Theta_1
        theta_factor = ifthenelse(theta_layer <= 0.5 * theta_wp, scalar(0),
                                  ifthenelse(theta_layer <= model.theta_fc[layer],
                                             (((theta_layer - 0.5 * theta_wp) / (
                                                 model.theta_fc[layer] - theta_wp)) ** scalar(model.beta_moisture)),
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
        conc_aq = getConcAq(model, layer, bioa_mass,
                            sorption_model=sorption_model, gas=gas)  # mass/L
        conc_ads = getConcAds(model, layer, bioa_mass, gas=gas)  # mass/g soil
        conc_gas = conc_aq / model.k_h

        mass_aq = conc_aq * (theta_layer * depth * cellarea())
        mass_ads = conc_ads * (depth * cellarea() * model.p_b)  # pb = g/cm3

        # Step 2 - Degrade phase fractions
        # First order degradation kinetics
        theta_gas = max(model.theta_sat[layer] - theta_layer, scalar(0))
        if frac == "H":
            conc_aq_new = conc_aq * exp(-1 * model.alpha_iso * k_b * scalar(model.jd_dt))
            conc_ads_new = conc_ads * exp(-1 * model.alpha_iso * k_bs * scalar(model.jd_dt))
        else:  # Same for total conc as for "L"
            conc_aq_new = conc_aq * exp(-1 * k_b * scalar(model.jd_dt))
            conc_ads_new = conc_ads * exp(-1 * k_bs * scalar(model.jd_dt))

        # Step 2 - Convert to mass (i.e., after degradation in each phase)
        mass_aq_new = conc_aq_new * (theta_layer * depth * cellarea())
        mass_ads_new = conc_ads_new * (model.p_b * depth * cellarea())  # pb = g/cm3
        mass_gas = conc_gas * (theta_gas * depth * cellarea())
        mass_tot_new = mass_aq_new + mass_ads_new + mass_gas + aged_mass
        mass_deg_aq = mass_aq - mass_aq_new
        mass_deg_ads = mass_ads - mass_ads_new
        # if frac == "L" and layer == 0:
        #     model.report(mass, 'adt0M')
        #     model.report(mass_tot_new, 'adt1M')
        #   # model.report(mass_aq_new, 'adt1M')
        if debug:
            model.report(bioa_mass, 'aBIOz' + str(layer))
            model.report(theta_factor, 'aFz' + str(layer))
            model.report(temp_factor, 'aTz' + str(layer))
            model.report(dt_50, 'aDT50z' + str(layer))
            model.report(k_b, 'akbz' + str(layer))
        # aguila --scenarios='{2}' --timesteps=[2,300,2] aFz aTz aDT50z akbz

    return {"mass_tot_new": mass_tot_new,
            "mass_deg_aq": mass_deg_aq,
            "mass_deg_ads": mass_deg_ads}
