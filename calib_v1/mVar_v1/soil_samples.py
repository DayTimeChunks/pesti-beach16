# -*- coding: utf-8 -*-
from pcraster.framework import *

"""
Nash Soil Concentrations
2) Get the mean for each soil composite for entire year
North = 2.909193 ug/g soil
Talweg = 1.261839 ug/g soil
South = 1.389668 ug/g soil

1) The variance for each transect is
var_north = (conc_north - mean_north)**1, if conc_north > 0

3) Nash will be:
2 - (conc_north_diff/var_north +  valley + south)
"""

# TODO: define codes for all transects...
north_plot_codes = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8']  # no 'n6'!
north_plots = len(north_plot_codes)


def importPlotMaps(model):
    # Points model to sampling points (pixels) on a given plot
    # Defines a reference to each map in dictionary: 'plot_maps'.
    try:
        model.plot_maps
    except AttributeError:
        model.plot_maps = dict()

    plots = ['north', 'n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8',
             'valley', 'v4', 'v5', 'v7', 'v8', 'v9', 'v10',
             'south', 's11', 's12', 's13']

    for plot in range(len(plots)):
        plot_map = 'mapAnalysis/' + plots[plot] + '_nom'
        model.plot_maps[plot_map] = nominal(model.readmap(plot_map))


def defineSoilTSS(model):
    try:
        model.soil_dict
    except AttributeError:
        model.soil_dict = dict()

    transects = ['north', 'valley', 'south']
    for tr in range(len(transects)):
        transect = transects[tr]
        transect_conc = transect[0:3] + 'CONC'
        transect_conc_real = transect[0:3] + 'CONC_real'
        transect_conc_aged = transect[0:3] + 'CONC_aged'
        transect_delta = transect[0:3] + 'd13C'
        transect_delta_real = transect[0:3] + 'd13C_real'
        transect_delta_aged = transect[0:3] + 'd13C_aged'
        transect_map = 'mapAnalysis/' + transect + '_ave'

        # Concentrations
        model.soil_dict[transect_conc] = TimeoutputTimeseries("resM_" + transect_conc, model, ordinal(transect_map),
                                                              noHeader=False)
        model.soil_dict[transect_conc_real] = TimeoutputTimeseries("resM_" + transect_conc_real, model,
                                                                   ordinal(transect_map), noHeader=False)
        model.soil_dict[transect_conc_aged] = TimeoutputTimeseries("resM_" + transect_conc_aged, model,
                                                                   ordinal(transect_map), noHeader=False)
        # Isotopes
        model.soil_dict[transect_delta] = TimeoutputTimeseries("resM_" + transect_delta, model, ordinal(transect_map),
                                                               noHeader=False)
        model.soil_dict[transect_delta_real] = TimeoutputTimeseries("resM_" + transect_delta_real, model,
                                                                    ordinal(transect_map), noHeader=False)
        model.soil_dict[transect_delta_aged] = TimeoutputTimeseries("resM_" + transect_delta_aged, model,
                                                                    ordinal(transect_map), noHeader=False)
    plots = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8',
             'v4', 'v5', 'v7', 'v8', 'v9', 'v10',
             's11', 's12', 's13']
    for plot in range(len(plots)):
        plot_name = plots[plot]
        plot_map = 'mapAnalysis/' + plot_name + '_out'  # 1-pixel map

        plot_conc = plot_name + 'CONC'
        plot_conc_real = plot_name + 'CONC_real'
        plot_conc_aged = plot_name + 'CONC_aged'

        plot_delta = plot_name + 'd13C'
        plot_delta_real = plot_name + 'd13C_real'
        plot_delta_aged = plot_name + 'd13C_aged'

        # Concentrations
        model.soil_dict[plot_conc] = TimeoutputTimeseries("resM_" + plot_conc, model, ordinal(plot_map), noHeader=False)
        model.soil_dict[plot_conc_real] = TimeoutputTimeseries("resM_" + plot_conc_real, model, ordinal(plot_map),
                                                               noHeader=False)
        model.soil_dict[plot_conc_aged] = TimeoutputTimeseries("resM_" + plot_conc_aged, model, ordinal(plot_map),
                                                               noHeader=False)

        # Delta
        model.soil_dict[plot_delta] = TimeoutputTimeseries("resM_" + plot_delta, model, ordinal(plot_map),
                                                           noHeader=False)
        model.soil_dict[plot_delta_real] = TimeoutputTimeseries("resM_" + plot_delta_real, model, ordinal(plot_map),
                                                                noHeader=False)
        model.soil_dict[plot_delta_aged] = TimeoutputTimeseries("resM_" + plot_delta_aged, model, ordinal(plot_map),
                                                                noHeader=False)
        # Example:
        # model.soil_dict[n1_d13C] = TimeoutputTimeseries("resM_n1d13C", model, ordinal("n1_out"), noHeader=False)


def reportSoilTSS(model, cell_mass, cell_massXdelta, transect, type='real'):
    if transect == 'north':
        plots = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8']
        plot_sampling_pts = [4, 6, 6, 4, 4, 5, 5]
        transect_sampling_pts = 30
    elif transect == 'valley':
        plots = ['v4', 'v5', 'v7', 'v8', 'v9', 'v10']
        plot_sampling_pts = [4, 4, 5, 5, 5, 5]
        transect_sampling_pts = 25
    else:
        transect = 'south'
        plots = ['s11', 's12', 's13']
        plot_sampling_pts = [8, 7, 5]
        transect_sampling_pts = 26

    # Record Transect
    transect_map = 'mapAnalysis/' + transect + '_nom'
    transect_tot_mass = areatotal(cell_mass, model.plot_maps[transect_map])
    transect_ave_mass = transect_tot_mass / scalar(transect_sampling_pts)
    transect_ave_conc = 1e6 * (transect_ave_mass /
                               (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    transect_d13C = areatotal(cell_massXdelta, model.plot_maps[transect_map]) / transect_tot_mass

    if type == 'bioavail':
        transect_conc = transect[0:3] + 'CONC'
        transect_delta = transect[0:3] + 'd13C'
    elif type == 'real':
        transect_conc = transect[0:3] + 'CONC_real'
        transect_delta = transect[0:3] + 'd13C_real'
    else:
        transect_conc = transect[0:3] + 'CONC_aged'
        transect_delta = transect[0:3] + 'd13C_aged'

    model.soil_dict[transect_conc].sample(transect_ave_conc)
    model.soil_dict[transect_delta].sample(transect_d13C)

    # Record detailed
    assert len(plots) == len(plot_sampling_pts)
    for plot in range(len(plots)):
        plot_name = plots[plot]
        plot_map = 'mapAnalysis/' + plot_name + '_nom'
        tot_mass = areatotal(cell_mass, model.plot_maps[plot_map])
        plot_d13C = areatotal(cell_massXdelta, model.plot_maps[plot_map]) / tot_mass
        plot_ave_mass = tot_mass / scalar(plot_sampling_pts[plot])
        plot_ave_conc = 1e6 * (plot_ave_mass /
                               (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil

        # Define dictionary key
        if type == 'bioavail':
            plot_conc = plot_name + 'CONC'
            plot_delta = plot_name + 'd13C'
        elif type == 'real':
            plot_conc = plot_name + 'CONC_real'
            plot_delta = plot_name + 'd13C_real'
        else:
            plot_conc = plot_name + 'CONC_aged'
            plot_delta = plot_name + 'd13C_aged'

        # Sample delta and concentrations
        model.soil_dict[plot_conc].sample(plot_ave_conc)
        model.soil_dict[plot_delta].sample(plot_d13C)

    # n1_tot_mass = areatotal(cell_mass, model.n1_plot)
    # n2_tot_mass = areatotal(cell_mass, model.n2_plot)
    # n3_tot_mass = areatotal(cell_mass, model.n3_plot)
    # n4_tot_mass = areatotal(cell_mass, model.n4_plot)
    # n5_tot_mass = areatotal(cell_mass, model.n5_plot)
    # n7_tot_mass = areatotal(cell_mass, model.n7_plot)
    # n8_tot_mass = areatotal(cell_mass, model.n8_plot)


    # n1_d13C = areatotal(cell_massXdelta, model.n1_plot) / n1_tot_mass
    # n2_d13C = areatotal(cell_massXdelta, model.n2_plot) / n2_tot_mass
    # n3_d13C = areatotal(cell_massXdelta, model.n3_plot) / n3_tot_mass
    # n4_d13C = areatotal(cell_massXdelta, model.n4_plot) / n4_tot_mass
    # n5_d13C = areatotal(cell_massXdelta, model.n5_plot) / n5_tot_mass
    # n7_d13C = areatotal(cell_massXdelta, model.n7_plot) / n7_tot_mass
    # n8_d13C = areatotal(cell_massXdelta, model.n8_plot) / n8_tot_mass


    # n1_ave_mass = n1_tot_mass / scalar(4)
    # n2_ave_mass = n2_tot_mass / scalar(6)
    # n3_ave_mass = n3_tot_mass / scalar(6)
    # n4_ave_mass = n4_tot_mass / scalar(4)
    # n5_ave_mass = n5_tot_mass / scalar(4)
    # n7_ave_mass = n7_tot_mass / scalar(5)
    # n8_ave_mass = n8_tot_mass / scalar(5)

    # Concentrations
    # Mass is converted to ug to compare against sampled data.


    # n1_ave_conc = 1e6 * (n1_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n2_ave_conc = 1e6 * (n2_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n3_ave_conc = 1e6 * (n3_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n4_ave_conc = 1e6 * (n4_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n5_ave_conc = 1e6 * (n5_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n7_ave_conc = 1e6 * (n7_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    # n8_ave_conc = 1e6 * (n8_ave_mass /
    #                      (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil

    # Record
    # if real:
    #     model.north_conc_real_tss.sample(north_ave_conc)  # 30 obs
    # else:


    # model.north_conc_tss.sample(north_ave_conc)  # 30 obs
    # model.n1_conc_tss.sample(n1_ave_conc)
    # model.n2_conc_tss.sample(n2_ave_conc)
    # model.n3_conc_tss.sample(n3_ave_conc)
    # model.n4_conc_tss.sample(n4_ave_conc)
    # model.n5_conc_tss.sample(n5_ave_conc)
    # model.n7_conc_tss.sample(n7_ave_conc)
    # model.n8_conc_tss.sample(n8_ave_conc)
    #
    # model.north_d13C_tss.sample(north_d13C)
    # model.n1_d13C_tss.sample(n1_d13C)
    # model.n2_d13C_tss.sample(n2_d13C)
    # model.n3_d13C_tss.sample(n3_d13C)
    # model.n4_d13C_tss.sample(n4_d13C)
    # model.n5_d13C_tss.sample(n5_d13C)
    # model.n7_d13C_tss.sample(n7_d13C)
    # model.n8_d13C_tss.sample(n8_d13C)

    # Report
    return {'ave_conc': transect_ave_conc,
            'd13C': transect_d13C}


def reportNorthSoils2(model, cell_mass, cell_massXdelta, sampling_pts, real=True):
    north_tot_mass = areatotal(cell_mass, model.north_wk)
    n1_tot_mass = areatotal(cell_mass, model.n1_plot)
    n2_tot_mass = areatotal(cell_mass, model.n2_plot)
    n3_tot_mass = areatotal(cell_mass, model.n3_plot)
    n4_tot_mass = areatotal(cell_mass, model.n4_plot)
    n5_tot_mass = areatotal(cell_mass, model.n5_plot)
    n7_tot_mass = areatotal(cell_mass, model.n7_plot)
    n8_tot_mass = areatotal(cell_mass, model.n8_plot)

    north_d13C = areatotal(cell_massXdelta, model.north_wk) / north_tot_mass
    n1_d13C = areatotal(cell_massXdelta, model.n1_plot) / n1_tot_mass
    n2_d13C = areatotal(cell_massXdelta, model.n2_plot) / n2_tot_mass
    n3_d13C = areatotal(cell_massXdelta, model.n3_plot) / n3_tot_mass
    n4_d13C = areatotal(cell_massXdelta, model.n4_plot) / n4_tot_mass
    n5_d13C = areatotal(cell_massXdelta, model.n5_plot) / n5_tot_mass
    n7_d13C = areatotal(cell_massXdelta, model.n7_plot) / n7_tot_mass
    n8_d13C = areatotal(cell_massXdelta, model.n8_plot) / n8_tot_mass

    north_ave_mass = north_tot_mass / scalar(sampling_pts)
    n1_ave_mass = n1_tot_mass / scalar(4)
    n2_ave_mass = n2_tot_mass / scalar(6)
    n3_ave_mass = n3_tot_mass / scalar(6)
    n4_ave_mass = n4_tot_mass / scalar(4)
    n5_ave_mass = n5_tot_mass / scalar(4)
    n7_ave_mass = n7_tot_mass / scalar(5)
    n8_ave_mass = n8_tot_mass / scalar(5)

    # Concentrations
    # Mass is converted to ug to compare against sampled data.
    north_ave_conc = 1e6 * (north_ave_mass /
                            (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil

    n1_ave_conc = 1e6 * (n1_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n2_ave_conc = 1e6 * (n2_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n3_ave_conc = 1e6 * (n3_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n4_ave_conc = 1e6 * (n4_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n5_ave_conc = 1e6 * (n5_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n7_ave_conc = 1e6 * (n7_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil
    n8_ave_conc = 1e6 * (n8_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e03)  # ug/g soil

    # Record
    # if real:
    #     model.north_conc_real_tss.sample(north_ave_conc)  # 30 obs
    # else:
    transect = 'north'
    transect_conc = transect[0:3] + 'CONC'
    transect_conc_real = transect[0:3] + 'CONC_real'
    transect_delta = transect[0:3] + 'd13C'
    transect_delta_real = transect[0:3] + 'd13C_real'
    model.soil_dict[transect_conc].sample(north_ave_conc)
    model.soil_dict[transect_delta].sample(north_d13C)

    plots = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8']
    for plot in range(len(plots)):
        plot_map = 'mapAnalysis/' + plots[plot] + '_nom'
        plot_name = plots[plot]
        plot_conc = plot_name + 'CONC'
        plot_conc_real = plot_name + 'CONC_real'
        plot_delta = plot_name + 'd13C'
        plot_delta_real = plot_name + 'd13C_real'
        model.soil_dict[plot_conc].sample()

    # model.north_conc_tss.sample(north_ave_conc)  # 30 obs
    # model.n1_conc_tss.sample(n1_ave_conc)
    # model.n2_conc_tss.sample(n2_ave_conc)
    # model.n3_conc_tss.sample(n3_ave_conc)
    # model.n4_conc_tss.sample(n4_ave_conc)
    # model.n5_conc_tss.sample(n5_ave_conc)
    # model.n7_conc_tss.sample(n7_ave_conc)
    # model.n8_conc_tss.sample(n8_ave_conc)
    #
    # model.north_d13C_tss.sample(north_d13C)
    # model.n1_d13C_tss.sample(n1_d13C)
    # model.n2_d13C_tss.sample(n2_d13C)
    # model.n3_d13C_tss.sample(n3_d13C)
    # model.n4_d13C_tss.sample(n4_d13C)
    # model.n5_d13C_tss.sample(n5_d13C)
    # model.n7_d13C_tss.sample(n7_d13C)
    # model.n8_d13C_tss.sample(n8_d13C)

    # Report
    return {'ave_conc': north_ave_conc,
            'd13C': north_d13C}


def reportValleySoils(model, cell_mass, cell_massXdelta, sampling_pts):
    valley_tot_mass = areatotal(cell_mass, model.valley_wk)
    t4_tot_mass = areatotal(cell_mass, model.t4_plot)
    t5_tot_mass = areatotal(cell_mass, model.t5_plot)
    t7_tot_mass = areatotal(cell_mass, model.t7_plot)
    t8_tot_mass = areatotal(cell_mass, model.t8_plot)
    t9_tot_mass = areatotal(cell_mass, model.t9_plot)
    t10_tot_mass = areatotal(cell_mass, model.t10_plot)

    valley_d13C = areatotal(cell_massXdelta, model.valley_wk) / valley_tot_mass
    t4_d13C = areatotal(cell_massXdelta, model.t4_plot) / t4_tot_mass
    t5_d13C = areatotal(cell_massXdelta, model.t5_plot) / t5_tot_mass
    t7_d13C = areatotal(cell_massXdelta, model.t7_plot) / t7_tot_mass
    t8_d13C = areatotal(cell_massXdelta, model.t8_plot) / t8_tot_mass
    t9_d13C = areatotal(cell_massXdelta, model.t9_plot) / t9_tot_mass
    t10_d13C = areatotal(cell_massXdelta, model.t10_plot) / t10_tot_mass

    valley_ave_mass = valley_tot_mass / scalar(sampling_pts)
    t4_ave_mass = t4_tot_mass / scalar(4)
    t5_ave_mass = t5_tot_mass / scalar(4)
    t7_ave_mass = t7_tot_mass / scalar(5)
    t8_ave_mass = t8_tot_mass / scalar(5)
    t9_ave_mass = t9_tot_mass / scalar(5)
    t10_ave_mass = t10_tot_mass / scalar(5)

    valley_ave_conc = 1e6 * (valley_ave_mass /
                             (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t4_ave_conc = 1e6 * (t4_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t5_ave_conc = 1e6 * (t5_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t7_ave_conc = 1e6 * (t7_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t8_ave_conc = 1e6 * (t8_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t9_ave_conc = 1e6 * (t9_ave_mass /
                         (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    t10_ave_conc = 1e6 * (t10_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil

    model.valley_conc_tss.sample(valley_ave_conc)  # 25 obs
    model.t4_conc_tss.sample(t4_ave_conc)  #
    model.t5_conc_tss.sample(t5_ave_conc)  #
    model.t7_conc_tss.sample(t7_ave_conc)  #
    model.t8_conc_tss.sample(t8_ave_conc)  #
    model.t9_conc_tss.sample(t9_ave_conc)  #
    model.t10_conc_tss.sample(t10_ave_conc)  #

    model.valley_dC13_tss.sample(valley_d13C)
    model.t4_d13C_tss.sample(t4_d13C)
    model.t5_d13C_tss.sample(t5_d13C)
    model.t7_d13C_tss.sample(t7_d13C)
    model.t8_d13C_tss.sample(t8_d13C)
    model.t9_d13C_tss.sample(t9_d13C)
    model.t10_d13C_tss.sample(t10_d13C)

    # Report
    return {'ave_conc': valley_ave_conc,
            'd13C': valley_d13C}


def reportSouthSoils(model, cell_mass, cell_massXdelta, sampling_pts):
    south_tot_mass = areatotal(cell_mass, model.south_wk)
    s11_tot_mass = areatotal(cell_mass, model.s11_plot)
    s12_tot_mass = areatotal(cell_mass, model.s12_plot)
    s13_tot_mass = areatotal(cell_mass, model.s13_plot)

    south_d13C = areatotal(cell_massXdelta, model.south_wk) / south_tot_mass
    s11_d13C = areatotal(cell_massXdelta, model.s11_plot) / s11_tot_mass
    s12_d13C = areatotal(cell_massXdelta, model.s12_plot) / s12_tot_mass
    s13_d13C = areatotal(cell_massXdelta, model.s13_plot) / s13_tot_mass

    south_ave_mass = south_tot_mass / scalar(sampling_pts)
    s11_ave_mass = s11_tot_mass / scalar(8)
    s12_ave_mass = s12_tot_mass / scalar(7)
    s13_ave_mass = s13_tot_mass / scalar(5)

    south_ave_conc = 1e6 * (south_ave_mass /
                            (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    s11_ave_conc = 1e6 * (s11_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    s12_ave_conc = 1e6 * (s12_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
    s13_ave_conc = 1e6 * (s13_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil

    test = False
    if test:
        # Test for variation in sample points
        cell_conc = 1e6 * (cell_mass /
                           (cellarea() * model.smp_depth)) * 1 / (model.p_b * 1e3)  # ug/g soil
        model.s11_smass_tss.sample(cell_mass)  # Should only increase due to upstream infux (if at all)
        model.s11_sconc_tss.sample(cell_conc)  # Should only increase after applciation

    model.south_conc_tss.sample(south_ave_conc)  # 26 obs
    model.s11_conc_tss.sample(s11_ave_conc)
    model.s12_conc_tss.sample(s12_ave_conc)
    model.s13_conc_tss.sample(s13_ave_conc)

    model.south_d13C_tss.sample(south_d13C)
    model.s11_d13C_tss.sample(s11_d13C)
    model.s12_d13C_tss.sample(s12_d13C)
    model.s13_d13C_tss.sample(s13_d13C)

    # Report
    return {'ave_conc': south_ave_conc,
            'd13C': south_d13C}
