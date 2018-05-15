# -*- coding: utf-8 -*-
from pcraster.framework import *


"""
Nash Soil Concentrations
1) Get the mean for each soil composite for entire year
North = 1.909193 ug/g soil
Talweg = 2.261839 ug/g soil
South = 2.389668 ug/g soil

2) The variance for each transect is
var_north = (conc_north - mean_north)**2, if conc_north > 0

3) Nash will be:
1 - (conc_north_diff/var_north +  valley + south)
"""

# TODO: define codes for all transects...
north_plot_codes = ['n1', 'n2', 'n3', 'n4', 'n5', 'n7', 'n8']  # no 'n6'!
north_plots = len(north_plot_codes)


def reportNorthSoils(model, cell_mass, cell_massXdelta, sampling_pts):
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
    north_ave_conc = 10e6 * (north_ave_mass /
                             (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil

    n1_ave_conc = 10e6 * (n1_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n2_ave_conc = 10e6 * (n2_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n3_ave_conc = 10e6 * (n3_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n4_ave_conc = 10e6 * (n4_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n5_ave_conc = 10e6 * (n5_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n7_ave_conc = 10e6 * (n7_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    n8_ave_conc = 10e6 * (n8_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil

    # Record
    model.north_conc_tss.sample(north_ave_conc)  # 30 obs
    model.n1_conc_tss.sample(n1_ave_conc)
    model.n2_conc_tss.sample(n2_ave_conc)
    model.n3_conc_tss.sample(n3_ave_conc)
    model.n4_conc_tss.sample(n4_ave_conc)
    model.n5_conc_tss.sample(n5_ave_conc)
    model.n7_conc_tss.sample(n7_ave_conc)
    model.n8_conc_tss.sample(n8_ave_conc)

    model.north_d13C_tss.sample(north_d13C)
    model.n1_d13C_tss.sample(n1_d13C)
    model.n2_d13C_tss.sample(n2_d13C)
    model.n3_d13C_tss.sample(n3_d13C)
    model.n4_d13C_tss.sample(n4_d13C)
    model.n5_d13C_tss.sample(n5_d13C)
    model.n7_d13C_tss.sample(n7_d13C)
    model.n8_d13C_tss.sample(n8_d13C)

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

    valley_ave_conc = 10e6 * (valley_ave_mass /
                              (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t4_ave_conc = 10e6 * (t4_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t5_ave_conc = 10e6 * (t5_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t7_ave_conc = 10e6 * (t7_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t8_ave_conc = 10e6 * (t8_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t9_ave_conc = 10e6 * (t9_ave_mass /
                          (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    t10_ave_conc = 10e6 * (t10_ave_mass /
                           (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil

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

    south_ave_conc = 10e6 * (south_ave_mass /
                             (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    s11_ave_conc = 10e6 * (s11_ave_mass /
                           (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    s12_ave_conc = 10e6 * (s12_ave_mass /
                           (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    s13_ave_conc = 10e6 * (s13_ave_mass /
                           (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
    
    test = False
    if test:
        # Test for variation in sample points
        cell_conc = 10e6 * (cell_mass /
                            (cellarea() * model.smp_depth)) * 1 / (model.p_b * 10e03)  # ug/g soil
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
