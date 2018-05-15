# -*- coding: utf-8 -*-
from pcraster.framework import *


def reportNashConcComposites(model, north_ave_conc, valley_ave_conc, south_ave_conc):
    # Nash
    conc_north_obs = timeinputscalar('northConc.tss', ordinal("north_ave"))
    model.northConc_diff += ifthenelse(conc_north_obs > 0, (north_ave_conc - conc_north_obs) ** 2, scalar(0))
    model.northConc_var += ifthenelse(conc_north_obs > 0, (north_ave_conc - scalar(1.909193)) ** 2,
                                     scalar(0))  # ug/g
    
    # Nash
    conc_valley_obs = timeinputscalar('valleyConc.tss', ordinal("valley_ave"))
    model.valleyConc_diff += ifthenelse(conc_valley_obs > 0, (valley_ave_conc - conc_valley_obs) ** 2, scalar(0))
    model.valleyConc_var += ifthenelse(conc_valley_obs > 0, (valley_ave_conc - scalar(2.261839)) ** 2,
                                      scalar(0))  # ug/g

    conc_south_obs = timeinputscalar('southConc.tss', ordinal("south_ave"))
    model.southConc_diff += ifthenelse(conc_south_obs > 0, (south_ave_conc - conc_south_obs) ** 2, scalar(0))
    model.southConc_var += ifthenelse(conc_south_obs > 0, (south_ave_conc - scalar(2.389668)) ** 2,
                                     scalar(0))  # ug/g
    nash_compConc_L = 1 - ((model.northConc_diff / model.northConc_var) * 1 / 3 +
                           (model.valleyConc_diff / model.valleyConc_var) * 1 / 3 +
                           (model.southConc_diff / model.southConc_var) * 1 / 3)
    model.nash_compConc_L_tss.sample(nash_compConc_L)


def reportNashDeltaComposites():
    pass
