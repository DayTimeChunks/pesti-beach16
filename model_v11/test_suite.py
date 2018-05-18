# -*- coding: utf-8 -*-
from pcraster.framework import *

# interval = 2
# def getPrintInterval(model, interval):
#     if model.currentTimeStep() % interval == 0:


# Layer depths
def checkLayerDepths(model, layer):
    model.report(model.layer_depth[layer], 'aDepth' + str(layer))
    # print('Mapminium, z' + str(layer) + ' ' + str(mapminimum(model.layer_depth[layer])))

#  aguila --scenarios='{1}' aDepth0 aDepth1 aDepth2 aDepth3

def checkRootDepths(model, root_depth_arr):
    for layer in range(len(root_depth_arr)):
        root_length = root_depth_arr[layer] / 10 ** 3  # Convert back to m
        model.report(root_length, 'aRDz'+str(layer))

def checkMoistureProps(model, prop_list, name):
    for i in range(len(prop_list)):
        mapname = str(name) + str(i)
        model.report(prop_list[i], mapname)

#  aguila --scenarios='{1}' --timesteps=[2,300,2] aSATz0 aSATz1 aSATz2 aSATz3
#  aguila --scenarios='{1}' --timesteps=[2,300,2] aFCz0 aFCz1 aFCz2 aFCz3

def reportKsatEvolution(model, ksat_list):
    for i in range(len(ksat_list)):
        name = str('Ksatz') + str(i)
        model.report(ksat_list[i], name)

#  aguila --scenarios='{1}' --timesteps=[2,300,2] Ksatz0 Ksatz1 Ksatz2 Ksatz3

def report_CN(model, CN2, cumulative_rain_mm):
    model.report(CN2, 'CN2')
    model.report(cumulative_rain_mm, 'CumPmm')

#  aguila --scenarios='{1}' --timesteps=[2,300,2] CN2 CumPmm landuse2016

def checkMoisture(model, moisture_maps, name):
    for layer in range(len(moisture_maps)):
        model.report(moisture_maps[layer], name + str(layer))

#  aguila --scenarios='{1}' --timesteps=[1,300,1] athz0 athz1 athz2 athz3
# aguila d14_theta_z0

# Water balance and reporting
def recordInfiltration(model, infiltration, layer):
    model.water_balance[layer] += infiltration
    model.report(infiltration, 'aINFz' + str(layer))

#  aguila --scenarios='{1}' --timesteps=[2,300,2] aINFz0 aINFz1 aINFz2 aINFz3
#  aguila --scenarios='{1}' --timesteps=[1,300,1] aROm3 outlet_multi_nom_v3.map


def recordPercolation(model, percolation, layer):
    model.water_balance[layer] -= percolation
    model.report(percolation, 'aPERz' + str(layer))

#  aguila --scenarios='{1}' --timesteps=[2,300,2] aPERz0 aPERz1 aPERz2 aPERz3


def recordRunOff(model, runoff, unit='m3'):
    if unit == 'm3':
        ro_m3 = runoff * cellarea() / 1000  # m3
        accu_ro_m3 = accuflux(model.ldd_subs, ro_m3)
        model.report(accu_ro_m3, 'aROm3')
        model.report(areatotal(accu_ro_m3, model.outlet_multi), 'totROm3')
    else:
        model.report(runoff, 'aROmm')


# aguila --scenarios='{1}' --timesteps=[2,300,2] aROm3
# aguila --scenarios='{1}' --timesteps=[2,300,2] aROmm
