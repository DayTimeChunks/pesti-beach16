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

def checkMoisture(model, moisture_maps, name):
    for layer in range(len(moisture_maps)):
        model.report(moisture_maps[layer], name + str(layer))

#  aguila --scenarios='{1}' --timesteps=[2,300,2] athz0 athz1 athz2 athz3


# Water balance and reporting
def recordInfiltration(model, infiltration, layer):
    model.water_balance[layer] += infiltration
    model.report(infiltration, 'aINFz' + str(layer))

#  aguila --scenarios='{1}' --timesteps=[2,300,2] aINFz0 aINFz1 aINFz2 aINFz3
#  aguila --scenarios='{1}' --timesteps=[2,300,2] aROm3

def recordPercolation(model, percolation, layer):
    model.water_balance[layer] -= percolation
    model.report(percolation, 'aPERz' + str(layer))

#  aguila --scenarios='{1}' --timesteps=[2,300,2] aPERz0 aPERz1 aPERz2 aPERz3


def recordRunOff(model, runoff, unit='m3'):
    if unit == 'mm':
        model.report(runoff, 'aROmm')
    else:
        ro_m3 = runoff * cellarea() / 1000  # m3
        model.report(ro_m3, 'aROm3')

# aguila --scenarios='{1}' --timesteps=[2,300,2] aROm3
# aguila --scenarios='{1}' --timesteps=[2,300,2] aROmm