from pcraster import *
from pcraster._pcraster import *
from pcraster.framework import *
import os
import sys

print(os.getcwd())



dem = readmap("demslope")
ldd_subs = readmap('ldd')
demldd = readmap('demldd')
theta = readmap("thetiZ0")

out_ditch_ldd = lddcreate(demldd, 1E35, 1E35, 1E35, 1E35)
""" Error: out_ditch_ldd is not exactly equivalent as ldd 

Doubt: 
demldd was created to route into the ditch, 
"""
outlet = readmap("outlet")

aguila(outlet, out_ditch_ldd)