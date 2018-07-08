# -*- coding: utf-8 -*-

# from time import *
import time
from datetime import datetime


from pcraster._pcraster import *
from pcraster.framework import *
import os


class MyFirstModel(DynamicModel):

    def __init__(self, cloneMap):
        DynamicModel.__init__(self)
        setclone(cloneMap)

    def initial(self):
        self.landuse = self.readmap("fields_cover")
        print 'running the initial'

    def dynamic(self):
        timeStep = self.currentTimeStep()
        setglobaloption('matrixtable')
        fields = timeinputscalar('landuse.tss', nominal(self.landuse))
        crop_type = lookupscalar('croptable.tbl', 2, fields)

        ksat1Beet = timeinputscalar('rain.tss', nominal(self.landuse))
        ksat1Corn = timeinputscalar('ET0.tss', nominal(self.landuse))
        ksat = timeinputscalar('ksats.tss', nominal(self.landuse))

        if timeStep == 3:
            aguila(self.landuse, ksat)



nrOfTimeSteps=5

myModel = MyFirstModel("clone_nom.map") # an instance of the model, which inherits from class: DynamicModel
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps) # an instance of the Dynamic Framework
dynamicModel.run()

  




