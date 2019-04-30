from dolfin import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6

### Initialize WindSE ###
windse.initialize("params.yaml")

### Generate Domain ###
dom = windse.BoxDomain()
dom.Save()

### Generate Wind Farm ###
farm = windse.ImportedWindFarm(dom)
farm.Plot()

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

### Save the Turbine Force ##
# farm.RotateFarm(3.14159/4.0)
u = farm.TurbineForce(fs,dom.mesh)

filename = windse.windse_parameters.folder+"/functions/turbine_force.pvd"
File(filename) << farm.tf