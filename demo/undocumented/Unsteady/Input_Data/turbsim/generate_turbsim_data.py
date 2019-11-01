"""
This script provides an example usage of the PyTurbSim API.
"""

# Begin by importing the PyTurbSim API:
import pyts.api as pyts_api
import pyts.io.main as pyts
import numpy as np

# Define some variables for the PyTurbSim run:
refht = 100
ustar = .15
Uref = 8

# First we initialize a PyTurbSim 'run' object:
tsr = pyts_api.tsrun()

# Next we give this run object a grid:
tsr.grid = pyts_api.tsGrid(
    center = refht, ny = 15, nz = 15, height = 200, width = 160, time_sec = 1000, dt = 0.5)

# Now we define a mean 'profile model',
prof_model = pyts_api.profModels.h2l(Uref, refht, ustar)
# and assign it to the run object,
tsr.prof = prof_model
# These two steps can be completed in one as:
#tsr.profModel=pyts.profModels.h2l(U,refht,ustar)

# Next we define and assign a 'spectral model' to the run object,
tsr.spec = pyts_api.specModels.tidal(ustar, refht)

# ... and define/assign a 'coherence model',
tsr.cohere = pyts_api.cohereModels.nwtc()

# ... and define/assign a 'stress model',
tsr.stress = pyts_api.stressModels.tidal(ustar, refht)

# Now simply 'call' the run oject to produce the TurbSim output.
turbsim_output = tsr()

# We can save the output in 'bladed' format,
turbsim_output.write_bladed('ExampleOutput.bl')


a = pyts.readModel('ExampleOutput.wnd')

fn_u = 'turb_u.npy'
fn_v = 'turb_v.npy'
fn_w = 'turb_w.npy'
fn_y = 'turb_y.npy'
fn_z = 'turb_z.npy'

turb_u_fp = open(fn_u, 'w+')
np.save(turb_u_fp, a.u)
turb_u_fp.close()

turb_v_fp = open(fn_v, 'w+')
np.save(turb_v_fp, a.v)
turb_v_fp.close()

turb_w_fp = open(fn_w, 'w+')
np.save(turb_w_fp, a.w)
turb_w_fp.close()


turb_y_fp = open(fn_y, 'w+')
np.save(turb_y_fp, a.y)
turb_y_fp.close()


turb_z_fp = open(fn_z, 'w+')
np.save(turb_z_fp, a.z)
turb_z_fp.close()
