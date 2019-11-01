"""
This script provides an example usage of the PyTurbSim API.
"""

# Begin by importing the PyTurbSim API:
import pyts.api as pyts
import numpy as np

# Define some variables for the PyTurbSim run:
refht = 10
ustar = 0.03
Uref = 8

# First we initialize a PyTurbSim 'run' object:
tsr = pyts.tsrun()

# Next we give this run object a grid:
tsr.grid = pyts.tsGrid(
    center = refht, ny = 30, nz = 5, height = 500.0, width = 3000.0, time_sec = 6000, dt = 20.0)

# Now we define a mean 'profile model',
prof_model = pyts.profModels.h2l(Uref, refht, ustar)
# and assign it to the run object,
tsr.prof = prof_model
# These two steps can be completed in one as:
#tsr.profModel=pyts.profModels.h2l(U,refht,ustar)

# Next we define and assign a 'spectral model' to the run object,
tsr.spec = pyts.specModels.tidal(ustar, refht)

# ... and define/assign a 'coherence model',
tsr.cohere = pyts.cohereModels.nwtc()

# ... and define/assign a 'stress model',
tsr.stress = pyts.stressModels.tidal(ustar, refht)

# Now simply 'call' the run oject to produce the TurbSim output.
turbsim_output = tsr()

print(type(turbsim_output))
print(turbsim_output.shape)


# # We can save the output in 'bladed' format,
# turbsim_output.write_bladed('ExampleOutput.bl')


# a = pyts.readModel('ExampleOutput.wnd')

# fn_u = '../windse/pyturbsim_outputs/turb_u.npy'
# fn_v = '../windse/pyturbsim_outputs/turb_v.npy'
# fn_w = '../windse/pyturbsim_outputs/turb_w.npy'
# fn_y = '../windse/pyturbsim_outputs/turb_y.npy'
# fn_z = '../windse/pyturbsim_outputs/turb_z.npy'

# turb_u_fp = open(fn_u, 'w+')
# np.save(turb_u_fp, a.u)
# turb_u_fp.close()

# turb_v_fp = open(fn_v, 'w+')
# np.save(turb_v_fp, a.v)
# turb_v_fp.close()

# turb_w_fp = open(fn_w, 'w+')
# np.save(turb_w_fp, a.w)
# turb_w_fp.close()


# turb_y_fp = open(fn_y, 'w+')
# np.save(turb_y_fp, a.y)
# turb_y_fp.close()


# turb_z_fp = open(fn_z, 'w+')
# np.save(turb_z_fp, a.z)
# turb_z_fp.close()
