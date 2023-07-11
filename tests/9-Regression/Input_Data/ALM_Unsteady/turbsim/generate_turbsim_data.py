"""
This script provides an example usage of the PyTurbSim API.
"""

# Begin by importing the PyTurbSim API:
import pyts.api as pyts_api
import pyts.io.main as pyts
import numpy as np

# Domain Params:
H   = 200  # Height of the domain (m)
L   = 160  # Width of the domain (m)
ny  = 15   # Grid point is the y direction (<30)
nz  = 15   # Grid point is the z direction (<30)
T   = 1000 # Final Time (s)
dt  = 0.5  # Time step size (s)
tol = 0.1  # Buffer between ground and grid (>0)

# Physics Params:
Zref = 32.1 # Hub height (m)
Uref = 8.0  # Velocity at hub height (m/s)
Z0   = 0.04 # Roughness height (m)
k    = 0.4  # von Karman constant (-)
Ri   = 0.0  # Richardson Number stability parameter (-)

# Calculate the friction velocity
Ustar = Uref*k/np.log(Zref/Z0)

# First we initialize a PyTurbSim 'run' object:
tsr = pyts_api.tsrun()

# Next we give this run object a grid:
tsr.grid = pyts_api.tsGrid(
    center = H/2.0+tol, ny = ny, nz = nz, height = H, width = L, time_sec = T, dt = dt)

# Now we define a mean 'profile model',
prof_model = pyts_api.profModels.h2l(Uref, Zref, Ustar)
# and assign it to the run object,
tsr.prof = prof_model
# These two steps can be completed in one as:
#tsr.profModel=pyts.profModels.h2l(U,refht,ustar)

# Next we define and assign a 'spectral model' to the run object,
tsr.spec = pyts_api.specModels.tidal(Ustar, H/2.0)

# ... and define/assign a 'coherence model',
tsr.cohere = pyts_api.cohereModels.nwtc()

# ... and define/assign a 'stress model',
tsr.stress = pyts_api.stressModels.tidal(Ustar, H/2.0)

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






import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt

mag = np.sqrt(np.power(a.u,2.0)+np.power(a.v,2.0),np.power(a.w,2.0))

p = plt.contourf(a.y,a.z,np.mean(mag,2),32)
plt.colorbar(p)
plt.show()

p = plt.contour(a.y,a.z,mag[:,:,-25],32)
plt.colorbar(p)
plt.show()
