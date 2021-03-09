from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import glob

from alm_class import UpdateActuatorLineForce
from simpler_alm_class import SimplerUpdateActuatorLineForce
import openmdao.api as om
from full_alm_openmdao import FullALM
from alm_group import ALMGroup

# Create a 3D box mesh from (-100, -100, 0) to (100, 100, 200)
# with RES x RES x RES total nodes
res = 3

mesh = BoxMesh(Point(-100, -100, 0),
              Point(100, 100, 200),
              res, res, res)

# Create a vector function space (used by fluid velocity and turbine force)
V = VectorFunctionSpace(mesh, 'P', 2)


# ================================================================
# Create a problem object (normally from ProblemManager) and the nested objects it inherits
# that represent what a WindSE run would normally look like at the point of calling
# the actuator line functions
# ================================================================
class BlankObject:
    pass

problem = BlankObject()


# ================================================================
# dom, from DomainManager
# ================================================================
problem.dom = BlankObject()
# ndim, the dimension of the mesh, either 2 or 3
problem.dom.dim = mesh.geometry().dim()


# ================================================================
# farm, from WindFarmManager
# ================================================================
problem.farm = BlankObject()
# radius, a list with num_turbs elements, radius[0] = blade radius of turbine 0
problem.farm.radius = [65.0]
# myaw, a list with num_turbs elements, myaw[0] = yaw angle of turbine 0
yaw = 0.
problem.farm.myaw = [yaw]
# mx, a list with num_turbs elements, mx[0] = x-position turbine 0
problem.farm.mx = [0]
problem.farm.x = [0]
# my, a list with num_turbs elements, my[0] = y-position turbine 0
problem.farm.my = [0]
problem.farm.y = [0]
# z, a list with num_turbs elements, z[0] = z-position turbine 0 (hub height, if ground level = 0)
problem.farm.z = [110.0]


# ================================================================
# fs, from FunctionSpaceManager
# ================================================================
problem.fs = BlankObject()
# V, a vector function space associated with the mesh
problem.fs.V = V


# ================================================================
# Variable set directly in ProblemManager
# ================================================================
# gaussian_width, the characteristic size of the gaussian, big enough to ensure it intersects mesh points,
# small enough to ensure a relatively "sharp" application of the blade force
hmin = mesh.hmin()/np.sqrt(3.0)
problem.gaussian_width = 3.0*hmin
# num_blade_segments, the number of actuator nodes (separate gaussian distributions), used to represent 
# the blade from root to tip
problem.num_blade_segments = 4

def build_lift_and_drag_tables(airfoil_data_path):

    # Determine the number of files in airfoil_data_path
    num_files = len(glob.glob('%s/*.txt' % (airfoil_data_path)))

    for file_id in range(num_files):
        # print('Reading Airfoil Data #%d' % (file_id))
        data = np.genfromtxt('%s/af_station_%d.txt' % (airfoil_data_path, file_id), skip_header=1, delimiter=' ')

        if file_id == 0:
            # If this is the first file, store the angle data
            interp_angles = data[:, 0]
            num_angles = np.size(interp_angles)

            # If this is the first file, allocate space for the tables        
            lift_table = np.zeros((num_angles, num_files))
            drag_table = np.zeros((num_angles, num_files))

        # Store all the lift and drag data in the file_id column
        lift_table[:, file_id] = data[:, 1]
        drag_table[:, file_id] = data[:, 2]

    return lift_table, drag_table, interp_angles

# Get the lift and drag values that will later be used in interpolation
lift_table, drag_table, interp_angles = build_lift_and_drag_tables('airfoil_polars')
problem.lift_table = lift_table
problem.drag_table = drag_table
problem.interp_angles = interp_angles

# Get the chord, twist, cL, and cD values from IEA RWT data
actual_turbine_data = np.genfromtxt('Input_Data/baseline.csv', delimiter = ',', skip_header = 1)

actual_x = actual_turbine_data[:, 0]
actual_chord = actual_turbine_data[:, 1]
actual_twist = actual_turbine_data[:, 2]/180.0*np.pi
actual_cl = actual_turbine_data[:, 3]
actual_cd = actual_turbine_data[:, 4]

chord_interp = interp.interp1d(actual_x, actual_chord)
twist_interp = interp.interp1d(actual_x, actual_twist)
cl_interp = interp.interp1d(actual_x, actual_cl)
cd_interp = interp.interp1d(actual_x, actual_cd)

interp_points = np.linspace(0.0, 1.0, problem.num_blade_segments)

# Generate the interpolated values
chord = chord_interp(interp_points)
twist = twist_interp(interp_points)
cl = cl_interp(interp_points)
cd = cd_interp(interp_points)


# mtwist, the twist of the sections along the blade length
problem.mtwist = [twist]
# mcl, the coefficient of lift along the blade
problem.mcl = [cl]
# mcl, the coefficient of drag along the blade
problem.mcd = [cd]
# mcl, the chord length along the blade
problem.mchord = [chord]

# rpm, the rotational speed of the blade
problem.rpm = 10.0
# The coordinates at each point of the mesh, indexed using the same ordering as the vector function space
coords = problem.fs.V.tabulate_dof_coordinates()
coords = np.copy(coords[0::problem.dom.dim, :])
problem.coords = coords
problem.coordsLinear = np.copy(coords.reshape(-1, 1))
# min_dist, a cutoff distance, if the min_distance < threshold, this section of mesh is close enough to
# the current turbine that the forcing effect is non-negligible and must be counted
problem.min_dist = [0.0]
# rotor_torque, a list with num_turbs elements that holds the rotor torque, each gaussian contributes a bit of force
problem.rotor_torque = [0]
problem.rotor_torque_count = [0]
# cyld_expr_list, a list with num_turbs elements where each entry is a cylinder 
# encapsulating the turbine (used for later calculting the power) 
problem.cyld_expr_list = [0]



# ================================================================
# Variables passed in directly as part of looping through all turbines and timesteps
# ================================================================

# Create a dummy fluid velocity field
u_local = Function(problem.fs.V)

# Set the values of the dummy velocity field to be a linear shear in the +x direction
inflow_profile = Expression(('9.0', '0', '0'), degree=2)
# inflow_profile = Expression(('cos(x[0])*sin(x[1])', 'cos(x[1])*2.2', 'sin(x[0])'), degree=6)
u_local.interpolate(inflow_profile)

# Currently, only a single turbine
turb_i = 0

# Not currently being used, but left in for easily merging your notebook changes to WindSE
mpi_u_fluid = 1

# The number of timesteps to simulate
tSteps = 1

# The timestep size
dt = 0.1

# The accumulated simulation time
simTime = 0.0

# simTime_list, a list containing all simulation timesteps accumulated
problem.simTime_list = []

# Initialize file pointers for writing
fp_vel = File('output/velocity.pvd')
fp_tf = File('output/turbine_force.pvd')

# ================================================================
# Call the actuator line function and return the turbine force being projected on the computational grid
# ================================================================

for k in range(tSteps):
    print('Timestep %d of %d, t = %.2f' % (k+1, tSteps, simTime))
    
    # Append the simTime list and set a convenience variable
    problem.simTime_list.append(simTime)
    simTime_id = k
    
    num_blades = 1
    prob = om.Problem()
    prob.model.add_subsystem('ALMGroup', ALMGroup(problem=problem,
        simTime_id=simTime_id,
        dt=dt,
        turb_i=turb_i,
        num_blades=num_blades,
        ), promotes=['*'])
    prob.setup()
    prob['u_local'] = [9., 0., 0.]
    prob['yaw'] = yaw
    prob.run_model()
    om_forces = prob['turbine_forces']
    
    totals = prob.compute_totals('turbine_forces', 'u_local')
    print(totals)
    
    # check_partials_data = prob.check_partials(compact_print=True, includes='*ComputeLiftDragForces')
    # 
    # om.partial_deriv_plot('lift_force', 'u_unit_vec', check_partials_data, binary=False)
    # om.partial_deriv_plot('lift_force', 'blade_unit_vec', check_partials_data, binary=False)
    