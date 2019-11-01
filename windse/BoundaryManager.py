""" 
The BoundaryManager submodule contains the classes required for 
defining the boundary conditions. The boundaries need to be numbered
as follows:

    * 1: Top
    * 2: Bottom
    * 3: Front
    * 4: Back
    * 5: Left
    * 6: Right

Currently, the inflow direction is the Left boundary.
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *  

    import math 
    from scipy.interpolate import RegularGridInterpolator


class GenericBoundary(object):
    def __init__(self,dom,fs,farm):
        self.params = windse_parameters
        self.dom = dom
        self.fs = fs
        self.farm = farm
        self.ig_first_save = True
        self.height_first_save = True
        self.fprint = self.params.fprint

        ### Check if boundary information is given in the params file ###
        bcs_params = self.params.get("boundary_condition",{})
        self.boundary_names = bcs_params.get("boundary_names", dom.boundary_names)
        self.boundary_types = bcs_params.get("boundary_types", dom.boundary_types)
        self.HH_vel = bcs_params.get("HH_vel", 8.0)
        self.power = bcs_params.get("power", 0.25)
        self.k = bcs_params.get("k", 0.4)
        self.init_wind = dom.init_wind

        ### Define the zero function based on domain dimension ###
        self.zeros = Constant(dom.mesh.topology().dim()*(0.0,))
        self.zero  = Constant(0.0)

    def SetupBoundaries(self):
        ### Create the equations need for defining the boundary conditions ###
        ### this is sloppy and will be cleaned up.
        ### Inflow is always from the front

        self.fprint("Applying Boundary Conditions",offset=1)

        ### Assemble boundary conditions ###
        bcu_eqns = []
        bcp_eqns = []
        for bc_type, bs in self.boundary_types.items():
            if bc_type == "inflow":
                for b in bs:
                    bcu_eqns.append([self.fs.V, self.fs.W.sub(0), self.bc_velocity, self.boundary_names[b]])

            elif bc_type == "no_slip":
                for b in bs:
                    bcu_eqns.append([self.fs.V, self.fs.W.sub(0), self.zeros, self.boundary_names[b]])

            elif bc_type == "free_slip":
                temp_list = list(self.boundary_names.keys()) # get ordered list
                for b in bs:
                    dim = math.floor(temp_list.index(b)/2.0) # get dim based on order
                    bcu_eqns.append([self.fs.V.sub(dim), self.fs.W.sub(0).sub(dim), self.zero, self.boundary_names[b]])

            elif bc_type == "no_stress":
                for b in bs:
                    bcu_eqns.append([None, None, None, self.boundary_names[b]])
                    bcp_eqns.append([self.fs.Q, self.fs.W.sub(1), self.zero, self.boundary_names[b]])

            else:
                raise ValueError(bc_type+" is not a recognized boundary type")
        bcs_eqns = bcu_eqns#+bcp_eqns

        ### Set the boundary conditions ###
        self.bcu = []
        for i in range(len(bcu_eqns)):
            if bcu_eqns[i][0] is not None:
                self.bcu.append(DirichletBC(bcu_eqns[i][0], bcu_eqns[i][2], self.dom.boundary_markers, bcu_eqns[i][3]))

        self.bcp = []
        for i in range(len(bcp_eqns)):
            if bcp_eqns[i][0] is not None:
                self.bcp.append(DirichletBC(bcp_eqns[i][0], bcp_eqns[i][2], self.dom.boundary_markers, bcp_eqns[i][3]))

        self.bcs = []
        for i in range(len(bcs_eqns)):
            if bcs_eqns[i][0] is not None:
                self.bcs.append(DirichletBC(bcs_eqns[i][1], bcs_eqns[i][2], self.dom.boundary_markers, bcs_eqns[i][3]))

        self.fprint("Boundary Conditions Applied",offset=1)
        self.fprint("")

    def PrepareVelocity(self,theta):
        length = len(self.unit_reference_velocity)
        ux_com = np.zeros(length)
        uy_com = np.zeros(length)
        uz_com = np.zeros(length)

        for i in range(length):
            v = self.HH_vel * self.unit_reference_velocity[i]
            ux_com[i] = math.cos(theta)*v
            uy_com[i] = math.sin(theta)*v
            if self.dom.dim == 3:
                uz_com[i] = 0.0   
        return [ux_com,uy_com,uz_com]

    def RecomputeVelocity(self,theta):
        self.fprint("Recomputing Velocity")
        ux_com, uy_com, uz_com = self.PrepareVelocity(theta)

        self.ux = Function(self.fs.V0)
        self.uy = Function(self.fs.V1)
        if self.dom.dim == 3:
            self.uz = Function(self.fs.V2)

        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        if self.dom.dim == 3:
            self.uz.vector()[:] = uz_com

        ### Assigning Velocity
        self.bc_velocity = Function(self.fs.V)
        if self.dom.dim == 3:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        else:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy])
        
        ### Create Pressure Boundary Function
        self.bc_pressure = Function(self.fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(self.fs.W)
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        self.SetupBoundaries()

    def UpdateVelocity(self, simTime):
        pass

    def SaveInitialGuess(self,val=0):
        """
        This function saves the turbine force if exists to output/.../functions/
        """
        if self.ig_first_save:
            self.u0_file = self.params.Save(self.bc_velocity,"u0",subfolder="functions/",val=val)
            self.p0_file = self.params.Save(self.bc_pressure,"p0",subfolder="functions/",val=val)
            self.ig_first_save = False
        else:
            self.params.Save(self.bc_velocity,"u0",subfolder="functions/",val=val,file=self.u0_file)
            self.params.Save(self.bc_pressure,"p0",subfolder="functions/",val=val,file=self.p0_file)

    def SaveHeight(self,val=0):
        """
        This function saves the turbine force if exists to output/.../functions/
        """
        if self.height_first_save:
            self.height_file = self.params.Save(self.height,"height",subfolder="functions/",val=val)
            self.depth_file = self.params.Save(self.depth,"depth",subfolder="functions/",val=val)
            self.height_first_save = False
        else:
            self.params.Save(self.height,"height",subfolder="functions/",val=val,file=self.height_file)
            self.params.Save(self.depth,"depth",subfolder="functions/",val=val,file=self.depth_file)

    def CalculateHeights(self):
        ### Calculate the distance to the ground for the Q function space ###
        # self.z_dist_Q = Function(fs.Q)
        self.height = Function(self.fs.Q)
        self.depth = Function(self.fs.Q)
        Q_coords = self.fs.Q.tabulate_dof_coordinates()
        height_vals = self.height.vector()[:]
        for i in range(len(Q_coords)):
            height_vals[i] = self.dom.Ground(Q_coords[i,0],Q_coords[i,1])
        z_dist_Q = Q_coords[:,2]-height_vals
        self.height.vector()[:]=height_vals
        self.depth.vector()[:]=z_dist_Q

        ### Calculate the distance to the ground for the V function space ###
        self.depth_V = Function(self.fs.V)
        V_coords = self.fs.V.tabulate_dof_coordinates()
        z_dist_V_val = np.zeros(len(V_coords))
        for i in range(len(V_coords)):
            z_dist_V_val[i] = V_coords[i,2]-self.dom.Ground(V_coords[i,0],V_coords[i,1])
        self.depth_V.vector()[:]=z_dist_V_val

        self.V0_coords = self.fs.V0.tabulate_dof_coordinates()


class UniformInflow(GenericBoundary):
    def __init__(self,dom,fs,farm):
        super(UniformInflow, self).__init__(dom,fs,farm)
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Uniform Inflow")
        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        ### Create the Velocity Function ###
        self.ux = Function(fs.V0)
        self.uy = Function(fs.V1)
        if self.dom.dim == 3:
            self.uz = Function(fs.V2)
        self.unit_reference_velocity = np.full(len(self.ux.vector()[:]),1.0)
        self.ux.vector()[:] = self.unit_reference_velocity

        ux_com, uy_com, uz_com = self.PrepareVelocity(self.init_wind)
        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        if self.dom.dim == 3:
            self.uz.vector()[:] = uz_com

        ### Compute distances ###
        if self.dom.dim == 3:
            self.fprint("Computing Distance to Ground")
            self.CalculateHeights()

        ### Assigning Velocity
        self.fprint("Computing Velocity Vector")
        self.bc_velocity = Function(fs.V)
        if self.dom.dim == 3:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        else:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(fs.W)
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Finished",special="footer")

class PowerInflow(GenericBoundary):
    """
    PowerInflow creates a set of boundary conditions where the x-component
    of velocity follows a power law. Currently the function is 

    .. math::

        u_x=8.0 \\left( \\frac{z-z_0}{z_1-z_0} \\right)^{0.15}.
        
    where :math:`z_0` is the ground and :math:`z_1` is the top of the domain.

    Args:
        dom (:class:`windse.DomainManager.GenericDomain`): A windse domain object.
        fs (:class:`windse.FunctionSpaceManager.GenericFunctionSpace`): 
            A windse function space object

    Todo:
        * Make the max velocity an input
        * Make the power an input
    """
    def __init__(self,dom,fs,farm):
        super(PowerInflow, self).__init__(dom,fs,farm)

        if self.dom.dim != 3:
            raise ValueError("PowerInflow can only be used with 3D domains.")

        ### Setup Boundary Conditions
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Power Law Inflow")
        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        self.fprint("")

        ### Compute distances ###
        self.fprint("Computing Distance to Ground")
        self.CalculateHeights()
        depth_v0,depth_v1,depth_v2 = self.depth_V.split(deepcopy=True)

        ### Create the Velocity Function ###
        self.fprint("Computing Velocity Vector")
        self.ux = Function(fs.V0)
        self.uy = Function(fs.V1)
        self.uz = Function(fs.V2)
        
        #################
        #################
        #################
        #################
        #################
        #################
        scaled_depth = np.abs(np.divide(depth_v0.vector()[:],(np.mean(farm.HH)-dom.z_range[0])))
        # scaled_depth = np.abs(np.divide(depth_v0.vector()[:],(np.mean(farm.HH)-0.0)))
        #################
        #################
        #################
        #################
        #################

        self.unit_reference_velocity = np.power(scaled_depth,self.power)
        # self.reference_velocity = np.multiply(self.HH_vel,np.power(scaled_depth,self.power))
        ux_com, uy_com, uz_com = self.PrepareVelocity(self.init_wind)

        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        self.uz.vector()[:] = uz_com

        ### Assigning Velocity
        self.bc_velocity = Function(self.fs.V)
        if self.dom.dim == 3:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        else:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(self.fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(self.fs.W)
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Setup",special="footer")

class LogLayerInflow(GenericBoundary):
    def __init__(self,dom,fs,farm):
        super(LogLayerInflow, self).__init__(dom,fs,farm)

        if self.dom.dim != 3:
            raise ValueError("LogLayerInflow can only be used with 3D domains.")

        ### Setup Boundary Conditions
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Power Law Inflow")

        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        self.fprint("")

        ### Compute distances ###
        self.fprint("Computing Distance to Ground")
        self.CalculateHeights()
        depth_v0,depth_v1,depth_v2 = self.depth_V.split(deepcopy=True)



        ### Create the Velocity Function ###
        self.fprint("Computing Velocity Vector")
        self.ux = Function(fs.V0)
        self.uy = Function(fs.V1)
        self.uz = Function(fs.V2)
        if dom.z_range[0] == 0:
            scaled_depth = np.abs(np.divide(depth_v0.vector()[:]+0.01,0.01))
            ustar = self.k/np.log(np.mean(farm.HH)/0.01)
        elif dom.z_range[0] <= 0:
            raise ValueError("Log profile cannot be used with negative z values")
        else:
            scaled_depth = np.abs(np.divide(depth_v0.vector()[:]+dom.z_range[0],(dom.z_range[0])))
            ustar = self.k/np.log(np.mean(farm.HH)/dom.z_range[0])
        self.unit_reference_velocity = np.multiply(ustar/self.k,np.log(scaled_depth))
        ux_com, uy_com, uz_com = self.PrepareVelocity(self.init_wind)

        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        self.uz.vector()[:] = uz_com

        ### Assigning Velocity
        self.bc_velocity = Function(self.fs.V)
        if self.dom.dim == 3:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        else:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(self.fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(self.fs.W)
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Setup",special="footer")

class TurbSimInflow(LogLayerInflow):
    def __init__(self,dom,fs,farm):
        super(TurbSimInflow, self).__init__(dom,fs,farm)

        ### Get the path for turbsim data ###
        self.turbsim_path = self.params["boundary_condition"].get("turbsim_path",None)
        if self.turbsim_path is None:
            raise ValueError("Please provide the path to the turbsim data")
        self.tFinal = self.params["solver"].get("final_time",1)

        ### Load Turbsim Data ###
        uTotal = np.load(self.turbsim_path+'turb_u.npy')
        vTotal = np.load(self.turbsim_path+'turb_v.npy')
        wTotal = np.load(self.turbsim_path+'turb_w.npy')

        ### Extract number of data points ###
        ny = np.shape(uTotal)[1]
        nz = np.shape(uTotal)[0]
        nt = np.shape(uTotal)[2]

        ### Create the data bounds ###
        y = np.linspace(self.dom.y_range[0], self.dom.y_range[1], ny)
        z = np.linspace(self.dom.z_range[0], self.dom.z_range[1], nz)
        t = np.linspace(0.0, self.tFinal, nt)

        ### Build interpolating functions ###
        self.interp_u = RegularGridInterpolator((z, y, t), uTotal)
        self.interp_v = RegularGridInterpolator((z, y, t), vTotal)
        self.interp_w = RegularGridInterpolator((z, y, t), wTotal)

        ### Locate Boundary DOFS indexes ###
        # Define tolerance
        tol = 1e-6

        ##### FIX MAKE WORK FOR ALL BOUNDARY INFLOW ####
        # Iterate and fine the boundary IDs
        self.boundaryIDs = []
        for k, pos in enumerate(self.V0_coords):
            if pos[0] < self.dom.x_range[0] + tol:
                self.boundaryIDs.append(k)

        self.UpdateVelocity(0.0)

    def UpdateVelocity(self, simTime):

        # Define tolerance
        tol = 1e-6

        # Interpolate a value at each boundary coordinate
        for k in self.boundaryIDs:
            # Get the position corresponding to this boundary id
            pos = self.V0_coords[k, :]

            # The interpolation point specifies a 3D (z, y, time) point
            xi = np.array([pos[2], pos[1], simTime])

            # Get the interpolated value at this point
            self.ux.vector()[k] = self.interp_u(xi)
            self.uy.vector()[k] = self.interp_v(xi)
            self.uz.vector()[k] = self.interp_w(xi)

        ### Assigning Velocity
        self.bc_velocity = Function(self.fs.V)
        if self.dom.dim == 3:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        else:
            self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy])

        self.SetupBoundaries()












