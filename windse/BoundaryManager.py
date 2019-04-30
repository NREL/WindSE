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

class GenericBoundary(object):
    def __init__(self,dom,fs):
        self.params = windse_parameters
        self.dom = dom
        self.fs = fs
        self.ig_first_save = True
        self.height_first_save = True
        self.fprint = self.params.fprint

        ### Check if boundary information is given in the params file ###
        bcs_params = self.params.get("boundary_conditions",{})
        self.boundary_names = bcs_params.get("names", dom.boundary_names)
        self.boundary_types = bcs_params.get("types", dom.boundary_types)
        self.wind_direction = dom.wind_direction

        ### Define the zero function based on domain dimension ###
        self.zeros = Constant(dom.mesh.topology().dim()*(0.0,))

    def SetupBoundaries(self):
        ### Create the equations need for defining the boundary conditions ###
        ### this is sloppy and will be cleaned up.
        ### Inflow is always from the front

        self.fprint("Applying Boundary Conditions",offset=1)

        ### Assemble boundary conditions ###
        bcs_eqns = []
        for bc_type, bs in self.boundary_types.items():
            if bc_type == "inflow":
                for b in bs:
                    bcs_eqns.append([self.bc_velocity,self.boundary_names[b]])

            elif bc_type == "no_slip":
                for b in bs:
                    bcs_eqns.append([self.zeros,self.boundary_names[b]])


            elif bc_type == "no_stress":
                for b in bs:
                    bcs_eqns.append([None,self.boundary_names[b]])

            else:
                raise ValueError(bc_type+"is not a recognized boundary type")


        ### Set the boundary conditions ###
        self.bcs = []
        for i in range(len(bcs_eqns)):
            if bcs_eqns[i][0] is not None:
                self.bcs.append(DirichletBC(self.fs.W.sub(0), bcs_eqns[i][0], self.dom.boundary_markers, bcs_eqns[i][1]))
        self.fprint("Boundary Conditions Applied",offset=1)
        self.fprint("")

    # def RotateVelocity(self,ux,uy,uz,theta):
    #     ux_hat = ux.vector()[:]
    #     uy_hat = uy.vector()[:]
    #     uz_hat = uz.vector()[:]

    #     for i in range(len(ux_hat)):
    #         x = ux_hat[i]
    #         y = uy_hat[i]
    #         z = uz_hat[i]
    #         ux_hat[i] = math.cos(theta)*(x) - math.sin(theta)*(y)
    #         uy_hat[i] = math.sin(theta)*(x) + math.cos(theta)*(y)
    #         uz_hat[i] = z

    #     ux.vector()[:] = ux_hat
    #     uy.vector()[:] = uy_hat
    #     uz.vector()[:] = uz_hat    
    #     return [ux,uy,uz]

    def RotateVelocity(self,theta):
        length = len(self.reference_velocity)
        ux_com = np.zeros(length)
        uy_com = np.zeros(length)
        uz_com = np.zeros(length)

        for i in range(length):
            v = self.reference_velocity[i]
            ux_com[i] = math.cos(theta)*v
            uy_com[i] = math.sin(theta)*v
            uz_com[i] = 0.0   
        return [ux_com,uy_com,uz_com]

    def RecomputeVelocity(self,theta):
        self.fprint("Recomputing Velocity")
        ux_com, uy_com, uz_com = self.RotateVelocity(theta)
        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        self.uz.vector()[:] = uz_com

        ### Assigning Velocity
        self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])
        
        self.SetupBoundaries()

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

class UniformInflow(GenericBoundary):
    def __init__(self,dom,fs):
        super(UniformInflow, self).__init__(dom,fs)
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Uniform Inflow")
        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        ### Create the Velocity Function ###
        ux = Function(fs.V0)
        uy = Function(fs.V1)
        uz = Function(fs.V2)
        ux.vector()[:] = np.full(len(ux.vector()[:]),8)

        ### Assigning Velocity
        self.fprint("Computing Velocity Vector")
        self.bc_velocity = Function(fs.V)
        fs.VelocityAssigner.assign(self.bc_velocity,[ux,uy,uz])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(fs.W)
        fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Finished",special="footer")

class UniformInflow2D(GenericBoundary):
    def __init__(self,dom,fs):
        super(UniformInflow2D, self).__init__(dom,fs)
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Uniform Inflow")
        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        ### Create the Velocity Function ###
        ux = Function(fs.V0)
        uy = Function(fs.V1)
        ux.vector()[:] = np.full(len(ux.vector()[:]),8)

        ### Assigning Velocity
        self.fprint("Computing Velocity Vector")
        self.bc_velocity = Function(fs.V)
        fs.VelocityAssigner.assign(self.bc_velocity,[ux,uy])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(fs.W)
        fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Setup",special="footer")


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
    def __init__(self,dom,fs):
        super(PowerInflow, self).__init__(dom,fs)
        self.fprint("Setting Up Boundary Conditions",special="header")
        self.fprint("Type: Power Law Inflow")

        for key, values in self.boundary_types.items():
            self.fprint("Boundary Type: {0}, Applied to:".format(key))
            for value in values:
                self.fprint(value,offset=1)
        self.fprint("")
        ### Calculate the distance to the ground for the Q function space ###
        # self.z_dist_Q = Function(fs.Q)
        self.fprint("Computing Distance to Ground")
        self.height = Function(fs.Q)
        self.depth = Function(fs.Q)
        Q_coords = fs.Q.tabulate_dof_coordinates()
        height_vals = self.height.vector()[:]
        for i in range(len(Q_coords)):
            height_vals[i] = dom.Ground(Q_coords[i,0],Q_coords[i,1])
        self.z_dist_Q = Q_coords[:,2]-height_vals
        self.height.vector()[:]=height_vals
        self.depth.vector()[:]=self.z_dist_Q
        scaled_z_dist_Q_val = np.abs(np.divide(self.z_dist_Q,(dom.z_range[1]-dom.z_range[0])))

        ### Calculate the distance to the ground for the Q function space ###
        self.z_dist_V = Function(fs.V)
        V_coords = fs.V.tabulate_dof_coordinates()
        z_dist_V_val = np.zeros(len(V_coords))
        for i in range(len(V_coords)):
            z_dist_V_val[i] = V_coords[i,2]-dom.Ground(V_coords[i,0],V_coords[i,1])
        self.z_dist_V.vector()[:]=z_dist_V_val
        z_dist_V0,z_dist_V1,z_dist_V2 = self.z_dist_V.split(True)

        ### Create the Velocity Function ###
        self.fprint("Computing Velocity Vector")
        self.ux = Function(fs.V0)
        self.uy = Function(fs.V1)
        self.uz = Function(fs.V2)
        scaled_z_dist_val = np.abs(np.divide(z_dist_V0.vector()[:],(dom.z_range[1]-dom.z_range[0])))
        self.reference_velocity = np.multiply(8.0,np.power(scaled_z_dist_Q_val,.15))
        # ux.vector()[:] = np.multiply(16.0,np.power(scaled_z_dist_val,1./4.))

        ux_com, uy_com, uz_com = self.RotateVelocity(self.wind_direction)
        self.ux.vector()[:] = ux_com
        self.uy.vector()[:] = uy_com
        self.uz.vector()[:] = uz_com


        ### Assigning Velocity
        self.bc_velocity = Function(self.fs.V)
        self.fs.VelocityAssigner.assign(self.bc_velocity,[self.ux,self.uy,self.uz])

        ### Create Pressure Boundary Function
        self.bc_pressure = Function(self.fs.Q)

        ### Create Initial Guess
        self.fprint("Assigning Initial Guess")
        self.u0 = Function(fs.W)
        self.fs.SolutionAssigner.assign(self.u0,[self.bc_velocity,self.bc_pressure])

        ### Setup the boundary Conditions ###
        self.SetupBoundaries()
        self.fprint("Boundary Condition Setup",special="footer")

# class LogLayerInflow(object):
#     def __init__(self,dom,fs):