"""
The windfarm manager contains everything required to set up a 
windfarm.
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    from sys import platform
    import math
    import time
    import shutil
    from scipy.special import gamma

    ### Import the cumulative parameters ###
    from windse import windse_parameters, BaseHeight, CalculateDiskTurbineForces

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

class GenericWindFarm(object):
    """
    A GenericProblem contains on the basic functions required by all problem objects.
    
    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, dom):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.turbine_method = self.params["wind_farm"].get("turbine_method","numpy")
        self.force = self.params["wind_farm"].get("force","sine")

        self.analytic = self.params["domain"].get("analytic",False)
        self.dom = dom
        self.rd_first_save = True
        self.fprint = self.params.fprint
        self.extra_kwarg = {}
        if self.params["general"].get("dolfin_adjoint", False):
            self.control_types = self.params["optimization"]["controls"]
            self.extra_kwarg["annotate"] = False

    def Plot(self,show=False,filename="wind_farm",power=None):
        """
        This function plots the locations of each wind turbine and
        saves the output to output/.../plots/

        :Keyword Arguments:
            * **show** (*bool*): Default: True, Set False to suppress output but still save.
        """
        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/"+filename+".pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        ### Create a list that outlines the extent of the farm ###
        if self.params.get("optimization",{}):
            Sx = self.dom.xscale
            layout_bounds = np.array(self.params["optimization"].get("layout_bounds",[np.array(self.ex_x)/Sx,np.array(self.ex_y)/Sx]))*Sx
            ex_list_x = [layout_bounds[0][0],layout_bounds[0][1],layout_bounds[0][1],layout_bounds[0][0],layout_bounds[0][0]]
            ex_list_y = [layout_bounds[1][0],layout_bounds[1][0],layout_bounds[1][1],layout_bounds[1][1],layout_bounds[1][0]]
        else:
            ex_list_x = [self.ex_x[0],self.ex_x[1],self.ex_x[1],self.ex_x[0],self.ex_x[0]]
            ex_list_y = [self.ex_y[0],self.ex_y[0],self.ex_y[1],self.ex_y[1],self.ex_y[0]]

        ### Generate and Save Plot ###
        plt.figure()
        if hasattr(self.dom,"boundary_line"):
            plt.plot(*self.dom.boundary_line/self.dom.xscale,c="k")
        plt.plot(np.array(ex_list_x)/self.dom.xscale, np.array(ex_list_y)/self.dom.xscale,c="r")

        ### Plot Blades
        for i in range(self.numturbs):
            blade_n = [np.cos(self.yaw[i]),np.sin(self.yaw[i])]
            rr = self.RD[i]/2.0
            blade_x = np.array([self.x[i]+rr*blade_n[1],self.x[i]-rr*blade_n[1]])/self.dom.xscale
            blade_y = np.array([self.y[i]-rr*blade_n[0],self.y[i]+rr*blade_n[0]])/self.dom.xscale
            plt.plot(blade_x,blade_y,c='k',linewidth=2,zorder=1)

        ### Choose coloring for the turbines ###
        if isinstance(power,(list,np.ndarray)):
            coloring = power
        else:
            coloring = np.array(self.z)/self.dom.xscale

        ### Plot Hub Locations
        p=plt.scatter(self.x/self.dom.xscale,self.y/self.dom.xscale,c=coloring,cmap="coolwarm",edgecolors=(0, 0, 0, 1),s=20,zorder=2)
        # p=plt.scatter(self.x,self.y,c="k",s=70)
        plt.xlim(self.dom.x_range[0]/self.dom.xscale,self.dom.x_range[1]/self.dom.xscale)
        plt.ylim(self.dom.y_range[0]/self.dom.xscale,self.dom.y_range[1]/self.dom.xscale)
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('Hub Elevation')

        if power is None:
            plt.title("Location of the Turbines")
        elif isinstance(power,(list,np.ndarray)):
            plt.title("Power Output: {: 5.2f}".format(sum(power)))
        else:
            plt.title("Power Output: {: 5.2f}".format(power))

        plt.savefig(file_string, transparent=True)

        if show:
            plt.show()

        plt.close()

    def SaveWindFarm(self,val=None,filename="wind_farm"):

        ### Create the path names ###
        folder_string = self.params.folder+"/data/"
        if val is not None:
            file_string = self.params.folder+"/data/"+filename+"_"+repr(val)+".txt"
        else:
            file_string = self.params.folder+"/data/"+filename+".txt"

        ### Check if folder exists ###
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        ### Define the header string ###
        head_str="#    x    y    HH    Yaw    Diameter    Thickness    Axial_Induction"


        ### Save text file ###
        Sx = self.dom.xscale
        output = np.array([self.x/Sx, self.y/Sx, self.HH/Sx, self.yaw, self.RD/Sx, self.W/Sx, self.a])
        np.savetxt(file_string,output.T,header=head_str)

    def SaveActuatorDisks(self,val=0):
        """
        This function saves the turbine force if exists to output/.../functions/
        """

        self.dom.mesh.coordinates()[:]=self.dom.mesh.coordinates()[:]/self.dom.xscale
        if hasattr(self.actuator_disks,"_cpp_object"):
            if self.rd_first_save:
                self.rd_file = self.params.Save(self.actuator_disks,"actuator_disks",subfolder="functions/",val=val)
                self.rd_first_save = False
            else:
                self.params.Save(self.actuator_disks,"actuator_disks",subfolder="functions/",val=val,file=self.rd_file)

        self.dom.mesh.coordinates()[:]=self.dom.mesh.coordinates()[:]*self.dom.xscale

    def CalculateExtents(self):
        """
        This functions takes into consideration the turbine locations, diameters, 
        and hub heights to create lists that describe the extent of the windfarm.
        These lists are append to the parameters object.
        """
        ### Locate the extreme turbines ### 
        x_min = np.argmin(self.x)
        x_max = np.argmax(self.x)
        y_min = np.argmin(self.y)
        y_max = np.argmax(self.y)
        z_min = np.argmin(self.z)
        z_max = np.argmax(self.z)

        ### Calculate the extent of the farm ###
        self.ex_x = [self.x[x_min]-self.RD[x_min]/2.0,self.x[x_max]+self.RD[x_max]/2.0]
        self.ex_y = [self.y[y_min]-self.RD[y_min]/2.0,self.y[y_max]+self.RD[y_max]/2.0]
        self.ex_z = [min(self.ground),self.z[z_max]+self.RD[z_max]/2.0]

        ### Update the options ###
        self.params["wind_farm"]["ex_x"] = self.ex_x
        self.params["wind_farm"]["ex_y"] = self.ex_y
        self.params["wind_farm"]["ex_z"] = self.ex_z

    def CalculateFarmRegion(self,region_type,factor=1.0,length=None):


        if region_type == "custom":
            return length

        elif region_type == "square":
            x_center = (self.ex_x[1]+self.ex_x[0])/2.0
            y_center = (self.ex_y[1]+self.ex_y[0])/2.0
            z_center = (self.ex_z[1]+self.ex_z[0])/2.0

            if length is None:
                x = factor*(np.subtract(self.ex_x,x_center))+x_center
                y = factor*(np.subtract(self.ex_y,y_center))+y_center
                z = factor*(self.ex_z - z_center)+z_center
            else:
                x = [-length/2.0+x_center,length/2.0+x_center]
                y = [-length/2.0+x_center,length/2.0+x_center]
                # z = 1.2*(self.ex_z - z_center)+z_center
                z = [self.dom.z_range[0], np.mean(self.HH)+np.mean(self.RD)]

            return [x,y,z]

        elif "circle" in region_type:
            if region_type == "farm_circle":
                center = [sum(self.x)/float(self.numturbs),sum(self.y)/float(self.numturbs)]
            else:
                center = [sum(self.dom.x_range)/2.0,sum(self.dom.y_range)/2.0]

            z_center = (self.ex_z[1]+self.ex_z[0])/2.0

            if length is None:
                length = factor*max(np.sqrt(np.power(self.x-center[0],2.0)+np.power(self.y-center[1],2.0)))
            
            # z = 1.2*(self.ex_z - z_center)+z_center
            z = [self.dom.z_range[0], np.mean(self.HH)+np.mean(self.RD)]

            return [[length],center,z]

        else:
            raise ValueError("Not a valid region type: "+region_type)

    def CreateConstants(self):
        """
        This functions converts lists of locations and axial inductions
        into dolfin.Constants. This is useful in optimization.
        """
        self.mx = []
        self.my = []
        self.ma = []
        self.myaw = []
        for i in range(self.numturbs):
            self.mx.append(Constant(self.x[i]))
            self.my.append(Constant(self.y[i]))
            self.ma.append(Constant(self.a[i]))
            self.myaw.append(Constant(self.yaw[i]))

        for i in range(self.numturbs):
            self.mx[i].rename("x"+repr(i),"x"+repr(i))
            self.my[i].rename("y"+repr(i),"y"+repr(i))
            self.myaw[i].rename("yaw"+repr(i),"yaw"+repr(i))
            self.ma[i].rename("a"+repr(i),"a"+repr(i))

    def UpdateControls(self,x=None,y=None,yaw=None,a=None):

        if x is not None:
            self.x = np.array(x,dtype=float)
        if y is not None:
            self.y = np.array(y,dtype=float)
        if yaw is not None:
            self.yaw = np.array(yaw,dtype=float)
        if a is not None:
            self.a = np.array(a,dtype=float)

        for i in range(self.numturbs):
            self.mx[i]=Constant(self.x[i])
            self.my[i]=Constant(self.y[i])
            # if self.analytic:
            #     self.mz[i] = self.dom.Ground(self.mx[i],self.my[i])+float(self.HH[i])
            # else:
            self.mz[i] = BaseHeight(self.mx[i],self.my[i],self.dom.Ground)+float(self.HH[i])
            self.z[i] = float(self.mz[i])
            self.ground[i] = self.z[i] - self.HH[i]
            self.ma[i]=Constant(self.a[i])
            self.myaw[i]=Constant(self.yaw[i])

    def CreateLists(self):
        """
        This function creates lists from single values. This is useful
        when the params.yaml file defines only one type of turbine.
        """
        self.HH = np.full(self.numturbs,self.HH)
        self.RD = np.full(self.numturbs,self.RD)
        self.W = np.full(self.numturbs,self.W)
        self.radius = np.full(self.numturbs,self.radius)
        self.yaw = np.full(self.numturbs,self.yaw)
        self.a = np.full(self.numturbs,self.axial)

    def CalculateHeights(self):
        """
        This function calculates the absolute heights of each turbine.
        """
        self.mz = [] 
        self.z = np.zeros(self.numturbs)
        self.ground = np.zeros(self.numturbs)
        for i in range(self.numturbs):
            # if self.analytic:
            #     self.mz.append(self.dom.Ground(self.mx[i],self.my[i])+float(self.HH[i]))
            # else:
            self.mz.append(BaseHeight(self.mx[i],self.my[i],self.dom.Ground)+float(self.HH[i]))
            self.z[i] = float(self.mz[i])
            self.ground[i] = self.z[i] - self.HH[i]

    def RefineTurbines(self,num,radius_multiplyer):

        self.fprint("Refining Near Turbines",special="header")
        mark_start = time.time()

        for i in range(num):
            if num>1:
                step_start = time.time()
                self.fprint("Refining Mesh Step {:d} of {:d}".format(i+1,num), special="header")

            cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)


            expand_turbine_radius = False

            if expand_turbine_radius:
                radius = (num-i)*radius_multiplyer*np.array(self.RD)/2.0
            else:
                radius = radius_multiplyer*np.array(self.RD)/2.0


            if self.dom.dim == 3:
                turb_x = np.array(self.x)
                turb_x = np.tile(turb_x,(4,1)).T
                turb_y = np.array(self.y)
                turb_y = np.tile(turb_y,(4,1)).T
                turb_z0 = self.dom.z_range[0]-radius
                turb_z1 = self.z+radius
            else:
                turb_x = np.array(self.x)
                turb_x = np.tile(turb_x,(3,1)).T
                turb_y = np.array(self.y)
                turb_y = np.tile(turb_y,(3,1)).T
            n = self.numturbs

            self.fprint("Marking Near Turbine")
            for cell in cells(self.dom.mesh):

                pt = cell.get_vertex_coordinates()
                if self.dom.dim == 3:
                    x = pt[0:-2:3]
                    x = np.tile(x,(n,1))
                    y = pt[1:-1:3]
                    y = np.tile(y,(n,1))
                    z = pt[2::3]
                else:
                    x = pt[0:-1:2]
                    x = np.tile(x,(n,1))
                    y = pt[1::2]
                    y = np.tile(y,(n,1))

                ### For each turbine, find which vertex is closet using squared distance
                min_r = np.min(np.power(turb_x-x,2.0)+np.power(turb_y-y,2.0),axis=1)


                downstream_teardrop_shape = False

                if downstream_teardrop_shape:
                    min_arg = np.argmin(np.power(turb_x-x,2.0)+np.power(turb_y-y,2.0),axis=1)
                    min_arg = np.argmin(min_arg)

                    if np.any(turb_x[min_arg] < x):
                        min_r = np.min(np.power(turb_x-x/2.0,2.0)+np.power(turb_y-y,2.0),axis=1)
                    else:
                        min_r = np.min(np.power(turb_x-x*2.0,2.0)+np.power(turb_y-y,2.0),axis=1)


                in_circle = min_r<=radius**2.0
                if self.dom.dim == 3:
                    in_z = np.logical_and(turb_z0 <= max(z), turb_z1 >= min(z))
                    near_turbine = np.logical_and(in_circle, in_z)
                else:
                    near_turbine = in_circle

                if any(near_turbine):
                    cell_f[cell] = True
            mark_stop = time.time()
            self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

            self.dom.Refine(1,cell_markers=cell_f)

            if num>1:
                step_stop = time.time()
                self.fprint("Step {:d} of {:d} Finished: {:1.2f} s".format(i+1,num,step_stop-step_start), special="footer")

        self.CalculateHeights()
        self.fprint("Turbine Refinement Finished",special="footer")

    def YawTurbine(self,x,x0,yaw):
        """
        This function yaws the turbines when creating the turbine force.

        Args:
            x (dolfin.SpacialCoordinate): the space variable, x
            x0 (list): the location of the turbine to be yawed
            yaw (float): the yaw value in radians
        """
        xrot =   cos(yaw)*(x[0]-x0[0]) + sin(yaw)*(x[1]-x0[1])
        yrot = - sin(yaw)*(x[0]-x0[0]) + cos(yaw)*(x[1]-x0[1])
        if self.dom.dim == 3:
            zrot = x[2]-x0[2]
        else:
            zrot = 0.0
        
        return [xrot,yrot,zrot]

    def NumpyTurbineForce(self,fs,mesh,inflow_angle=0.0):
        tf_start = time.time()
        self.fprint("Calculating Turbine Force",special="header")
        self.fprint("Using a Numpy Representation")

        self.inflow_angle = inflow_angle
        x = fs.tf_V0.tabulate_dof_coordinates().T
        [tf1, tf2, tf3], sparse_ids, actuator_array = CalculateDiskTurbineForces(x, self, fs, save_actuators=True)

        self.fprint("Turbine Force Space:  {}".format(fs.tf_space))
        self.fprint("Turbine Force Degree: {:d}".format(fs.tf_degree))
        self.fprint("Quadrature DOFS:      {:d}".format(fs.tf_V.dim()))
        self.fprint("Turbine DOFs:         {:d}".format(len(sparse_ids)))
        self.fprint("Compression:          {:1.4f} %".format(len(sparse_ids)/fs.tf_V.dim()*100))

        ### Rename for Identification ###
        tf1.rename("tf1","tf1")
        tf2.rename("tf2","tf2")
        tf3.rename("tf3","tf3")

        ### Construct the actuator disks for post processing ###
        # self.actuator_disks_list = actuator_disks
        self.actuator_disks = Function(fs.tf_V)
        self.actuator_disks.vector()[:] = np.sum(actuator_array,axis=1)
        self.fprint("Projecting Turbine Force")
        self.actuator_disks = project(self.actuator_disks,fs.V,solver_type='mumps',form_compiler_parameters={'quadrature_degree': fs.tf_degree},**self.extra_kwarg)
        
        self.actuator_disks_list = []
        for i in range(self.numturbs):
            temp = Function(fs.tf_V)
            temp.vector()[:] = np.array(actuator_array[:,i])
            self.actuator_disks_list.append(temp)

        tf_stop = time.time()
        self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
        return (tf1, tf2, tf3)

    def DolfinTurbineForce(self,fs,mesh,inflow_angle=0.0):
        """
        This function creates a turbine force by applying 
        a spacial kernel to each turbine. This kernel is 
        created from the turbines location, yaw, thickness, diameter,
        and force density. Currently, force density is limit to a scaled
        version of 

        .. math::

            r\\sin(r),

        where :math:`r` is the distance from the center of the turbine.

        Args:
            V (dolfin.FunctionSpace): The function space the turbine force will use.
            mesh (dolfin.mesh): The mesh

        Returns:
            tf (dolfin.Function): the turbine force.

        Todo:
            * Setup a way to get the force density from file
        """
        tf_start = time.time()
        self.fprint("Calculating Turbine Force",special="header")
        self.fprint("Using a Dolfin Representation")

        x=SpatialCoordinate(mesh)
        tf=0
        rd=0
        tf1=0
        tf2=0
        tf3=0
        self.actuator_disks_list = []
        for i in range(self.numturbs):
            x0 = [self.mx[i],self.my[i],self.mz[i]]
            yaw = self.myaw[i]+inflow_angle
            W = self.W[i]*1.0
            R = self.RD[i]/2.0
            ma = self.ma[i]
            C_tprime = 4*ma/(1-ma)

            ### Set up some dim dependent values ###
            S_norm = (2.0+pi)/(2.0*pi)
            T_norm = 2.0*gamma(7.0/6.0)
            if self.dom.dim == 3:
                WTGbase = as_vector((cos(yaw),sin(yaw),0.0))
                A = pi*R**2.0 
                D_norm = pi*gamma(4.0/3.0)
            else:
                WTGbase = as_vector((cos(yaw),sin(yaw)))
                A = 2*R 
                D_norm = 2.0*gamma(7.0/6.0)

            ### Rotate and Shift the Turbine ###
            xs = self.YawTurbine(x,x0,yaw)

            ### Create the function that represents the Thickness of the turbine ###
            T = exp(-pow((xs[0]/W),6.0))

            ### Create the function that represents the Disk of the turbine
            r = sqrt(xs[1]**2.0+xs[2]**2.0)/R
            D = exp(-pow(r,6.0))

            ### Create the function that represents the force ###
            if self.force == "constant":
                F = 0.5*A*C_tprime
            elif self.force == "sine":
                F = 0.5*A*C_tprime*(r*sin(pi*r)+0.5)/S_norm

            ### Calculate normalization constant ###
            volNormalization = T_norm*D_norm*W*R**(self.dom.dim-1)
            # volNormalization = assemble(T*D*dx)
            # print(volNormalization_a,volNormalization,volNormalization/(W*R**(self.dom.dim-1)),T_norm*D_norm)

            # compute disk averaged velocity in yawed case and don't project
            self.actuator_disks_list.append(F*T*D*WTGbase/volNormalization)
            rd  += F*T*D*WTGbase/volNormalization
            tf1 += F*T*D*WTGbase/volNormalization * cos(yaw)**2
            tf2 += F*T*D*WTGbase/volNormalization * sin(yaw)**2
            tf3 += F*T*D*WTGbase/volNormalization * 2.0 * cos(yaw) * sin(yaw)
        # exit()

        ### Save the actuator disks for post processing ###
        self.fprint("Projecting Turbine Force")
        self.actuator_disks = project(rd,fs.V,solver_type='mumps',**self.extra_kwarg)

        tf_stop = time.time()
        self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
        return (tf1, tf2, tf3)

class GridWindFarm(GenericWindFarm):
    """
    A GridWindFarm produces turbines on a grid. The params.yaml file determines
    how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                #                     # Description              | Units
                HH: 90                # Hub Height               | m
                RD: 126.0             # Turbine Diameter         | m
                thickness: 10.5       # Effective Thickness      | m
                yaw: 0.0              # Yaw                      | rads
                axial: 0.33           # Axial Induction          | -
                ex_x: [-1500, 1500]   # x-extent of the farm     | m
                ex_y: [-1500, 1500]   # y-extent of the farm     | m
                grid_rows: 6          # Number of rows           | -
                grid_cols: 6          # Number of columns        | -

        This will produce a 6x6 grid of turbines equally spaced within the 
        region [-1500, 1500]x[-1500, 1500].

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):
        super(GridWindFarm, self).__init__(dom)
        Sx = self.dom.xscale

        self.fprint("Generating Grid Wind Farm",special="header")

        ### Initialize Values from Options ###
        self.grid_rows = self.params["wind_farm"]["grid_rows"]
        self.grid_cols = self.params["wind_farm"]["grid_cols"]
        self.numturbs = self.grid_rows * self.grid_cols
        self.params["wind_farm"]["numturbs"] = self.numturbs


        self.HH = [self.params["wind_farm"]["HH"]*Sx]*self.numturbs
        self.RD = [self.params["wind_farm"]["RD"]*Sx]*self.numturbs
        self.W = [self.params["wind_farm"]["thickness"]*Sx]*self.numturbs
        self.yaw = [self.params["wind_farm"]["yaw"]]*self.numturbs
        self.axial = [self.params["wind_farm"]["axial"]]*self.numturbs
        self.radius = self.RD[0]/2.0
        self.jitter = self.params["wind_farm"].get("jitter",0.0)*Sx
        self.seed = self.params["wind_farm"].get("seed",None)

        self.ex_x = np.array(self.params["wind_farm"]["ex_x"])*Sx
        self.ex_y = np.array(self.params["wind_farm"]["ex_y"])*Sx

        ### Print some useful stats ###
        self.fprint("Force Type:         {0}".format(self.force))
        self.fprint("Number of Rows:     {:d}".format(self.grid_rows))
        self.fprint("Number of Columns:  {:d}".format(self.grid_cols))
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))
        if self.jitter > 0.0:
            self.fprint("Amount of Jitter:   {: 1.2f}".format(self.jitter))
            self.fprint("Random Seed: " + repr(self.seed))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0]/Sx,self.ex_x[1]/Sx))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0]/Sx,self.ex_y[1]/Sx))

        ### Create the x and y coords ###
        self.grid_x = np.linspace(self.ex_x[0]+self.radius,self.ex_x[1]-self.radius,self.grid_cols)
        self.grid_y = np.linspace(self.ex_y[0]+self.radius,self.ex_y[1]-self.radius,self.grid_rows)

        ### Use the x and y coords to make a mesh grid ###
        self.x, self.y = np.meshgrid(self.grid_x,self.grid_y)
        self.x = self.x.flatten()
        self.y = self.y.flatten()

        ### Apply Jitter ###
        if self.jitter > 0.0:
            if self.seed is not None:
                np.random.seed(self.seed)
            self.x += np.random.randn(self.numturbs)*self.jitter
            self.y += np.random.randn(self.numturbs)*self.jitter

        ### Convert the constant parameters to lists ###
        self.CreateLists()

        ### Convert the lists into lists of dolfin Constants ###
        self.CreateConstants() 

        ### Calculate Ground Heights ###
        self.CalculateHeights()

        ### Update the extent in the z direction ###
        self.ex_z = [min(self.ground),max(self.z+self.RD)]
        self.params["wind_farm"]["ex_z"] = self.ex_z

        self.fprint("Wind Farm Generated",special="footer")

class RandomWindFarm(GenericWindFarm):
    """
    A RandomWindFarm produces turbines located randomly with a defined 
    range. The params.yaml file determines how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                #                     # Description              | Units
                HH: 90                # Hub Height               | m
                RD: 126.0             # Turbine Diameter         | m
                thickness: 10.5       # Effective Thickness      | m
                yaw: 0.0              # Yaw                      | rads
                axial: 0.33           # Axial Induction          | -
                ex_x: [-1500, 1500]   # x-extent of the farm     | m
                ex_y: [-1500, 1500]   # y-extent of the farm     | m
                numturbs: 36          # Number of Turbines       | -
                seed: 15              # Random Seed for Numpy    | -

        This will produce a 36 turbines randomly located within the 
        region [-1500, 1500]x[-1500, 1500]. The seed is optional but 
        useful for reproducing test.

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):
        super(RandomWindFarm, self).__init__(dom)
        Sx = self.dom.xscale
        self.fprint("Generating Random Farm",special="header")

        ### Initialize Values from Options ###
        self.numturbs = self.params["wind_farm"]["numturbs"]
        
        self.HH = [self.params["wind_farm"]["HH"]*Sx]*self.numturbs
        self.RD = [self.params["wind_farm"]["RD"]*Sx]*self.numturbs
        self.W = [self.params["wind_farm"]["thickness"]*Sx]*self.numturbs
        self.yaw = [self.params["wind_farm"]["yaw"]]*self.numturbs
        self.axial = [self.params["wind_farm"]["axial"]]*self.numturbs
        self.radius = self.RD[0]/2.0

        self.ex_x = np.array(self.params["wind_farm"]["ex_x"])*Sx
        self.ex_y = np.array(self.params["wind_farm"]["ex_y"])*Sx

        self.seed = self.params["wind_farm"].get("seed",None)
        

        ### Print some useful stats ###
        self.fprint("Force Type:         {0}".format(self.force))
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0]/Sx,self.ex_x[1]/Sx))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0]/Sx,self.ex_y[1]/Sx))
        self.fprint("Random Seed: " + repr(self.seed))

        ### Check if random seed is set ###
        if self.seed is not None:
            np.random.seed(self.seed)

        ### Create the x and y coords ###
        self.x = np.random.uniform(self.ex_x[0]+self.radius,self.ex_x[1]-self.radius,self.numturbs)
        self.y = np.random.uniform(self.ex_y[0]+self.radius,self.ex_y[1]-self.radius,self.numturbs)


        ### Convert the constant parameters to lists ###
        self.CreateLists()
        
        ### Convert the lists into lists of dolfin Constants ###
        self.CreateConstants() 

        ### Calculate Ground Heights ###
        self.CalculateHeights()

        ### Update the extent in the z direction ###
        self.ex_z = [min(self.ground),max(self.z+self.RD)]
        self.params["wind_farm"]["ex_z"] = self.ex_z


        self.fprint("Wind Farm Generated",special="footer")


class ImportedWindFarm(GenericWindFarm):
    """
    A ImportedWindFarm produces turbines located based on a text file.
    The params.yaml file determines how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                imported: true
                path: "inputs/wind_farm.txt"

        The "wind_farm.txt" needs to be set up like this::

            #    x      y     HH           Yaw Diameter Thickness Axial_Induction
            200.00 0.0000 80.000  0.0000000000      126      10.5            0.33
            800.00 0.0000 80.000  0.0000000000      126      10.5            0.33

        The first row isn't necessary. Each row defines a different turbine.

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):
        super(ImportedWindFarm, self).__init__(dom)
        Sx = self.dom.xscale
        self.fprint("Importing Wind Farm",special="header")
        
        ### Import the data from path ###
        self.path = self.params["wind_farm"]["path"]
        raw_data = np.loadtxt(self.path,comments="#")

        ### Copy Files to input folder ###
        shutil.copy(self.path,self.params.folder+"input_files/")

        ### Parse the data ###
        if len(raw_data.shape) > 1:
            self.x     = raw_data[:,0]*Sx 
            self.y     = raw_data[:,1]*Sx
            self.HH    = raw_data[:,2]*Sx
            self.yaw   = raw_data[:,3]
            self.RD    = raw_data[:,4]*Sx
            self.radius = self.RD/2.0
            self.W     = raw_data[:,5]*Sx
            self.a     = raw_data[:,6]
            self.numturbs = len(self.x)

        else:
            self.x     = np.array((raw_data[0],))*Sx
            self.y     = np.array((raw_data[1],))*Sx
            self.HH    = np.array((raw_data[2],))*Sx
            self.yaw   = np.array((raw_data[3],))
            self.RD    = np.array((raw_data[4],))*Sx
            self.radius = np.array((raw_data[4]/2.0,))*Sx
            self.W     = np.array((raw_data[5],))*Sx
            self.a     = np.array((raw_data[6],))
            self.numturbs = 1

        ### Update the options ###
        self.params["wind_farm"]["numturbs"] = self.numturbs
        self.fprint("Force Type:         {0}".format(self.force))
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))

        ### Convert the lists into lists of dolfin Constants ###
        self.CreateConstants() 

        ### Calculate Ground Heights ###
        self.CalculateHeights()

        ### Calculate the extent of the farm ###
        self.CalculateExtents()
    

        self.fprint("Wind Farm Imported",special="footer")
