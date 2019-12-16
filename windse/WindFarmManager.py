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

    ### Import the cumulative parameters ###
    from windse import windse_parameters, BaseHeight

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

def TurbineForceNumpy(locs,yaws,axials,Ws,RDs,fs,inflow_angle=0.0,force_type="constant"):
    """

    """
    ### Set up some useful values ###
    dim = fs.V.mesh().topology().dim()
    locs = np.array(locs).T
    numturbs = locs.shape[0]

    ### Get the coordinates for a single component of the velocity function space ###
    coordinates = fs.V0.tabulate_dof_coordinates()
    x = coordinates[:,0]
    y = coordinates[:,1]
    if dim == 3:
        z = coordinates[:,2]

    ### Initialize Functions ###
    rotor_disks = Function(fs.V0)
    tf_x = Function(fs.V0)
    tf_y = Function(fs.V1)
    tf_z = Function(fs.V2)
    tf1  = Function(fs.V)
    tf2  = Function(fs.V)
    tf3  = Function(fs.V)
    tf_temp = Function(fs.V)

    ### For each turbine compute the space kernel and force ###
    for i in range(numturbs):

        ### Collect Turbine Specfic Data ###
        x0 = locs[i,0]
        y0 = locs[i,1]
        z0 = locs[i,2]
        yaw = yaws[i]+inflow_angle
        a = axials[i]
        W = Ws[i]
        R = RDs[i]/2.0
        A = pi*R**2.0

        ### Yaw and translate the turbine ###
        xrot =   math.cos(yaw)*(x-x0) + math.sin(yaw)*(y-y0)
        yrot = - math.sin(yaw)*(x-x0) + math.cos(yaw)*(y-y0)
        if dim == 3:
            zrot = z-z0
        else:
            zrot = 0.0

        ### Create the Thickness Kernel ###
        T_norm = 1.855438667500383
        T = np.exp(-np.power((xrot/W),6.0))/(T_norm*W)

        ### Create the Disk Kernel ###
        D_norm = 2.914516237206873
        r = np.power(yrot,2.0)+np.power(zrot,2.0)
        D = np.exp(-np.power(r/R**2.0,6.0))/(D_norm*R**2.0)

        ### Create Radial Force ###
        if force_type == "constant":
            C_t = 4./3.
            F = 0.5*A*C_t*a/(1.-a)
        elif force_type == "sine":
            F = 4.*0.5*A*a/(1.-a)*(r/R*np.sin(pi*r/R)+0.5) #* 1/(.81831)

        ### Assemble the turbine force components ##
        tf_x.vector()[:] = F*T*D*cos(yaw)
        tf_y.vector()[:] = F*T*D*sin(yaw)
        rotor_disks.vector()[:] += T*D 

        ### Create a tempory turbine force vector ###
        if dim == 3:
            fs.VelocityAssigner.assign(tf_temp,[tf_x,tf_y,tf_z])
        else:
            fs.VelocityAssigner.assign(tf_temp,[tf_x,tf_y])

        ### Create the 3 Turbine force function to completely reconstruct tf*(u.n)^2 without velocity ###
        tf1.vector()[:] += tf_temp.vector()[:] * cos(yaw)**2
        tf2.vector()[:] += tf_temp.vector()[:] * sin(yaw)**2
        tf3.vector()[:] += tf_temp.vector()[:] * cos(yaw) * sin(yaw)

    return (tf1,tf2,tf3,rotor_disks)

class GenericWindFarm(object):
    """
    A GenericProblem contains on the basic functions required by all problem objects.
    
    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, dom):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.force = self.params["wind_farm"].get("force","sine")
        self.analytic = self.params["domain"].get("analytic",False)
        self.dom = dom
        self.rd_first_save = True
        self.fprint = self.params.fprint
        self.extra_kwarg = {}
        if self.params["general"].get("dolfin_adjoint", False):
            self.extra_kwarg["annotate"] = False

    def Plot(self,show=True,filename="wind_farm"):
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
        ex_list_x = [self.ex_x[0],self.ex_x[1],self.ex_x[1],self.ex_x[0],self.ex_x[0]]
        ex_list_y = [self.ex_y[0],self.ex_y[0],self.ex_y[1],self.ex_y[1],self.ex_y[0]]

        ### Generate and Save Plot ###
        plt.figure()
        if hasattr(self.dom,"boundary_line"):
            plt.plot(*self.dom.boundary_line,c="k")
        plt.plot(ex_list_x,ex_list_y,c="r")
        p=plt.scatter(self.x,self.y,c=self.z,cmap="coolwarm",edgecolors=(0, 0, 0, 1))
        # p=plt.scatter(self.x,self.y,c="k",s=70)
        plt.xlim(self.dom.x_range[0],self.dom.x_range[1])
        plt.ylim(self.dom.y_range[0],self.dom.y_range[1])
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('Hub Elevation')

        plt.title("Location of the Turbines")
        plt.savefig(file_string, transparent=True)
        if show:
            plt.show()

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
        output = np.array([self.x, self.y, self.HH, self.yaw, self.RD, self.W, self.a])
        np.savetxt(file_string,output.T,header=head_str)


    def SaveActuatorDisks(self,val=0):
        """
        This function saves the turbine force if exists to output/.../functions/
        """
        if isinstance(self.actuator_disks,Function):
            if self.rd_first_save:
                self.rd_file = self.params.Save(self.actuator_disks,"actuator_disks",subfolder="functions/",val=val)
                self.rd_first_save = False
            else:
                self.params.Save(self.actuator_disks,"actuator_disks",subfolder="functions/",val=val,file=self.rd_file)

    def GetLocations(self):
        """
        This function gets three lists containing the x, y, and z locations
        of each turbine.

        Returns:
            x, y, z (lists): turbine coordinates
        """
        return self.x, self.y, self.z

    def PrintLocations(self):
        """
        This function prints out the  locations of each turbine.
        """
        for i in range(self.numturbs):
            print("Turbine "+repr(i+1)+": "+repr([self.x[i],self.y[i],self.z[i]]))

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

    def UpdateConstants(self,m=None,indexes=None,control_types=None):

        if m is not None:
            m = np.array(m)
            if "layout" in control_types:
                self.x = m[indexes[0]]
                self.y = m[indexes[1]]
            if "yaw" in control_types:
                self.yaw = m[indexes[2]]
            if "axial" in control_types:
                self.a = m[indexes[3]]

        for i in range(self.numturbs):
            self.mx[i].assign(self.x[i])
            self.my[i].assign(self.y[i])
            if self.analytic:
                self.mz[i] = self.dom.Ground(self.mx[i],self.my[i])+float(self.HH[i])
            else:
                self.mz[i] = BaseHeight(self.mx[i],self.my[i],self.dom.Ground)+float(self.HH[i])
            self.z[i] = float(self.mz[i])
            self.ground[i] = self.z[i] - self.HH[i]
            self.ma[i].assign(self.a[i])
            self.myaw[i].assign(self.yaw[i])

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
            if self.analytic:
                self.mz.append(self.dom.Ground(self.mx[i],self.my[i])+float(self.HH[i]))
            else:
                self.mz.append(BaseHeight(self.mx[i],self.my[i],self.dom.Ground)+float(self.HH[i]))
            self.z[i] = float(self.mz[i])
            self.ground[i] = self.z[i] - self.HH[i]

    def RotateFarm(self,angle):
        """
        This function rotates the position of each turbine. It does not change 
        the yaw of the turbines. 

        Args:
            angle (float): the rotation angle in radians
        """

        center = [sum(self.dom.x_range)/2.0,sum(self.dom.y_range)/2.0,sum(self.dom.z_range)/2.0]
        for i in range(self.numturbs):
            x = [self.x[i],self.y[i],self.z[i]]
            self.x[i] = math.cos(angle)*(x[0]-center[0]) - math.sin(angle)*(x[1]-center[1])+center[0]
            self.y[i] = math.sin(angle)*(x[0]-center[0]) + math.cos(angle)*(x[1]-center[1])+center[1]
            self.z[i] = self.HH[i]+self.dom.ground(self.x[i],self.y[i])[0]
            # self.angle[i] -= rot
        self.CalculateExtents()
        self.UpdateConstants()

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

    def CreateRotorDiscs(self,fs,mesh,delta_yaw=0.0):
        tf_start = time.time()
        self.fprint("Creating Rotor Discs",special="header")
        x=SpatialCoordinate(mesh)

        self.discs = Function(fs.Q)
        # self.discs = Function(fs.V)

        for i in range(self.numturbs):
            x0 = [self.mx[i],self.my[i],self.mz[i]]
            yaw = self.myaw[i]+delta_yaw
            W = self.W[i]/2.0
            R = self.RD[i]/2.0
            
            ### Calculate the normal vector ###
            n = Constant((cos(yaw),sin(yaw),0.0))

            ### Rotate and Shift the Turbine ###
            xs = self.YawTurbine(x,x0,yaw)

            ### Create the function that represents the Thickness of the turbine ###
            T_norm = 1.902701539733748
            T = exp(-pow((xs[0]/W),10.0))/(T_norm*W)

            ### Create the function that represents the Disk of the turbine
            D_norm = 2.884512175878827
            D = exp(-pow((pow((xs[1]/R),2)+pow((xs[2]/R),2)),5.0))/(D_norm*R**2.0)

            ### Combine and add to the total ###
            self.discs = self.discs + T*D
            # self.discs = self.discs + T*D*n

        self.fprint("Projecting Rotor Discs")
        self.discs = project(self.discs,fs.Q,solver_type='mumps',**self.extra_kwarg)

        tf_stop = time.time()
        self.fprint("Rotor Discs Created: {:1.2f} s".format(tf_stop-tf_start),special="footer")



    def YawTurbine(self,x,x0,yaw):
        """
        This function yaws the turbines when creating the turbine force.

        Args:
            x (dolfin.SpacialCoordinate): the space variable, x
            x0 (list): the location of the turbine to be yawed
            yaw (float): the yaw value in radians
        """
        xrot =   math.cos(yaw)*(x[0]-x0[0]) + math.sin(yaw)*(x[1]-x0[1])
        yrot = - math.sin(yaw)*(x[0]-x0[0]) + math.cos(yaw)*(x[1]-x0[1])
        if self.dom.dim == 3:
            zrot = x[2]-x0[2]
        else:
            zrot = 0.0
        
        return [xrot,yrot,zrot]

    def TurbineForce_numpy(self,fs,mesh,u_next,delta_yaw=0.0):
        tf_start = time.time()
        self.fprint("Calculating Turbine Force",special="header")

        tf1, tf2, tf3, rotor_disks = TurbineForceNumpy([self.x,self.y,self.z],self.yaw,self.a,self.W,self.RD,fs,inflow_angle=delta_yaw,force_type=self.force)
        self.actuator_disks = rotor_disks

        tf_stop = time.time()
        self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
        return (tf1, tf2, tf3)

    def TurbineForce(self,fs,mesh,u_next,delta_yaw=0.0):
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

        x=SpatialCoordinate(mesh)
        tf=0
        tf1=0
        tf2=0
        tf3=0
        for i in range(self.numturbs):
            x0 = [self.mx[i],self.my[i],self.mz[i]]
            yaw = self.myaw[i]+delta_yaw
            W = self.W[i]*1.0
            R = self.RD[i]/2.0
            A = pi*R**2.0 
            ma = self.ma[i]
            if self.dom.dim == 3:
                WTGbase = Expression(("cos(yaw)","sin(yaw)","0.0"),yaw=float(yaw),degree=1)
            else:
                WTGbase = Expression(("cos(yaw)","sin(yaw)"),yaw=float(yaw),degree=1)

            ### Rotate and Shift the Turbine ###
            xs = self.YawTurbine(x,x0,yaw)

            ### Create the function that represents the Thickness of the turbine ###
            T_norm = 1.855438667500383
            T = exp(-pow((xs[0]/W),6.0))/(T_norm*W)

            ### Create the function that represents the Disk of the turbine
            D_norm = 2.914516237206873
            D = exp(-pow((pow((xs[1]/R),2)+pow((xs[2]/R),2)),6.0))/(D_norm*R**2.0)

            ### Create the function that represents the force ###
            if self.force == "constant":
                C_t = 4./3.
                F = 0.5*A*C_t*ma/(1.-ma)
            elif self.force == "sine":
                r = sqrt(xs[1]**2.0+xs[2]**2)
                F = 4.*0.5*A*ma/(1.-ma)*(r/R*sin(pi*r/R)+0.5)/(.81831)

            # compute disk averaged velocity in yawed case and don't project
            u_d = u_next[0]*cos(yaw) + u_next[1]*sin(yaw)
            tf  += F*T*D*WTGbase*u_d**2
            # tf1 += F*T*D*WTGbase * cos(yaw)**2#*u_d**2
            # tf2 += F*T*D*WTGbase * sin(yaw)**2#*u_d**2
            # tf3 += F*T*D*WTGbase * cos(yaw) * sin(yaw)#*u_d**2

        ### Project Turbine Force to save on Assemble time ###
        self.fprint("Projecting Turbine Force")
        # self.actuator_disks = None
        self.actuator_disks = project(tf,fs.V,solver_type='mumps',**self.extra_kwarg)

        tf_stop = time.time()
        self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
        # return (tf1, tf2, tf3)
        return tf

    # def TurbineForce2D(self,fs,mesh):
    #     """
    #     This function creates a turbine force by applying 
    #     a spatial kernel to each turbine. This kernel is 
    #     created from the turbines location, yaw, thickness, diameter,
    #     and force density. Currently, force density is limit to a scaled
    #     version of 

    #     .. math::

    #         r\\sin(r),

    #     where :math:`r` is the distance from the center of the turbine.

    #     Args:
    #         V (dolfin.FunctionSpace): The function space the turbine force will use.
    #         mesh (dolfin.mesh): The mesh

    #     Returns:
    #         tf (dolfin.Function): the turbine force.

    #     Todo:
    #         * Setup a way to get the force density from file
    #     """

    #     tf_start = time.time()
    #     self.fprint("Calculating Turbine Force",special="header")
    #     x=SpatialCoordinate(mesh)

    #     tf_x=Function(fs.V0)
    #     tf_y=Function(fs.V1)

    #     for i in range(self.numturbs):
    #         x0 = [self.mx[i],self.my[i]]
    #         yaw = self.myaw[i]
    #         W = self.W[i]/2.0
    #         R = self.RD[i]/2.0 
    #         ma = self.ma[i]

    #         ### Rotate and Shift the Turbine ###
    #         xs = self.YawTurbine2D(x,x0,yaw)

    #         ### Create the function that represents the Thickness of the turbine ###
    #         T_norm = 1.902701539733748
    #         T = exp(-pow((xs[0]/W),10.0))/(T_norm*W)

    #         ### Create the function that represents the Disk of the turbine
    #         D_norm = 2.884512175878827
    #         D = exp(-pow((pow((xs[1]/R),2)),5.0))/(D_norm*R**2.0)

    #         ### Create the function that represents the force ###
    #         # F = 0.75*0.5*4.*A*self.ma[i]/(1.-self.ma[i])/beta
    #         r = sqrt(xs[1]**2.0)
    #         # F = 4.*0.5*(pi*R**2.0)*ma/(1.-ma)*(r/R*sin(pi*r/R)+0.5)
    #         F = 4.*0.5*(pi*R**2.0)*ma/(1.-ma)*(r/R*sin(pi*r/R)+0.5)

    #         ### Combine and add to the total ###
    #         tf_x = tf_x + F*T*D*cos(yaw)
    #         tf_y = tf_y + F*T*D*sin(yaw)

    #     ### Project Turbine Force to save on Assemble time ###
    #     self.fprint("Projecting X Force")
    #     tf_x = project(tf_x,fs.V0,solver_type='mumps')
    #     self.fprint("Projecting Y Force")
    #     tf_y = project(tf_y,fs.V1,solver_type='mumps')    

    #     ### Assign the components to the turbine force ###
    #     self.tf = Function(fs.V)
    #     fs.VelocityAssigner.assign(self.tf,[tf_x,tf_y])

    #     tf_stop = time.time()
    #     self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
    #     return as_vector((tf_x,tf_y))

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

        self.fprint("Generating Grid Wind Farm",special="header")

        ### Initialize Values from Options ###
        self.grid_rows = self.params["wind_farm"]["grid_rows"]
        self.grid_cols = self.params["wind_farm"]["grid_cols"]
        self.numturbs = self.grid_rows * self.grid_cols
        self.params["wind_farm"]["numturbs"] = self.numturbs

        self.HH = [self.params["wind_farm"]["HH"]]*self.numturbs
        self.RD = [self.params["wind_farm"]["RD"]]*self.numturbs
        self.W = [self.params["wind_farm"]["thickness"]]*self.numturbs
        self.yaw = [self.params["wind_farm"]["yaw"]]*self.numturbs
        self.axial = [self.params["wind_farm"]["axial"]]*self.numturbs
        self.radius = self.RD[0]/2.0
        self.jitter = self.params["wind_farm"].get("jitter",0.0)
        self.seed = self.params["wind_farm"].get("seed",None)

        self.ex_x = self.params["wind_farm"]["ex_x"]
        self.ex_y = self.params["wind_farm"]["ex_y"]

        ### Print some useful stats ###
        self.fprint("Force Type:         {0}".format(self.force))
        self.fprint("Number of Rows:     {:d}".format(self.grid_rows))
        self.fprint("Number of Columns:  {:d}".format(self.grid_cols))
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))
        if self.jitter > 0.0:
            self.fprint("Amount of Jitter:   {: 1.2f}".format(self.jitter))
            self.fprint("Random Seed: " + repr(self.seed))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0],self.ex_x[1]))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0],self.ex_y[1]))

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
        self.fprint("Generating Random Farm",special="header")

        ### Initialize Values from Options ###
        self.numturbs = self.params["wind_farm"]["numturbs"]
        
        self.HH = [self.params["wind_farm"]["HH"]]*self.numturbs
        self.RD = [self.params["wind_farm"]["RD"]]*self.numturbs
        self.W = [self.params["wind_farm"]["thickness"]]*self.numturbs
        self.yaw = [self.params["wind_farm"]["yaw"]]*self.numturbs
        self.axial = [self.params["wind_farm"]["axial"]]*self.numturbs
        self.radius = self.RD[0]/2.0

        self.ex_x = self.params["wind_farm"]["ex_x"]
        self.ex_y = self.params["wind_farm"]["ex_y"]

        self.seed = self.params["wind_farm"].get("seed",None)
        

        ### Print some useful stats ###
        self.fprint("Force Type:         {0}".format(self.force))
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0],self.ex_x[1]))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0],self.ex_y[1]))
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
        self.fprint("Importing Wind Farm",special="header")
        
        ### Import the data from path ###
        self.path = self.params["wind_farm"]["path"]
        raw_data = np.loadtxt(self.path,comments="#")

        ### Copy Files to input folder ###
        shutil.copy(self.path,self.params.folder+"input_files/")

        ### Parse the data ###
        if len(raw_data.shape) > 1:
            self.x     = raw_data[:,0] 
            self.y     = raw_data[:,1]
            self.HH    = raw_data[:,2]
            self.yaw   = raw_data[:,3]
            self.RD    = raw_data[:,4]
            self.radius = self.RD/2.0
            self.W     = raw_data[:,5]
            self.a     = raw_data[:,6]
            self.numturbs = len(self.x)

        else:
            self.x     = np.array((raw_data[0],))
            self.y     = np.array((raw_data[1],))
            self.HH    = np.array((raw_data[2],))
            self.yaw   = np.array((raw_data[3],))
            self.RD    = np.array((raw_data[4],))
            self.radius = np.array((raw_data[4]/2.0,))
            self.W     = np.array((raw_data[5],))
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
