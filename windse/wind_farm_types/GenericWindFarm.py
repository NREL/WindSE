from windse import windse_parameters
from windse.turbine_types import turbine_dict
import numpy as np
import time, os
from . import MeshFunction, CompiledSubDomain, Measure, cells, project, inner, FiniteElement, FunctionSpace, MixedElement, assemble, dx, parameters, Form, File
import matplotlib.pyplot as plt
from pyadjoint.tape import stop_annotating 
from pyadjoint import AdjFloat

class GenericWindFarm(object):
    """
    A GenericWindFarm contains on the basic functions and attributes required by all wind farm objects.
    
    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, dom):
        """
        Store anything needed prior to setup
        """

        # Store WindSE Objects
        self.dom = dom
        self.params = windse_parameters
        self.fprint = self.params.fprint
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        # Load any global parameter
        self.farm_type = self.params["wind_farm"]["type"]
        self.local_dx_scaling = self.params["wind_farm"]["local_dx_scaling"]
        self.use_local_tf_dx = self.params["wind_farm"]["use_local_tf_dx"]
        self.turbine_type = self.params["turbines"]["type"]

        # Set any global flags
        self.func_first_save = True
        self.power_first_save = True

        # Init blank turbine list
        self.turbines = []

        # Get Optimization Data for plotting
        self.optimizing = False
        if self.params.performing_opt_calc:
            self.layout_bounds = self.params["optimization"]["layout_bounds"]
            self.control_types = self.params["optimization"]["control_types"]
            self.optimizing = True

        ### extra_kwargs will reduce overhead on operations that we don't want dolfin_adjoint to track 
        self.extra_kwarg = {}
        if self.params.dolfin_adjoint:
            self.extra_kwarg["annotate"] = False

        ### Set the number of turbine when there are just too many and we need to project before solving
        self.too_many_turbines = 60

        # Setup the wind farm (eventually this will happen outside the init)
        self.fprint(f"Generating {self.name}",special="header")
        self.setup()
        self.fprint(f"{self.name} Generated",special="footer")

    def setup(self):
        """
        This function builds the wind farm as well as sets up the turbines
        """  
        self.load_parameters()      
        self.compute_parameters()      
        self.fprint("Number of Turbines: {:d}".format(self.numturbs))
        self.fprint("Type of Turbines: {}".format(self.turbine_type))
        self.initial_turbine_locations = self.initialize_turbine_locations()

        self.fprint("Setting Up Turbines: ",special="header")
        self.setup_turbines()
        self.fprint("Turbines Set up",special="footer")

        self.debug_output() 

    def load_parameters(self):
        """
        This function will parse the parameters from the yaml file
        """  
        raise NotImplementedError(type(self))

    def compute_parameters(self):
        """
        This function will compute any additional parameters
        """  
        raise NotImplementedError(type(self))

    def initialize_turbine_locations():
        """
        This function will compute the initial locations of the turbines in the farm
        It must return an Nx2 array where N is the number of turbines
        """
        raise NotImplementedError(type(self))

    def setup_turbines(self):
        """
        Using the parameters and initial locations this function will populate the list of turbines
        """
        turbine_method = turbine_dict[self.turbine_type]
        for i,(x,y) in enumerate(self.initial_turbine_locations):
            self.turbines.append(turbine_method(i,x,y,self.dom))

    def get_hub_locations(self):
        """
        returns a nx3 numpy array containing the x, y, and z location of the turbine hubs
        """
        temp = np.zeros((self.numturbs,3))
        for i,turb in enumerate(self.turbines):
            temp[i] = [turb.mx,turb.my,turb.mz]
        return temp

    def get_ground_heights(self):
        """
        returns a nx1 numpy array containing the hight of the ground directly below the turbine hub
        """
        temp = np.zeros(self.numturbs)
        for i,turb in enumerate(self.turbines):
            temp[i] = turb.ground
        return temp

    def get_rotor_diameters(self):
        """
        returns a nx1 numpy array containing the rotor diameter of the turbine hubs
        """
        temp = np.zeros(self.numturbs)
        for i,turb in enumerate(self.turbines):
            temp[i] = turb.RD
        return temp

    def get_yaw_angles(self):
        """
        returns a nx1 numpy array containing the rotor diameter of the turbine hubs
        """
        temp = np.zeros(self.numturbs)
        for i,turb in enumerate(self.turbines):
            temp[i] = turb.myaw
        return temp

    def calculate_farm_bounding_box(self):
        """
        This functions takes into consideration the turbine locations, diameters, 
        and hub heights to create lists that describe the extent of the windfarm.
        These lists are append to the parameters object.
        """
        x,y,z = self.get_hub_locations().T
        RD = self.get_rotor_diameters()
        ground = self.get_ground_heights()

        ### Locate the extreme turbines ### 
        x_min = np.argmin(x)
        x_max = np.argmax(x)
        y_min = np.argmin(y)
        y_max = np.argmax(y)
        z_min = np.argmin(z)
        z_max = np.argmax(z)

        ### Calculate the extent of the farm ###
        self.ex_x = [x[x_min]-RD[x_min]/2.0,x[x_max]+RD[x_max]/2.0]
        self.ex_y = [y[y_min]-RD[y_min]/2.0,y[y_max]+RD[y_max]/2.0]
        self.ex_z = [min(ground),z[z_max]+RD[z_max]/2.0]
        
        return [self.ex_x,self.ex_y,self.ex_z]

    def update_heights(self):
        """
        updates the hub height for all turbines in farm 
        """
        for turb in self.turbines:
            turb.calculate_heights()

    def compute_turbine_force(self,u,v,inflow_angle,fs,**kwargs):
        """
        Iterates over the turbines and adds up the turbine forces
        """
        # store the function space for later use
        self.fs = fs

        # create a list of turbine forces and integration domain
        self.tf_list = []

        # Create a cell function that will be used to store the subdomains
        self.turbine_subdomains = MeshFunction("size_t", self.dom.mesh, self.dom.mesh.topology().dim())
        self.turbine_subdomains.set_all(0)

        # Iterate over each turbine to build its individual subdomain
        for turb in self.turbines:
            self.tf_list.append(turb.turbine_force(u,inflow_angle,fs,**kwargs))

            # the subdomain is a sphere centered at the hub of the turbine with a height of scaling*RD and diameter of scaling*RD

            if self.dom.dim == 3:
                subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)+(x[2]-HH)*(x[2]-HH)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, HH=turb.z, x0=turb.x, y0=turb.y)
            else:
                subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, x0=turb.x, y0=turb.y)


            # add the subdomain to the marker
            subdomain_marker.mark(self.turbine_subdomains,turb.index+1)


        # create a measure for the turbines
        self.local_dx = Measure('dx', subdomain_data=self.turbine_subdomains)

        # check if the number of turbines might trigger a python recursion error.    
        if self.numturbs > self.too_many_turbines:
            final_tf_list, final_local_dx = self.compress_turbine_function(u,fs)
        else:
            final_tf_list = self.tf_list
            final_local_dx = self.local_dx

        # Combine the force and dx to create a single tf term
        if self.use_local_tf_dx:
            tf_term = 0
            for i in range(len(final_tf_list)):
                tf_term += inner(final_tf_list[i],v)*final_local_dx(i+1)

        # use the full domain of integration if local is not requested
        else:
            tf_term = inner(sum(final_tf_list),v)*dx

        return tf_term

    def compress_turbine_function(self, u, fs):
        """
        used to assemble the turbine force in chunks to avoid recursion errors
        """

        # evenly distribute the number of turbine per project
        n_groups = int((((self.numturbs-1)-((self.numturbs-1)%self.too_many_turbines))+self.too_many_turbines)/self.too_many_turbines)
        split_value = round(self.numturbs/n_groups)
        print(f"Oh dear, there are just too many turbines (>={self.too_many_turbines}). I'm just going to collect them into {n_groups} groups of around {split_value} turbines each and then project.")

        # init loop values
        combined_tf_list = []
        final_local_dx = []
        counter = 0
        group = 1
        
        # Create a cell function that will be used to store the subdomains
        temp_subdomain = MeshFunction("size_t", self.dom.mesh, self.dom.mesh.topology().dim())
        temp_subdomain.set_all(0)

        # disks need to be treated specially due to how tf is convoluted with u
        if "disk" in self.turbine_type:

            # init combined tf
            tf1 = 0
            tf2 = 0
            tf3 = 0


            # sum and project tf parts
            for turb in self.turbines:

                # get the turbine parts
                tf1 += turb.tf1
                tf2 += turb.tf2
                tf3 += turb.tf3

                # the subdomain is a sphere centered at the hub of the turbine with a height of scaling*RD and diameter of scaling*RD
                if self.dom.dim == 3:
                    subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)+(x[2]-HH)*(x[2]-HH)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, HH=turb.z, x0=turb.x, y0=turb.y)
                else:
                    subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, x0=turb.x, y0=turb.y)                

                # add the subdomain to the marker
                subdomain_marker.mark(temp_subdomain,group)

                # once we have collected the maximum number of turbines, combine them an reset the sums
                if counter >= split_value:
                    self.fprint(f"Projecting tf1 for group {group}")
                    proj_tf1 = project(tf1,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
                    self.fprint(f"Projecting tf2 for group {group}")
                    proj_tf2 = project(tf2,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
                    self.fprint(f"Projecting tf3 for group {group}")
                    proj_tf3 = project(tf3,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
                    combined_tf_list.append(proj_tf1*u[0]**2+proj_tf2*u[1]**2+proj_tf3*u[0]*u[1])
                    
                    # reset loop values
                    tf1 = 0
                    tf2 = 0
                    tf3 = 0
                    counter = 0
                    group += 1

                else:
                    counter += 1

            # project the final turbine group
            self.fprint(f"Projecting tf1 for group {group}")
            proj_tf1 = project(tf1,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
            self.fprint(f"Projecting tf2 for group {group}")
            proj_tf2 = project(tf2,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
            self.fprint(f"Projecting tf3 for group {group}")
            proj_tf3 = project(tf3,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
            combined_tf_list.append(proj_tf1*u[0]**2+proj_tf2*u[1]**2+proj_tf3*u[0]*u[1])

            # create a measure for the final turbine group
            final_local_dx = Measure('dx', subdomain_data=temp_subdomain)

        # combine turbine for every other style of turbine
        else:

            # init combined tf
            tf = 0

            # sum and project tf parts
            for turb in self.turbines:

                # get the turbine parts
                tf += turb.tf

                # the subdomain is a sphere centered at the hub of the turbine with a height of scaling*RD and diameter of scaling*RD
                if self.dom.dim == 3:
                    subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)+(x[2]-HH)*(x[2]-HH)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, HH=turb.z, x0=turb.x, y0=turb.y)
                else:
                    subdomain_marker = CompiledSubDomain("(x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0)<=r*r",r=self.local_dx_scaling*turb.RD/2.0, x0=turb.x, y0=turb.y)
                               
                # add the subdomain to the marker
                subdomain_marker.mark(temp_subdomain,group)

                # once we have collected the maximum number of turbines, combine them an reset the sums
                if counter >= split_value:
                    self.fprint(f"Projecting tf for group {group}")
                    proj_tf = project(tf,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
                    combined_tf_list.append(proj_tf)
                    
                    # reset loop values
                    tf = 0
                    counter = 0
                    group += 1
                    temp_subdomain = MeshFunction("size_t", self.dom.mesh, self.dom.mesh.topology().dim())
                    temp_subdomain.set_all(0)
                else:
                    counter += 1

            # project the final turbine group
            self.fprint(f"Projecting tf for group {group}")
            proj_tf = project(tf,fs.V,solver_type='gmres',preconditioner_type="hypre_amg")
            combined_tf_list.append(proj_tf)

            # create a measure for the final turbine group
            final_local_dx = Measure('dx', subdomain_data=temp_subdomain)

        # create the final turbine force          
        return combined_tf_list, final_local_dx

    def update_controls(self):
        """
        iterates over the controls list assigns self.<name> to the control, self.m<name> 
        """
        for turb in self.turbines:
            turb.update_controls()

    def update_turbine_force(self, u, inflow_angle, fs, **kwargs):
        """
        Updates the turbines
        """
        for turb in self.turbines:
            turb.update_turbine_force(u, inflow_angle, fs, **kwargs)

    def compute_power(self, u, inflow_angle):
        """
        Computes the power for the full farm
        """

        ### Build the integrand
        val = 0
        for turb in self.turbines:
            temp = turb.power(u, inflow_angle)
            if not isinstance(temp,(float,AdjFloat)):
                if self.use_local_tf_dx:
                    val += temp*self.local_dx(turb.index+1)
                else:
                    val += temp*dx
            else:
                val += temp

        ### Assemble if needed
        if not isinstance(val,(float,AdjFloat)):
            J = assemble(val)
        else:
            J = val

        return J

    def save_power(self, u, inflow_angle, iter_val = 0.0, simTime = 0.0):
        """
        saves the power output for each turbine
        """

        J_list=np.zeros(self.numturbs+3)
        J_list[0]=iter_val
        J_list[1]=simTime
        with stop_annotating():
            for i,turb in enumerate(self.turbines):
                val = turb.power(u,inflow_angle)

                ### Assemble if needed
                if not isinstance(val,(float,AdjFloat)):
                    if self.use_local_tf_dx:
                        J = assemble(val*self.local_dx(turb.index+1))
                    else:
                        J = assemble(val*dx)
                else:
                    J = val
                    
                J_list[i+2] = J

        J_list[-1]=sum(J_list[2:])

        folder_string = self.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if self.power_first_save:
            header = str("Iter_val    "+"Time    "+"Turbine_%d    "*self.numturbs % tuple(range(self.numturbs))+"Sum")
            self.params.save_csv("power_data",header=header,data=[J_list],subfolder=self.params.folder+"data/",mode='w')
            self.power_first_save = False
        else:
            self.params.save_csv("power_data",data=[J_list],subfolder=self.params.folder+"data/",mode='a')

    def debug_output(self):
        """
        This function computes and save any output needed for regression tests
        """
        if self.debug_mode:

            # Get useful values
            x, y, z = self.get_hub_locations().T
            yaw = self.get_yaw_angles()

            # tag global statistics
            self.tag_output("min_x", np.min(x))
            self.tag_output("max_x", np.max(x))
            self.tag_output("avg_x", np.mean(x))
            self.tag_output("min_y", np.min(y))
            self.tag_output("max_y", np.max(y))
            self.tag_output("avg_y", np.mean(y))
            self.tag_output("min_z", np.min(z))
            self.tag_output("max_z", np.max(z))
            self.tag_output("avg_z", np.mean(z))
            self.tag_output("min_yaw", np.min(yaw))
            self.tag_output("max_yaw", np.max(yaw))
            self.tag_output("avg_yaw", np.mean(yaw))

            # iterate through all turbines
            for turb in self.turbines:
                turb.debug_output()

    def plot_farm(self,show=False,filename="wind_farm",objective_value=None):
        """
        This function plots the locations of each wind turbine and
        saves the output to output/.../plots/

        :Keyword Arguments:
            * **show** (*bool*): Default: True, Set False to suppress output but still save.
        """
        if self.numturbs == 0:
            return

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/"+filename+".pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string) and self.params.rank == 0: os.makedirs(folder_string)

        ### Collect turbine data
        x, y, z = self.get_hub_locations().T
        yaw = self.get_yaw_angles()
        RD = self.get_rotor_diameters()

        ### Create a list that outlines the extent of the farm ###
        if self.optimizing and "layout" in self.control_types and self.layout_bounds != "wind_farm":
            ex_list_x = [self.layout_bounds[0][0],self.layout_bounds[0][1],self.layout_bounds[0][1],self.layout_bounds[0][0],self.layout_bounds[0][0]]
            ex_list_y = [self.layout_bounds[1][0],self.layout_bounds[1][0],self.layout_bounds[1][1],self.layout_bounds[1][1],self.layout_bounds[1][0]]
        else:
            ex_list_x = [self.ex_x[0],self.ex_x[1],self.ex_x[1],self.ex_x[0],self.ex_x[0]]
            ex_list_y = [self.ex_y[0],self.ex_y[0],self.ex_y[1],self.ex_y[1],self.ex_y[0]]

        ### Generate and Save Plot ###
        fig, ax = plt.subplots()
        if hasattr(self.dom,"boundary_line"):
            ax.plot(*self.dom.boundary_line/self.dom.xscale,c="k")
        ax.plot(np.array(ex_list_x)/self.dom.xscale, np.array(ex_list_y)/self.dom.xscale,c="r")

        ### Plot Blades
        for i in range(self.numturbs):
            blade_n = [np.cos(yaw[i]),np.sin(yaw[i])]
            rr = RD[i]/2.0
            blade_x = np.array([x[i]+rr*blade_n[1],x[i]-rr*blade_n[1]])/self.dom.xscale
            blade_y = np.array([y[i]-rr*blade_n[0],y[i]+rr*blade_n[0]])/self.dom.xscale
            ax.plot(blade_x,blade_y,c='k',linewidth=2,zorder=1)

        ### Choose coloring for the turbines ###
        if isinstance(objective_value,(list,np.ndarray)):
            coloring = objective_value
        else:
            coloring = np.array(z)/self.dom.xscale

        ### Plot Hub Locations
        p=ax.scatter(x/self.dom.xscale,y/self.dom.xscale,c=coloring,cmap="coolwarm",edgecolors=(0, 0, 0, 1),s=20,zorder=2)
        # p=plt.scatter(x,y,c="k",s=70)
        plt.xlim(self.dom.x_range[0]/self.dom.xscale,self.dom.x_range[1]/self.dom.xscale)
        plt.ylim(self.dom.y_range[0]/self.dom.xscale,self.dom.y_range[1]/self.dom.xscale)
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('Hub Elevation')

        ### Annotate ###
        for i in range(self.numturbs):
            ax.annotate(i, (x[i]/self.dom.xscale,y[i]/self.dom.xscale),(5,0),textcoords='offset pixels')

        if objective_value is None:
            plt.title("Location of the Turbines")
        elif isinstance(objective_value,(list,np.ndarray)):
            plt.title("Objective Value: {: 5.6f}".format(sum(objective_value)))
        else:
            plt.title("Objective Value: {: 5.6f}".format(float(objective_value)))

        plt.savefig(file_string, transparent=True)

        if show:
            plt.show()

        plt.close()

    def plot_chord(self,show=False,filename="chord_profiles",objective_value=None, bounds=None):
        """
        saves a plot of the chord for each turbine?
        """

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/"+filename+".pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string) and self.params.rank == 0: os.makedirs(folder_string)

        ### Calculate x values ###
        baseline_chord = self.turbines[0].get_baseline_chord()
        baseline_blade_segments = len(baseline_chord)
        x = np.linspace(0,1,baseline_blade_segments)

        ### Plot Chords ###
        plt.figure()
        plt.plot(x,baseline_chord,label="baseline",c="k")
        if bounds is None:
            lower=[]
            upper=[]
            c_avg = 0
            for k, seg_chord in enumerate(baseline_chord):
                    modifier = 2.0
                    max_chord = self.turbines[0].max_chord
                    lower.append(seg_chord/modifier)
                    upper.append(max(min(seg_chord*modifier,max_chord),c_avg))
                    c_avg = (c_avg*k+seg_chord)/(k+1)
            plt.plot(x,lower,"--r",label="Optimization Boundaries")
            plt.plot(x,upper,"--r")
        else:
            bounds = np.array(bounds,dtype=float)
            plt.plot(x,bounds[0][-baseline_blade_segments:],"--r",label="Optimization Boundaries")
            plt.plot(x,bounds[1][-baseline_blade_segments:],"--r")

        for turb in self.turbines:
            x = np.linspace(0,1,turb.num_actuator_nodes)
            y = np.array(turb.chord,dtype=float)
            plt.plot(x,y,'.-',label=turb.index)

        plt.xlim(0,1)
        if objective_value is None:
            plt.title("Chord along blade span")
        elif isinstance(objective_value,(list,np.ndarray)):
            plt.title("Objective Value: {: 5.6f}".format(sum(objective_value)))
        else:
            plt.title("Objective Value: {: 5.6f}".format(float(objective_value))) 
        plt.xlabel("Blade Span")      
        plt.ylabel("Chord")
        plt.legend()      

        plt.savefig(file_string, transparent=True)

        if show:
            plt.show()

        plt.close()

    def save_functions(self,val=0):
        """
        This function call the prepare_saved_functions from each turbine, combines the functions and saves them.
        It then check to see if it can save the function out right and if not it projects. 
        "val" can be the time, or angle, its just the iterator for saving mulitple steps
        Note: this function is way over engineered!
        """

        if self.numturbs >= self.too_many_turbines:
            print(f"Oh dear, there are just too many turbines so saving the turbine force is disabled")
            return

        # gather functions
        func_list = []
        # if self.numturbs>=self.too_many_turbines:
        #     print(f"Oh dear, there are just too many turbines (>={self.too_many_turbines}). Saving turbine function is disabled.")
        # else:
        for turb in self.turbines:
            func_list = turb.prepare_saved_functions(func_list)

        # prepare to store file pointer if needed
        if self.func_first_save:
            self.func_files = []

        # save gathered
        for i, func_to_save in enumerate(func_list):
            func = func_to_save[0]
            func_name = func_to_save[1]

            # Check if function is in a save able form, if not, project ###
            # Note: this is overkill
            if not hasattr(func,"_cpp_object"):

                # choose an appropriate function space 
                if func.geometric_dimension() == 1:
                    FS = self.fs.Q
                else:
                    FS = self.fs.V
                # project onto the function space
                func = project(func,FS,solver_type='cg',preconditioner_type="hypre_amg",**self.extra_kwarg)

            # save, if first time, store the file location pointers
            if self.func_first_save:
                out_file = self.params.Save(func,func_name,subfolder="functions/",val=val)
                self.func_files.append(out_file)
            else:
                self.params.Save(func, func_name,subfolder="functions/",val=val,file=self.func_files[i])
        
        # save the farm level turbine subdomain function
        if self.func_first_save:
            self.subdomain_file = self.params.Save(self.turbine_subdomains,"turbine_subdomains",subfolder="mesh/",val=val,filetype="pvd")
        else:
            self.params.Save(self.turbine_subdomains,"turbine_subdomains",subfolder="mesh/",val=val,file=self.subdomain_file,filetype="pvd")

        # Flip the flag!
        self.func_first_save = False

    def save_wind_farm(self,val=None,filename="wind_farm"):
        """
        saves a text file containing the wind_farm and turbine parameters
        """
        # TODO: update for csv mode
        folder_string = self.params.folder+"/data/"
        if val is not None:
            file_string = self.params.folder+"/data/"+filename+"_"+repr(val)+".txt"
        else:
            file_string = self.params.folder+"/data/"+filename+".txt"

        ### Check if folder exists ###
        if not os.path.exists(folder_string) and self.params.rank == 0: os.makedirs(folder_string)

        ### Define the header string ###
        head_str="#    x    y    HH    Yaw    Diameter    Thickness    Axial_Induction"

        ### Get quantities ###
        x, y, z = self.get_hub_locations().T
        HH = z-self.get_ground_heights()
        yaw = np.degrees(self.get_yaw_angles())
        RD = self.get_rotor_diameters()

        ### Might not need these anymore
        thickness = np.zeros(self.numturbs)
        axial = np.zeros(self.numturbs)
        for i,turb in enumerate(self.turbines):
            if 'line' in turb.type:
                thickness[i] = np.nan
                axial[i] = np.nan
            elif turb.type == 'disabled':
                thickness[i] = np.nan
                axial[i] = np.nan
            else:
                thickness[i] = turb.thickness
                axial[i] = turb.axial

        ### Save text file ###
        output = np.array([x, y, HH, yaw, RD, thickness, axial])
        np.savetxt(file_string,output.T,header=head_str)

    def finalize_farm(self):
        for k, turb in enumerate(self.turbines):
            turb.finalize_turbine()




















    ########################################################################################################
    ############## It would be smart to eventually move these to a separate refinement module ##############
    ########################################################################################################
 
    def SimpleRefine(self,radius,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Cylinder Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_locs = self.get_hub_locations()
        turb_x = turb_locs[:,0]
        turb_y = turb_locs[:,1]
        if self.dom.dim == 3:
            turb_z0 = self.dom.z_range[0]-radius
            turb_z1 = turb_locs[:,2]+radius

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Find the minimum distance for each turbine with the vertices ###
            x_diff = np.power(np.subtract.outer(x,turb_x),2.0)
            y_diff = np.power(np.subtract.outer(y,turb_y),2.0)
            min_r = np.min(x_diff+y_diff,axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_circle = min_r <= radius**2.0

            ### Determine if in cylinder ###
            if d == 3:
                in_z = np.logical_and(turb_z0 <= max(z), turb_z1 >= min(z))
                near_turbine = np.logical_and(in_circle, in_z)
            else:
                near_turbine = in_circle

            ### mark if cell is near any cylinder ###
            if any(near_turbine):
                cell_f[cell] = True

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.update_heights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")


    def WakeRefine(self,radius,length,theta=0.0,expand_factor=1,centered=False):
        if self.numturbs == 0:
            return

        self.fprint("Wake Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius/2.0
        length = length+2*(expand_factor-1)*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)

        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_locs = self.get_hub_locations()
        turb_x = turb_locs[:,0]
        turb_y = turb_locs[:,1]
        if self.dom.dim == 3:
            turb_z = turb_locs[:,2]

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Rotate the Cylinder about the turbine axis ###
            x_diff = np.subtract.outer(x,turb_x)
            y_diff = np.subtract.outer(y,turb_y)
            x = (np.cos(theta)*(x_diff)-np.sin(theta)*(y_diff) + turb_x)
            y = (np.sin(theta)*(x_diff)+np.cos(theta)*(y_diff) + turb_y)

            ### Determine if in wake ###
            if centered:
                # Center the refinement region around the turbine
                # upstream and downstream by length/2
                x0 = turb_x - length/2.0
                x1 = turb_x + length/2.0
            else:
                # Otherwise, refine the default amount upstream (R)
                # and the full length amount downstreeam
                x0 = turb_x - radius
                x1 = turb_x + length

            in_wake = np.logical_and(np.greater(x,x0),np.less(x,x1))
            in_wake = np.any(in_wake,axis=0)

            ### Find the minimum distance for each turbine with the vertices ###
            y_diff = y-turb_y
            if d == 3:
                z_diff = np.subtract.outer(z,turb_z)
                min_r = np.min(np.power(y_diff,2.0)+np.power(z_diff,2.0),axis=0)
            else:
                min_r = np.min(np.power(y_diff,2.0),axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_circle = min_r <= radius**2.0
            near_turbine = np.logical_and(in_circle, in_wake)

            ### mark if cell is near any cylinder ###
            if any(near_turbine):
                cell_f[cell] = True


        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.update_heights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def TearRefine(self,radius,theta=0.0,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Tear Drop Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_locs = self.get_hub_locations()
        turb_x = turb_locs[:,0]
        turb_y = turb_locs[:,1]
        if self.dom.dim == 3:
            turb_z0 = self.dom.z_range[0]-radius
            turb_z1 = turb_locs[:,2]+radius

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Rotate the Cylinder about the turbine axis ###
            x_diff = np.subtract.outer(x,turb_x)
            y_diff = np.subtract.outer(y,turb_y)
            x = (np.cos(theta)*(x_diff)-np.sin(theta)*(y_diff) + turb_x)
            y = (np.sin(theta)*(x_diff)+np.cos(theta)*(y_diff) + turb_y)
            x_diff = x - turb_x
            y_diff = y - turb_y

            ### Find Closest Turbine ###
            min_dist = np.min(np.power(x_diff,2.0)+np.power(y_diff,2.0),axis=0)
            min_turb_id = np.argmin(min_dist)

            # ### Do something based on upstream or downstream ###
            # if min(x[:,min_turb_id]-turb_x[min_turb_id]) <= 0:
            #     val = -min_turb_id-1
            # else:
            #     val = min_turb_id+1

            # ### Determine if in z_range ###
            # if d == 3:
            #     in_z = turb_z0 <= max(z) and turb_z1[min_turb_id] >= min(z)
            #     if in_z:
            #         near_turbine = val
            #     else:
            #         near_turbine = 0
            # else:
            #     near_turbine = val

            ### Check if upstream or downstream and adjust accordingly ###
            cl_x = turb_x[min_turb_id]
            cl_y = turb_y[min_turb_id]
            dx = min(x[:,min_turb_id]-cl_x)
            dy = min(y[:,min_turb_id]-cl_y)
            if dx <= 0:
                dist = (dx*1.5)**2 + dy**2
            else:
                dist = (dx/2)**2 + dy**2

            ### determine if in z-range ###
            if d == 3:
                in_z = turb_z0 <= max(z) and turb_z1[min_turb_id] >= min(z)
                near_turbine = in_z and dist <= radius**2
            else:
                near_turbine = dist <= radius**2
            ### mark if cell is near any cylinder ###
            if near_turbine:
                cell_f[cell] = True


        # File("test.pvd") << cell_f
        # exit()

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.update_heights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def SphereRefine(self,radius,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Sphere Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_locs = self.get_hub_locations()
        turb_x = turb_locs[:,0]
        turb_y = turb_locs[:,1]
        if self.dom.dim == 3:
            turb_z = turb_locs[:,2]
            

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Find the minimum distance for each turbine with the vertices ###
            min_r  = np.power(np.subtract.outer(x,turb_x),2.0)
            min_r += np.power(np.subtract.outer(y,turb_y),2.0)
            if d == 3:
                min_r += np.power(np.subtract.outer(z,turb_z),2.0)
            min_r = np.min(min_r,axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_sphere = min_r <= radius**2.0

            ### mark if cell is near any cylinder ###
            if any(in_sphere):
                cell_f[cell] = True

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.update_heights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

        