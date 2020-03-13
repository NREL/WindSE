"""
The DomainManager submodule contains the various classes used for 
creating different types of domains
"""

import __main__
import os

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    from mshr import *
    import copy
    import time
    import warnings
    import os
    from sys import platform
    import numpy as np
    from scipy.interpolate import interp2d, interp1d,RectBivariateSpline
    import shutil

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *
        
    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

    ### This parameter allows for refining the mesh functions ###
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"


def Elliptical_Grid(x, y, z, radius):
    #x_hat = x sqrt(1 - y^2/2)
    #y_hat = y sqrt(1 - x^2/2)
    x_hat = np.multiply(radius*x,np.sqrt(1.0-np.power(y,2.0)/2.0))
    y_hat = np.multiply(radius*y,np.sqrt(1.0-np.power(x,2.0)/2.0))
    z_hat = z
    return [x_hat, y_hat, z_hat]

def FG_Squircular(x, y, z, radius):
    #x_hat = x sqrt(x^2 + y^2 - x^2y^2) / sqrt(x^2 + y^2)
    #y_hat = y sqrt(x^2 + y^2 - x^2y^2) / sqrt(x^2 + y^2)
    innerp = np.power(x,2.0)+np.power(y,2.0)
    prod  = np.multiply(np.power(x,2.0),np.power(y,2.0))
    innerp[innerp==0] = 1 #handle the point (0,0)
    ratio = np.divide(np.sqrt(np.subtract(innerp,prod)),np.sqrt(innerp))

    x_hat = np.multiply(radius*x,ratio)
    y_hat = np.multiply(radius*y,ratio)
    z_hat = z
    return [x_hat, y_hat, z_hat]

def Simple_Stretching(x, y, z, radius):
    radii = np.sqrt(np.power(x,2.0)+np.power(y,2.0))
    radii[radii==0] = 1 #handle the point (0,0)
    prod  = np.multiply(x,y)
    x2_ratio = np.divide(np.power(x,2.0),radii)
    y2_ratio = np.divide(np.power(y,2.0),radii)
    xy_ratio = np.divide(prod,radii)

    x2_gte_y2 = np.power(x,2.0)>=np.power(y,2.0)

    x_hat = np.zeros(len(x))
    y_hat = np.zeros(len(y))

    x_hat[x2_gte_y2]  = np.multiply(radius*np.sign(x),x2_ratio)[x2_gte_y2]
    x_hat[~x2_gte_y2] = np.multiply(radius*np.sign(y),xy_ratio)[~x2_gte_y2]

    y_hat[x2_gte_y2]  = np.multiply(radius*np.sign(x),xy_ratio)[x2_gte_y2]
    y_hat[~x2_gte_y2] = np.multiply(radius*np.sign(y),y2_ratio)[~x2_gte_y2]

    z_hat = z
    return [x_hat, y_hat, z_hat]




class GenericDomain(object):
    """
    A GenericDomain contains on the basic functions required by all domain objects
    """

    def __init__(self):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.first_save = True
        self.finalized = False
        self.fprint = self.params.fprint

    def Plot(self):
        """
        This function plots the domain using matplotlib and saves the 
        output to output/.../plots/mesh.pdf
        """

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/mesh.pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        p=plot(self.mesh)
        plt.savefig(file_string)
        plt.show()

    def Save(self,val=0):
        """
        This function saves the mesh and boundary markers to output/.../mesh/
        """

        # :Keyword Arguments:
        #     * **filename** (*str*): the file name that preappends the meshes
        #     * **filetype** (*str*): the file type to save (pvd or xml)
        # folder_string = self.params.folder+"/mesh/"
        # mesh_string = self.params.folder+"/mesh/"+filename+"_mesh."+filetype
        # bc_string = self.params.folder+"/mesh/"+filename+"_boundaries."+filetype

        # ### Check if folder exists ###
        # if not os.path.exists(folder_string): os.makedirs(folder_string)

        # print("Saving Mesh")
        # ### Save Mesh ###
        # file = File(mesh_string)
        # file << self.mesh

        # ### Save Boundary Function ###
        # file = File(bc_string)
        # file << self.boundary_markers
        # print("Mesh Saved")

        if self.first_save:
            self.mesh_file = self.params.Save(self.mesh,"mesh",subfolder="mesh/",val=val,filetype="pvd")
            self.bmesh_file   = self.params.Save(self.bmesh,"boundary_mesh",subfolder="mesh/",val=val,filetype="pvd")
            self.bc_file   = self.params.Save(self.boundary_markers,"facets",subfolder="mesh/",val=val,filetype="pvd")
            # self.mr_file   = self.params.Save(self.mesh_radius,"mesh_radius",subfolder="mesh/",val=val,filetype="pvd")
            self.first_save = False
        else:
            self.params.Save(self.mesh,"mesh",subfolder="mesh/",val=val,file=self.mesh_file,filetype="pvd")
            self.params.Save(self.bmesh,"boundary_mesh",subfolder="mesh/",val=val,file=self.bmesh_file,filetype="pvd")
            self.params.Save(self.boundary_markers,"facets",subfolder="mesh/",val=val,file=self.bc_file,filetype="pvd")
            # self.params.Save(self.mesh_radius,"mesh_radius",subfolder="mesh/",val=val,file=self.mr_file,filetype="pvd")

    def Refine(self,num,region=None,region_type=None,cell_markers=None):
        """
        This function can be used to refine the mesh. If a region is
        specified, the refinement is local

        Args:
            num (int): the number of times to refine

        :Keyword Arguments:
            * **region** (*list*): 
                                for square region use: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
                                for circle region use: [[radius],[c_x,c_y],[zmin,zmax]]
            * **region_type** (*str*): Either "circle" or "square" 
            * **cell_markers** (*:meth:dolfin.mesh.MeshFunction*): A cell function marking which cells to refine 
        """
        refine_start = time.time()

        ### Print some useful stats 
        self.fprint("Starting Mesh Refinement",special="header")
        if cell_markers is not None:
            self.fprint("Region Type: {0}".format("cell_markers"))
        elif region is not None:
            self.fprint("Region Type: {0}".format(region_type))
            if "circle" in region_type:
                self.fprint("Circle Radius: {:.2f}".format(region[0][0]))
                self.fprint("Circle Center: ({:.2f}, {:.2f})".format(region[1][0],region[1][1]))
                if self.dim == 3:
                    self.fprint("Z Range:       [{:.2f}, {:.2f}]".format(region[2][0],region[2][1]))
            else:
                self.fprint("X Range: [{: .2f}, {: .2f}]".format(region[0][0],region[0][1]))
                self.fprint("Y Range: [{: .2f}, {: .2f}]".format(region[1][0],region[1][1]))
                if self.dim == 3:
                    self.fprint("Z Range: [{: .2f}, {: .2f}]".format(region[2][0],region[2][1]))
        else:
            self.fprint("Region Type: {0}".format("full"))

        ### Mark cells for refinement
        for i in range(num):
            if num>1:
                step_start = time.time()
                self.fprint("Refining Mesh Step {:d} of {:d}".format(i+1,num), special="header")
            else:
                self.fprint("")


            ### Check if cell markers were provided ###
            if cell_markers is not None:
                cell_f = cell_markers
                self.fprint("Cells Marked for Refinement: {:d}".format(sum(cell_markers.array())))

            ### Check if a region was provided ###
            elif region is not None:
                self.fprint("Marking Cells")

                ### Create an empty cell marker function
                cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
                cells_marked = 0

                ### Check if we are refining in a circle ###
                if "circle" in region_type:
                    radius=region[0][0]
                    for cell in cells(self.mesh):
                        in_circle = (cell.midpoint()[0]-region[1][0])**2.0+(cell.midpoint()[1]-region[1][1])**2.0<=radius**2.0
                        if self.dim == 3:
                            # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                            g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                            in_z = between(cell.midpoint()[2],(region[2][0],region[2][1]+g_z))
                            in_region = in_circle and in_z
                        else:
                            in_region = in_circle
                        if in_region:
                            cell_f[cell] = True
                            cells_marked += 1

                ### or a rectangle ###
                else:
                    cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
                    for cell in cells(self.mesh):
                        in_square = between(cell.midpoint()[0],tuple(region[0])) and between(cell.midpoint()[1],tuple(region[1]))
                        if self.dim == 3:
                            # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                            g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                            in_z = between(cell.midpoint()[2],(region[2][0],region[2][1]+g_z))
                            in_region = in_square and in_z
                        else:
                            in_region = in_square
                        if in_region:
                            cell_f[cell] = True
                            cells_marked += 1
                self.fprint("Cells Marked for Refinement: {:d}".format(cells_marked))

            
            ### If neither a region or cell markers were provided, Refine everwhere ###
            else:
                cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),True)

            old_verts = self.mesh.num_vertices()
            old_cells = self.mesh.num_cells()
            self.fprint("Refining Mesh")
            self.mesh = refine(self.mesh,cell_f)
            self.bmesh = BoundaryMesh(self.mesh,"exterior")
            self.boundary_markers = adapt(self.boundary_markers,self.mesh)

            self.fprint("Original Mesh Vertices: {:d}".format(old_verts))
            self.fprint("Original Mesh Cells:    {:d}".format(old_cells))
            self.fprint("New Mesh Vertices:      {:d}".format(self.mesh.num_vertices()))
            self.fprint("New Mesh Cells:         {:d}".format(self.mesh.num_cells()))
            if num>1:
                step_stop = time.time()
                self.fprint("Step {:d} of {:d} Finished: {:1.2f} s".format(i+1,num,step_stop-step_start), special="footer")
        
        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")


    def WarpSplit(self,h,s):
        """
        This function warps the mesh to shift more cells towards the ground. 
        is achieved by spliting the domain in two and moving the cells so 
        that a percentage of them are below the split.

        Args:
            h (float): the height that split occurs
            s (float): the percent below split in the range [0,1)
        """

        warp_start = time.time()
        self.fprint("Starting Mesh Warping",special="header")
        self.fprint("Height of Split:     {:1.2f} m".format(h))
        self.fprint("Percent below Split: {:1.2f}%".format(s*100.0))

        if self.mesh.topology().dim() == 3:
            z0 = self.z_range[0]
            z1 = self.z_range[1]
            z_ind = 2
        else:
            z0 = self.y_range[0]
            z1 = self.y_range[1]
            z_ind = 1

        # h1 = 60
        # h2 = 160

        z = copy.deepcopy(self.mesh.coordinates()[:,z_ind])
        # cubic_spline = interp1d([z0,a-r,a+r,z1],[z0,a-(1-s)*r,a+(1-s)*r,z1])
        cubic_spline = interp1d([z0,z0+s*(z1-z0),z1],[z0,h,z1],fill_value=(z0,z1),bounds_error=False)
        # cubic_spline = interp1d([z0,h1-s*(z0+h1),h2+s*(z1-h2),z1],[z0,h1,h2,z1],fill_value=(z0,z1),bounds_error=False)

        # plt.figure()
        # x = np.linspace(z0,z1,100)
        # y = cubic_spline(x)
        # plt.plot(x,y)
        # plt.show()
        # exit()

        self.fprint("Moving Nodes")
        z_new = cubic_spline(z)
        self.mesh.coordinates()[:,z_ind]=z_new
        self.mesh.bounding_box_tree().build(self.mesh)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")

        warp_stop = time.time()
        self.fprint("Warping Finished: {:1.2f} s".format(warp_stop-warp_start),special="footer")


    def WarpSmooth(self,s):
        """
        This function warps the mesh to shift more cells towards the ground. 
        The cells are shifted based on the function:

        .. math::

            z_new = z_0 + (z_1-z_0) \\left( \\frac{z_old-z_0}{z_1-z_0} \\right)^{s}.

        where :math:`z_0` is the ground and :math:`z_1` is the top of the domain.

        Args:
            s (float): compression strength
        """
        warp_start = time.time()
        self.fprint("Starting Mesh Warping",special="header")
        self.fprint("Compression Strength: {:1.4f}".format(s))
        self.fprint("Moving Nodes")

        z=self.mesh.coordinates()[:,2].copy()
        z0 = self.z_range[0]
        z1 = self.z_range[1]
        z1 = z0 + (z1 - z0)*((z-z0)/(z1-z0))**s
        self.mesh.coordinates()[:,2]=z1
        self.mesh.bounding_box_tree().build(self.mesh)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")

        warp_stop = time.time()
        self.fprint("Warping Finished: {:1.2f} s".format(warp_stop-warp_start),special="footer")

    def Move(self,ground):
        move_start = time.time()
        self.fprint("Moving Mesh",special="header")

        def transform(x,y,z,z0,z1):
            new_z = np.zeros(len(z))
            for i in range(len(z)):
                new_z[i]=(z1-ground(x[i],y[i]))*(z[i]-z0)/(z1-z0)+ground(x[i],y[i])
            return new_z

        self.fprint("Moving Boundary Mesh")
        x_hd = copy.deepcopy(self.bmesh.coordinates()[:,0])
        y_hd = copy.deepcopy(self.bmesh.coordinates()[:,1])
        z_hd = copy.deepcopy(self.bmesh.coordinates()[:,2])
        z_hd = transform(x_hd,y_hd,z_hd,self.z_range[0],self.z_range[1])
        self.bmesh.coordinates()[:,0]=x_hd
        self.bmesh.coordinates()[:,1]=y_hd
        self.bmesh.coordinates()[:,2]=z_hd
        self.bmesh.bounding_box_tree().build(self.bmesh)
        
        self.fprint("Moving Mesh to New Boundary Using ALE")
        ALE.move(self.mesh,self.bmesh)

        move_stop = time.time()
        self.fprint("Mesh Moved: {:1.2f} s".format(move_stop-move_start),special="footer")

    def Finalize(self):
        # self.ComputeCellRadius()
        self.finalized = True

    # def ComputeCellRadius(self):
    #     self.mesh_radius = MeshFunction("double", self.mesh, self.mesh.topology().dim())
    #     for cell in cells(self.mesh):
    #         self.mesh_radius.set_value(cell.index(),cell.h())
    #         # self.mesh_radius.set_value(cell.index(),cell.inradius())
    #         # self.mesh_radius.set_value(cell.index(),cell.circumradius())

    def SetupInterpolatedGround(self):
        self.fprint("Ground Type: Interpolated From File")

        ### Import data from Options ###
        if "path" in self.params["domain"]:
            self.path = self.params["domain"]["path"]
            self.terrain_path  = self.path + "topography.txt"
        else:
            self.terrain_path  = self.params["domain"]["terrain_path"]

        ### Copy Files to input folder ###
        shutil.copy(self.terrain_path,self.params.folder+"input_files/")

        self.fprint("Path: {0}".format(self.terrain_path),offset=1)

        ### import ground data
        self.topography = np.loadtxt(self.terrain_path)
        x_data = self.topography[1:,0]
        y_data = self.topography[1:,1]
        z_data = self.topography[1:,2]

        ### generate interpolating function
        x_data = np.sort(np.unique(x_data))
        y_data = np.sort(np.unique(y_data))
        z_data = np.reshape(z_data,(int(self.topography[0,0]),int(self.topography[0,1])))
        self.topography_interpolated = RectBivariateSpline(x_data,y_data,z_data.T)

    def SetupAnalyticGround(self):
        self.hill_sigma_x = self.params["domain"]["gaussian"]["sigma_x"]
        self.hill_sigma_y = self.params["domain"]["gaussian"]["sigma_y"]
        self.hill_theta = self.params["domain"]["gaussian"].get("theta",0.0)
        self.hill_amp = self.params["domain"]["gaussian"]["amp"]
        self.hill_center = self.params["domain"]["gaussian"].get("center",[0.0,0.0])
        self.hill_x0 = self.hill_center[0]
        self.hill_y0 = self.hill_center[1]
        self.fprint("")
        self.fprint("Ground Type: Gaussian Hill")
        self.fprint("Hill Center:   ({: .2f}, {: .2f})".format(self.hill_x0,self.hill_y0),offset=1)
        self.fprint("Hill Rotation:  {: <7.2f}".format(self.hill_theta),offset=1)
        self.fprint("Hill Amplitude: {: <7.2f}".format(self.hill_amp),offset=1)
        self.fprint("Hill sigma_x:   {: <7.2f}".format(self.hill_sigma_x),offset=1)
        self.fprint("Hill sigma_y:   {: <7.2f}".format(self.hill_sigma_y),offset=1)
        self.hill_a = np.cos(self.hill_theta)**2/(2*self.hill_sigma_x**2) + np.sin(self.hill_theta)**2/(2*self.hill_sigma_y**2)
        self.hill_b = np.sin(2*self.hill_theta)/(4*self.hill_sigma_y**2) - np.sin(2*self.hill_theta)/(4*self.hill_sigma_x**2)
        self.hill_c = np.cos(self.hill_theta)**2/(2*self.hill_sigma_y**2) + np.sin(self.hill_theta)**2/(2*self.hill_sigma_x**2)

    def GaussianGroundFuncion(self,x,y,dx=0,dy=0):
        return self.hill_amp*exp( - (self.hill_a*(x-self.hill_x0)**2 + 2*self.hill_b*(x-self.hill_x0)*(y-self.hill_y0) + self.hill_c*(y-self.hill_y0)**2)**2)+self.z_range[0]

    def InterplatedGroundFunction(self,x,y,dx=0,dy=0):
        if dx == 0 and dy == 0:
            return float(self.topography_interpolated(x,y)[0]+self.z_range[0])
        else:
            return float(self.topography_interpolated(x,y,dx=dx,dy=dy)[0])


    def Ground(self,x,y,dx=0,dy=0):
        """
        Ground returns the ground height given an (*x*, *y*) coordinate.

        Args:
            x (float/list): *x* location within the domain
            y (float/list): *y* location within the domain

        Returns:
            float/list: corresponding z coordinates of the ground.

        """

        if isinstance(x,Constant):
            z = self.ground_function(x,y,dx=dx,dy=dy)
            return z
        else:
            if (isinstance(x,list) and isinstance(y,list)) or (isinstance(x,np.ndarray) and isinstance(y,np.ndarray)):
                nx = len(x)
                ny = len(y)
                if nx != ny:
                    raise ValueError("Length mismatch: len(x)="+repr(nx)+", len(y)="+repr(ny))
                else:
                    z = np.zeros(nx)
                    for i in range(nx):
                        z[i] = float(self.ground_function(x[i],y[i],dx=dx,dy=dy))
                    return z
            else:
                return float(self.ground_function(x,y,dx=dx,dy=dy))

class BoxDomain(GenericDomain):
    """
    A box domain is simply a 3D rectangular prism. This box is defined
    by 6 parameters in the param.yaml file. 

    Example:
        In the yaml file define::

            domain: 
                #                      # Description           | Units
                x_range: [-2500, 2500] # x-range of the domain | m
                y_range: [-2500, 2500] # y-range of the domain | m
                z_range: [0.04, 630]   # z-range of the domain | m
                nx: 10                 # Number of x-nodes     | -
                ny: 10                 # Number of y-nodes     | -
                nz: 2                  # Number of z-nodes     | -

        This will produce a box with corner points (-2500,-2500,0.04) 
        to (2500,2500,630). The mesh will have *nx* nodes in the *x*-direction,
        *ny* in the *y*-direction, and *nz* in the *z*-direction.
    """

    def __init__(self):
        super(BoxDomain, self).__init__()

        self.fprint("Generating Box Domain",special="header")

        ### Initialize values from Options ###
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)
        self.x_range = self.params["domain"]["x_range"]
        self.y_range = self.params["domain"]["y_range"]
        self.z_range = self.params["domain"]["z_range"]
        self.nx = self.params["domain"]["nx"]
        self.ny = self.params["domain"]["ny"]
        self.nz = self.params["domain"]["nz"]
        self.dim = 3

        ### Get the initial wind direction ###
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)

        ### Print Some stats ###
        self.fprint("X Range: [{: .2f}, {: .2f}]".format(self.x_range[0],self.x_range[1]))
        self.fprint("Y Range: [{: .2f}, {: .2f}]".format(self.y_range[0],self.y_range[1]))
        self.fprint("Z Range: [{: .2f}, {: .2f}]".format(self.z_range[0],self.z_range[1]))

        ### Create mesh ###
        mesh_start = time.time()
        self.fprint("")
        self.fprint("Generating Mesh")
        start = Point(self.x_range[0], self.y_range[0], self.z_range[0])
        stop  = Point(self.x_range[1], self.y_range[1], self.z_range[1])
        self.mesh = BoxMesh(start, stop, self.nx, self.ny, self.nz)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.fprint("Mesh Generated: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        top     = CompiledSubDomain("near(x[2], z1, tol) && on_boundary",z1 = self.z_range[1], tol = 1e-10)
        bottom  = CompiledSubDomain("near(x[2], z0, tol) && on_boundary",z0 = self.z_range[0], tol = 1e-10)
        front   = CompiledSubDomain("near(x[0], x0, tol) && on_boundary",x0 = self.x_range[0], tol = 1e-10)
        back    = CompiledSubDomain("near(x[0], x1, tol) && on_boundary",x1 = self.x_range[1], tol = 1e-10)
        left    = CompiledSubDomain("near(x[1], y0, tol) && on_boundary",y0 = self.y_range[0], tol = 1e-10)
        right   = CompiledSubDomain("near(x[1], y1, tol) && on_boundary",y1 = self.y_range[1], tol = 1e-10)
        self.boundary_subdomains = [top,bottom,front,back,left,right]
        self.boundary_names = {"top":1,"bottom":2,"front":3,"back":4,"left":5,"right":6}
        self.boundary_types = {"inflow":          ["front","left","right"],
                               "no_slip":         ["bottom"],
                               "horizontal_slip": ["top"],
                               "no_stress":       ["back"]}

        ### Generate the boundary markers for boundary conditions ###
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            self.boundary_subdomains[i].mark(self.boundary_markers, i+1)
        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))


        self.fprint("Initial Domain Setup",special="footer")

    def ground_function(self,x,y,dx=0,dy=0):
        if dx == 0 and dy == 0:
            #################
            #################
            #################
            #################
            #################
            return self.z_range[0]
            # return 0.0
            #################
            #################
            #################
            #################
            #################
        else:
            return 0.0

    def RecomputeBoundaryMarkers(self,theta):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        tol = 1e-3

        ### This function rounds to the nearest ordinal direction ###
        theta45 = round((theta-pi/4.0)/(pi/2))*(pi/2)+pi/4.0

        ### Check if the wind angle is a ordinal direction ###
        if   near(theta45, 1.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","left"]
            self.boundary_types["no_stress"] = ["back","right"]
        elif near(theta45, 3.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","left"]
            self.boundary_types["no_stress"] = ["front","right"]
        elif near(theta45, 5.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","right"]
            self.boundary_types["no_stress"] = ["front","left"]
        elif near(theta45, 7.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["top","front","right"]
            self.boundary_types["no_stress"] = ["back","left"]

        ### Check the special cases that the wind angle is a cardinal direction ###
        if   near(theta, 0.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","right","left"]
            self.boundary_types["no_stress"] = ["back"]
        elif near(theta, 1.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","back","left"]
            self.boundary_types["no_stress"] = ["right"]
        elif near(theta, 2.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","right","left"]
            self.boundary_types["no_stress"] = ["front"]
        elif near(theta, 3.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","back","right"]
            self.boundary_types["no_stress"] = ["left"]

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))

class CylinderDomain(GenericDomain):
    """
    A cylinder domain is a cylinder that is centered a c0 and has radius r. 
    This domain is defined by 6 parameters in the param.yaml file. The 
    center of the cylinder is assumed to be the z-axis.

    Example:
        In the yaml file define::

            domain: 
                #                      # Description           | Units
                z_range: [0.04, 630]   # z-range of the domain | m
                radius: 2500           # radius of base circle | m
                nt: 100                # Number of radial nodes| -
                nz: 10                 # Number of z nodes     | -

        This will produce a upright cylinder centered at (0.0,0.0) with a 
        radius of 2500 m and extends from z=0.04 to 630 m. The mesh will 
        have *nx* nodes in the *x*-direction, *ny* in the *y*-direction, and
        *nz* in the *z*-direction.
    """

    def __init__(self):
        super(CylinderDomain, self).__init__()

        self.fprint("Generating Cylinder Domain",special="header")

        ### Initialize values from Options ###
        self.radius   = self.params["domain"]["radius"]
        self.center   = self.params["domain"]["center"]
        self.z_range  = self.params["domain"]["z_range"]
        self.nt = self.params["domain"]["nt"]
        self.mesh_type = self.params["domain"].get("mesh_type","mshr")
        self.x_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.y_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.dim = 3

        ### Get the initial wind direction ###
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)

        ### Calculating the boundary of the shadow ###
        angles = np.linspace(0,2.0*np.pi,self.nt+1)
        self.boundary_line = (self.radius*np.cos(angles)+self.center[0],self.radius*np.sin(angles)+self.center[1])

        self.fprint("Radius:        {: .2f}".format(self.radius))
        self.fprint("Center:       ({: .2f}, {: .2f})".format(self.center[0],self.center[1]))
        self.fprint("Z Range:      [{: .2f}, {: .2f}]".format(self.z_range[0],self.z_range[1]))
        self.fprint("Meshing Type:  {0}".format(self.mesh_type))

        mesh_start = time.time()
        self.fprint("")
        if self.mesh_type == "mshr":
            self.fprint("Generating Mesh Using mshr")

            self.res = self.params["domain"]["res"]

            ### Create Mesh ###
            # mshr_circle = Circle(Point(self.center[0],self.center[1]), self.radius, self.nt)
            # mshr_domain = Extrude2D(mshr_circle,self.z_range[1]-self.z_range[0])
            top    = Point(self.center[0],self.center[1],self.z_range[1])
            bottom = Point(self.center[0],self.center[1],self.z_range[0])
            mshr_domain = Cylinder(top,bottom,self.radius,self.radius,self.nt)
            self.mesh = generate_mesh(mshr_domain,self.res)
            self.mesh = refine(self.mesh)
            self.mesh = refine(self.mesh)
            # self.mesh = refine(self.mesh)

            # z = self.mesh.coordinates()[:,2]#+self.z_range[0]
            # self.mesh.coordinates()[:,2] = z
            # self.mesh.bounding_box_tree().build(self.mesh)

        else:
            self.fprint("Generating Box Mesh")

            self.nz = self.params["domain"]["nz"]
            self.nxy = int(self.nt/4.0)

            ### Create mesh ###
            start = Point(-1.0, -1.0, self.z_range[0])
            stop  = Point( 1.0,  1.0, self.z_range[1])
            self.mesh = BoxMesh(start, stop, int(self.nxy/2.0), int(self.nxy/2.0), int(self.nz/2.0))
            self.mesh = refine(self.mesh)
            x = self.mesh.coordinates()[:,0]
            y = self.mesh.coordinates()[:,1]
            z = self.mesh.coordinates()[:,2]
            
            self.fprint("Morphing Mesh")
            if self.mesh_type == "elliptic":
                x_hat, y_hat, z_hat = Elliptical_Grid(x, y, z, self.radius)
            elif self.mesh_type == "squircular":
                x_hat, y_hat, z_hat = FG_Squircular(x, y, z, self.radius)
            elif self.mesh_type == "stretch":
                x_hat, y_hat, z_hat = Simple_Stretching(x, y, z, self.radius)

            x_hat += self.center[0]
            y_hat += self.center[1]

            xy_hat_coor = np.array([x_hat, y_hat, z_hat]).transpose()
            self.mesh.coordinates()[:] = xy_hat_coor
            self.mesh.bounding_box_tree().build(self.mesh)

        ### Create the boundary mesh ###
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.fprint("Mesh Generated: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Define Plane Normal ###
        nom_x = np.cos(self.init_wind)
        nom_y = np.sin(self.init_wind)

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1])
        inflow  = CompiledSubDomain("x[2] >= z0 && x[2] <= z1 && nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1], c0=self.center[0], c1=self.center[1])
        top     = CompiledSubDomain("near(x[2], z1) && on_boundary",z1 = self.z_range[1])
        bottom  = CompiledSubDomain("near(x[2], z0) && on_boundary",z0 = self.z_range[0])
        self.boundary_subdomains = [outflow,inflow,top,bottom]
        self.boundary_names = {"inflow":2,"outflow":1,"top":3,"bottom":4}
        self.boundary_types = {"inflow":          ["inflow"],
                               "no_stress":       ["outflow"],
                               "horizontal_slip": ["top"],
                               "no_slip":         ["bottom"]}

        ### Generate the boundary markers for boundary conditions ###
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            self.boundary_subdomains[i].mark(self.boundary_markers, i+1,check_midpoint=False)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def ground_function(self,x,y,dx=0,dy=0):
        if dx == 0 and dy == 0:
            return self.z_range[0]
        else:
            return 0.0

    def RecomputeBoundaryMarkers(self,theta):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Define Plane Normal ###
        nom_x = np.cos(theta)
        nom_y = np.sin(theta)

        # ### Define Boundary Subdomains ###
        # outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1])
        # inflow  = CompiledSubDomain("x[2] >= z0 && x[2] <= z1 && nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1], c0=self.center[0], c1=self.center[1])
        # top     = CompiledSubDomain("near(x[2], z1) && on_boundary",z1 = self.z_range[1])
        # bottom  = CompiledSubDomain("near(x[2], z0) && on_boundary",z0 = self.z_range[0])
        # self.boundary_subdomains = [outflow,inflow,top,bottom]

        # ### Generate the boundary markers for boundary conditions ###
        # self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        # self.boundary_markers.set_all(0)
        # for i in range(len(self.boundary_subdomains)):
        #     self.boundary_subdomains[i].mark(self.boundary_markers, i+1,check_midpoint=False)

        ### Define center ###
        c0 = self.center[0]
        c1 = self.center[1]

        ### Set Tol ###
        tol = 1e-5

        wall_facets = self.boundary_markers.where_equal(self.boundary_names["inflow"]) \
                    + self.boundary_markers.where_equal(self.boundary_names["outflow"])

        boundary_val_temp = self.boundary_markers.array()
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_values(boundary_val_temp)

        for facet_id in wall_facets:
            facet = Facet(self.mesh,facet_id)
            vert_ids = facet.entities(0)
            vert_coords = self.mesh.coordinates()[vert_ids]
            x = vert_coords[:,0]
            y = vert_coords[:,1]

            if all(nom_x*(x-c0)+nom_y*(y-c1)<=0+tol):
                self.boundary_markers.set_value(facet_id,self.boundary_names["inflow"])
            else:
                self.boundary_markers.set_value(facet_id,self.boundary_names["outflow"])

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))

class CircleDomain(GenericDomain):
    """
    ADD DOCUMENTATION
    """

    def __init__(self):
        super(CircleDomain, self).__init__()

        self.fprint("Generating Circle Domain",special="header")

        ### Initialize values from Options ###
        self.radius   = self.params["domain"]["radius"]
        self.center   = self.params["domain"]["center"]
        self.nt = self.params["domain"]["nt"]
        self.mesh_type = self.params["domain"].get("mesh_type","mshr")
        self.x_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.y_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.dim = 2

        ### Get the initial wind direction ###
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)

        ### Calculating the boundary of the shadow ###
        angles = np.linspace(0,2.0*np.pi,self.nt+1)
        self.boundary_line = (self.radius*np.cos(angles),self.radius*np.sin(angles))

        self.fprint("Radius:        {: .2f}".format(self.radius))
        self.fprint("Center:       ({: .2f}, {: .2f})".format(self.center[0],self.center[1]))
        self.fprint("Meshing Type:  {0}".format(self.mesh_type))

        mesh_start = time.time()
        self.fprint("")
        if self.mesh_type == "mshr":

            self.res = self.params["domain"]["res"]
            self.fprint("Generating Mesh Using mshr")

            ### Create Mesh ###
            mshr_circle = Circle(Point(self.center[0],self.center[1]), self.radius, self.nt)
            self.mesh = generate_mesh(mshr_circle,self.res)

        else:
            self.fprint("Generating Rectangle Mesh")

            self.nxy = int(self.nt/4.0)

            ### Create mesh ###
            start = Point(-1.0, -1.0)
            stop  = Point( 1.0,  1.0)
            self.mesh = RectangleMesh(start, stop, int(self.nxy/2.0), int(self.nxy/2.0))
            self.mesh = refine(self.mesh)
            x = self.mesh.coordinates()[:,0]
            y = self.mesh.coordinates()[:,1]
            
            self.fprint("Morphing Mesh")
            if self.mesh_type == "elliptic":
                x_hat, y_hat, z_hat = Elliptical_Grid(x, y, 0, self.radius)
            elif self.mesh_type == "squircular":
                x_hat, y_hat, z_hat = FG_Squircular(x, y, 0, self.radius)
            elif self.mesh_type == "stretch":
                x_hat, y_hat, z_hat = Simple_Stretching(x, y, 0, self.radius)
            else:
                raise ValueError("Mesh type: "+self.mesh_type+" not recognized")

            x_hat += self.center[0]
            y_hat += self.center[1]

            xy_hat_coor = np.array([x_hat, y_hat]).transpose()
            self.mesh.coordinates()[:] = xy_hat_coor
            self.mesh.bounding_box_tree().build(self.mesh)

        ### Create the boundary mesh ###
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.fprint("Mesh Generated: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Define Plane Normal ###
        nom_x = np.cos(self.init_wind)
        nom_y = np.sin(self.init_wind)

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y)
        inflow  = CompiledSubDomain("nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, c0=self.center[0], c1=self.center[1])
        self.boundary_subdomains = [outflow,inflow]
        self.boundary_names = {"inflow":2,"outflow":1}
        self.boundary_types = {"inflow":  ["inflow"],
                               "no_stress": ["outflow"]}

        ### Generate the boundary markers for boundary conditions ###
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            self.boundary_subdomains[i].mark(self.boundary_markers, i+1,check_midpoint=False)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def ground_function(self,x,y,dx=0,dy=0):
        return 0.0

    def RecomputeBoundaryMarkers(self,theta):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Define Plane Normal ###
        nom_x = np.cos(theta)
        nom_y = np.sin(theta)

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y)
        inflow  = CompiledSubDomain("nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, c0=self.center[0], c1=self.center[1])
        self.boundary_subdomains = [outflow,inflow]

        ### Generate the boundary markers for boundary conditions ###
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            self.boundary_subdomains[i].mark(self.boundary_markers, i+1,check_midpoint=False)

        # ### Define center ###
        # c0 = self.center[0]
        # c1 = self.center[1]

        # ### Set Tol ###
        # tol = 1e-5

        # wall_facets = self.boundary_markers.where_equal(self.boundary_names["inflow"]) \
        #             + self.boundary_markers.where_equal(self.boundary_names["outflow"])

        # boundary_val_temp = self.boundary_markers.array()
        # self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        # self.boundary_markers.set_values(boundary_val_temp)

        # for facet_id in wall_facets:
        #     facet = Facet(self.mesh,facet_id)
        #     vert_ids = facet.entities(0)
        #     vert_coords = self.mesh.coordinates()[vert_ids]
        #     x = vert_coords[:,0]
        #     y = vert_coords[:,1]

        #     if all(nom_x*(x-c0)+nom_y*(y-c1)<=0+tol):
        #         self.boundary_markers.set_value(facet_id,self.boundary_names["inflow"])
        #     else:
        #         self.boundary_markers.set_value(facet_id,self.boundary_names["outflow"])

        # mark_stop = time.time()
        # self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))

class RectangleDomain(GenericDomain):
    """
    A rectangle domain is simply a 2D rectangle. This mesh is defined
    by 4 parameters in the param.yaml file. 

    Example:
        In the yaml file define::

            domain: 
                #                      # Description           | Units
                x_range: [-2500, 2500] # x-range of the domain | m
                y_range: [-2500, 2500] # y-range of the domain | m
                nx: 10                 # Number of x-nodes     | -
                ny: 10                 # Number of y-nodes     | -

        This will produce a rectangle with corner points (-2500,-2500) 
        to (2500,2500). The mesh will have *nx* nodes in the *x*-direction,
        and *ny* in the *y*-direction.

    Todo:
        Properly implement a RectangleDomain and 2D in general.
    """

    def __init__(self):
        super(RectangleDomain, self).__init__()

        self.fprint("Generating Rectangle Domain",special="header")

        ### Initialize values from Options ###
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)
        self.x_range = self.params["domain"]["x_range"]
        self.y_range = self.params["domain"]["y_range"]
        self.nx = self.params["domain"]["nx"]
        self.ny = self.params["domain"]["ny"]
        self.fprint("X Range: [{: .2f}, {: .2f}]".format(self.x_range[0],self.x_range[1]))
        self.fprint("Y Range: [{: .2f}, {: .2f}]".format(self.y_range[0],self.y_range[1]))
        self.dim = 2

        ### Create mesh ###
        mesh_start = time.time()
        self.fprint("")
        self.fprint("Generating Mesh")
        start = Point(self.x_range[0], self.y_range[0])
        stop  = Point(self.x_range[1], self.y_range[1])
        self.mesh = RectangleMesh(start, stop, self.nx, self.ny)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.fprint("Mesh Generated: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        front   = CompiledSubDomain("near(x[0], x0) && on_boundary",x0 = self.x_range[0])
        back    = CompiledSubDomain("near(x[0], x1) && on_boundary",x1 = self.x_range[1])
        left    = CompiledSubDomain("near(x[1], y0) && on_boundary",y0 = self.y_range[0])
        right   = CompiledSubDomain("near(x[1], y1) && on_boundary",y1 = self.y_range[1])
        self.boundary_subdomains = [front,back,left,right]
        self.boundary_names = {"front":1,"back":2,"left":3,"right":4}
        self.boundary_types = {"inflow":    ["front","left","right"],
                               "no_stress": ["back"]}

        ### Generate the boundary markers for boundary conditions ###
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            self.boundary_subdomains[i].mark(self.boundary_markers, i+1)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def ground_function(self,x,y,dx=0,dy=0):
        return 0.0

    def RecomputeBoundaryMarkers(self,theta):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        tol = 1e-3

        ### This function rounds to the nearest ordinal direction ###
        theta45 = round((theta-pi/4.0)/(pi/2))*(pi/2)+pi/4.0

        ### Check if the wind angle is a ordinal direction ###
        if   near(theta45, 1.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","left"]
            self.boundary_types["no_stress"] = ["back","right"]
        elif near(theta45, 3.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","left"]
            self.boundary_types["no_stress"] = ["front","right"]
        elif near(theta45, 5.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","right"]
            self.boundary_types["no_stress"] = ["front","left"]
        elif near(theta45, 7.0*pi/4.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","right"]
            self.boundary_types["no_stress"] = ["back","left"]

        ### Check the special cases that the wind angle is a cardinal direction ###
        if   near(theta, 0.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","right","left"]
            self.boundary_types["no_stress"] = ["back"]
        elif near(theta, 1.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","back","left"]
            self.boundary_types["no_stress"] = ["right"]
        elif near(theta, 2.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["back","right","left"]
            self.boundary_types["no_stress"] = ["front"]
        elif near(theta, 3.0*pi/2.0, eps=tol):
            self.boundary_types["inflow"]    = ["front","back","right"]
            self.boundary_types["no_stress"] = ["left"]

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))

class ImportedDomain(GenericDomain):
    """
    This class generates a domain from imported files. This mesh is defined
    by 2 parameters in the param.yaml file. 

    Example:
        In the yaml file define::

            domain: 
                path: "Mesh_data/"
                filetype: "xml.gz"

        The supported filetypes are "xml.gz" and "h5". For "xml.gz" 3 files are
        required: 

            * mesh.xml.gz - this contains the mesh in a format dolfin can handle
            * boundaries.xml.gz - this contains the facet markers that define where the boundaries are
            * topology.txt - this contains the data for the ground topology. 
                It assumes that the coordinates are from a uniform mesh.
                It contains three column: x, y, z. The x and y columns contain 
                just the unique values. The z column contains the ground values
                for every combination of x and y. The first row must be the number
                of points in the x and y direction. Here is an example for z=x+y/10::

                    3 3 9
                    0 0 0.0
                    1 1 0.1
                    2 2 0.2
                        1.0
                        1.1
                        1.2
                        2.0
                        2.1
                        2.2

    """

    def __init__(self):
        super(ImportedDomain, self).__init__()

        self.fprint("Importing Domain",special="header")

        ### Get the file type for the mesh (h5, xml.gz) ###
        self.filetype = self.params["domain"].get("filetype", "xml.gz")
        self.init_wind = self.params.get("solver",{}).get("init_wind_angle",0.0)

        ### Import data from Options ###
        if "path" in self.params["domain"]:
            self.path = self.params["domain"]["path"]
            self.mesh_path  = self.path + "mesh." + self.filetype
            if self.filetype == "xml.gz":
                self.boundary_path = self.path + "boundaries." + self.filetype
            self.terrain_path  = self.path + "topology.txt"
        else:
            self.mesh_path = self.params["domain"]["mesh_path"]
            if self.filetype == "xml.gz":
                self.boundary_path = self.params["domain"]["bound_path"]
            self.terrain_path  = self.params["domain"]["terrain_path"]

        ### Copy Files to input folder ###
        shutil.copy(self.mesh_path,self.params.folder+"input_files/")
        if self.filetype == "xml.gz":
            shutil.copy(self.boundary_path,self.params.folder+"input_files/")

        ### Create the mesh ###
        mesh_start = time.time()
        self.fprint("")
        self.fprint("Importing Mesh")
        if self.filetype == "h5":
            self.mesh = Mesh()
            hdf5 = HDF5File(self.mesh.mpi_comm(), self.mesh_path, 'r')
            hdf5.read(self.mesh, '/mesh', False)
        elif self.filetype == "xml.gz":
            self.mesh = Mesh(self.mesh_path)
        else:
            raise ValueError("Supported mesh types: h5, xml.gz.")
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.dim = self.mesh.topology().dim()
        self.fprint("Mesh Imported: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Calculate the range of the domain and push to options ###
        self.x_range = [min(self.mesh.coordinates()[:,0]),max(self.mesh.coordinates()[:,0])]
        self.y_range = [min(self.mesh.coordinates()[:,1]),max(self.mesh.coordinates()[:,1])]
        self.z_range = [min(self.mesh.coordinates()[:,2]),max(self.mesh.coordinates()[:,2])]
        self.params["domain"]["x_range"] = self.x_range
        self.params["domain"]["y_range"] = self.y_range
        self.params["domain"]["z_range"] = self.z_range

        ### Load the boundary markers ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Importing Boundary Markers")
        if self.filetype == "h5":
            self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.geometry().dim()-1)
            hdf5.read(self.boundary_markers, "/boundaries")
        elif self.filetype == "xml.gz":
            self.boundary_markers = MeshFunction("size_t", self.mesh, self.boundary_path)
        print("Markers Imported")
        self.boundary_names = {"top":1,"bottom":2,"front":3,"back":4,"left":5,"right":6}
        self.boundary_types = {"inflow":          ["front","left","right"],
                               "no_slip":         ["bottom"],
                               "horizontal_slip": ["top"],
                               "no_stress":       ["back"]}
        mark_stop = time.time()
        self.fprint("Boundary Markers Imported: {:1.2f} s".format(mark_stop-mark_start))

        ### Create the interpolation function for the ground ###
        interp_start = time.time()
        self.fprint("")
        self.fprint("Building Interpolating Function")

        self.SetupInterpolatedGround()
        self.ground_function = self.InterplatedGroundFunction

        interp_stop = time.time()
        self.fprint("Interpolating Function Built: {:1.2f} s".format(interp_stop-interp_start))
        self.fprint("Initial Domain Setup",special="footer")

class InterpolatedCylinderDomain(CylinderDomain):
    def __init__(self):
        super(InterpolatedCylinderDomain, self).__init__()
        # self.original_refine = super(InterpolatedCylinderDomain, self).Refine
        # self.original_move = super(InterpolatedCylinderDomain, self).Move

        # if shape == "cylinder":
        #     base = CylinderDomain()
        # elif shape == "box":
        #     base = BoxDomain()
        # else:
        #     raise ValueError(shape+" is not a recognized base type")

        # print(dir(base))
        # self.__dict__.update(base.__dict__)

        interp_start = time.time()
        self.fprint("Building Ground Function",special="header")

        self.analytic = self.params["domain"].get("analytic",False)
        if self.analytic:
            self.SetupAnalyticGround()
            self.ground_function = self.GaussianGroundFuncion
        else:
            self.SetupInterpolatedGround()
            self.ground_function = self.InterplatedGroundFunction

        interp_stop = time.time()
        self.fprint("Interpolating Function Built: {:1.2f} s".format(interp_stop-interp_start),special="footer")

    def Finalize(self):
        self.Move(self.ground_function)
        self.finalized = True
        self.fprint("")
        self.fprint("Domain Finalized")

class InterpolatedBoxDomain(BoxDomain):
    def __init__(self):
        super(InterpolatedBoxDomain, self).__init__()
        # self.original_refine = super(InterpolatedCylinderDomain, self).Refine
        # self.original_move = super(InterpolatedCylinderDomain, self).Move

        # if shape == "cylinder":
        #     base = CylinderDomain()
        # elif shape == "box":
        #     base = BoxDomain()
        # else:
        #     raise ValueError(shape+" is not a recognized base type")

        # print(dir(base))
        # self.__dict__.update(base.__dict__)

        interp_start = time.time()
        self.fprint("Building Ground Function",special="header")

        self.analytic = self.params["domain"].get("analytic",False)
        if self.analytic:
            self.SetupAnalyticGround()
            self.ground_function = self.GaussianGroundFuncion
        else:
            self.SetupInterpolatedGround()
            self.ground_function = self.InterplatedGroundFunction

        interp_stop = time.time()
        self.fprint("Ground Function Built: {:1.2f} s".format(interp_stop-interp_start),special="footer")

    def Finalize(self):
        self.Move(self.ground_function)
        self.finalized = True
        self.fprint("")
        self.fprint("Domain Finalized")
