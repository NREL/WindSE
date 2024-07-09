"""
The DomainManager submodule contains the various classes used for
creating different types of domains

"""

import __main__
import os
from pyadjoint.tape import no_annotations

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    from dolfin import *
    import copy
    import time
    import warnings
    import os
    from sys import platform
    import numpy as np
    from scipy.interpolate import interp2d, interp1d,RectBivariateSpline
    import shutil
    from pyadjoint.tape import stop_annotating

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Import some helper functions
    from windse.helper_functions import meteor_to_math

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

    ### This parameter allows for refining the mesh functions ###
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"

else:
    class SubDomain:
        def __init__(self):
            pass

class SinglePeriodicBoundary(SubDomain):

    def __init__(self,a,b,axis=0):
        super(SinglePeriodicBoundary, self).__init__()
        self.axis   = axis # 0 for periodic x and 1 for periodic y
        self.target = a    # location of the primary boundary
        self.offset = b-a  # the offset from the primary to the secondary

    # returns true primary boundary
    def inside(self, x, on_boundary):
        return near(x[self.axis],self.target) and on_boundary

    # this maps the secondary boundary to the primary boundary
    def map(self, x, y):
        y[:] = x[:] # Set all equal
        y[self.axis] = x[self.axis] - self.offset # Modify the value that is periodic

class DoublePeriodicBoundary(SubDomain):

    def __init__(self,x_range,y_range):
        super(DoublePeriodicBoundary, self).__init__()
        self.x_range = x_range    # location of the primary boundary
        self.y_range = y_range    # location of the primary boundary
        self.x_offset = x_range[1]-x_range[0]  # the offset from the primary to the secondary
        self.y_offset = y_range[1]-y_range[0]  # the offset from the primary to the secondary
        self.eps = 1e-4

    # returns true primary boundary
    def inside(self, x, on_boundary):
        on_either  = near(x[0], self.x_range[0],self.eps)  or near(x[1], self.y_range[0],self.eps)
        on_corner1 = near(x[0], self.x_range[1],self.eps) and near(x[1], self.y_range[0],self.eps)
        on_corner2 = near(x[0], self.x_range[0],self.eps) and near(x[1], self.y_range[1],self.eps)
        return (on_either and not (on_corner1 or on_corner2)) and on_boundary

    # this maps the secondary boundary to the primary boundary
    def map(self, x, y):
        if near(x[0], self.x_range[1],self.eps) and near(x[1], self.y_range[1],self.eps):
            y[0] = x[0] - self.x_offset
            y[1] = x[1] - self.y_offset
            y[2] = x[2]
            # print(1)
        elif near(x[0], self.x_range[1],self.eps):
            y[0] = x[0] - self.x_offset
            y[1] = x[1]
            y[2] = x[2]
            # print(2)
        elif near(x[1], self.y_range[1],self.eps):
            y[0] = x[0]
            y[1] = x[1] - self.y_offset
            y[2] = x[2]
            # print(3)
        else:
            y[0] = -1000
            y[1] = -1000
            y[2] = -1000



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

def DefaultFinalize(dom):
    # self.ComputeCellRadius()
    dom.SetupPeriodicMapping()
    dom.ComputeGlobalHmin()
    dom.DebugOutput()
    dom.finalized = True
    dom.fprint("")
    dom.fprint("Domain Finalized")

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
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        ### Update attributes based on params file ###
        for key, value in self.params["domain"].items():
            if isinstance(value,list):
                setattr(self,key,np.array(value))
            else:
                setattr(self,key,value)

        ### Define the scaling if needed ###
        if self.scaled:
            self.xscale = 1.0e-3
        else:
            self.xscale = 1.0
        self.ground_reference = self.ground_reference*self.xscale

        ### Get the initial wind direction ###
        self.raw_inflow_angle = self.params["boundary_conditions"]["inflow_angle"]
        if self.raw_inflow_angle is None:
            inflow_angle = 270.0 # default inflow angle is from west to east
        elif isinstance(self.raw_inflow_angle,list):
            inflow_angle = self.raw_inflow_angle[0]
        else:
            inflow_angle = self.raw_inflow_angle
        print(inflow_angle)
        self.inflow_angle = meteor_to_math(inflow_angle)
        print(self.inflow_angle)
        self.initial_inflow_angle = self.inflow_angle

        ### Add a space to store periodic boundary conditions ###
        self.periodic_mapping = None
        self.bcs_to_remove = []

    @no_annotations
    def DebugOutput(self):
        if self.debug_mode:
            self.tag_output("dim", self.dim)
            self.tag_output("x_min", self.x_range[0])
            self.tag_output("x_max", self.x_range[1])
            self.tag_output("y_min", self.y_range[0])
            self.tag_output("y_max", self.y_range[1])
            if self.dim == 3:
                self.tag_output("z_min", self.z_range[0])
                self.tag_output("z_max", self.z_range[1])

            # Create a unit vector for integrating
            V = FiniteElement('Lagrange', self.mesh.ufl_cell(), 1)
            V = FunctionSpace(self.mesh,V)
            u = Function(V)
            u.vector()[:] = 1.0

            # Volume Check
            self.volume = assemble(u*dx)
            self.tag_output("volume", self.volume)

            # Check surface area
            ds = Measure("ds", subdomain_data=self.boundary_markers)
            for key, value in self.boundary_names.items():
                if value is not None:
                    self.tag_output("ID_"+key, value)
                    sa = assemble(u*ds(value))
                    self.tag_output("SA_"+key, sa)

    def Plot(self):
        """
        This function plots the domain using matplotlib and saves the
        output to output/.../plots/mesh.pdf
        """

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/mesh.pdf"

        ### Check if folder exists ###
        # if not os.path.exists(folder_string): os.makedirs(folder_string)
        if not os.path.exists(folder_string) and self.params.rank == 0: os.makedirs(folder_string)


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

        self.mesh.coordinates()[:]=self.mesh.coordinates()[:]/self.xscale

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
        self.mesh.coordinates()[:]=self.mesh.coordinates()[:]*self.xscale

    def ExportMesh(self):
        folder_string = self.params.folder+"mesh/exported_mesh/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)
        if self.filetype == "xml.gz":
            File(folder_string+"mesh.xml.gz") << self.mesh
            File(folder_string+"boundaries.xml.gz") << self.boundary_markers
        elif self.filetype == "h5":
            hdfile = HDF5File(self.mesh.mpi_comm(), folder_string+"mesh_data.h5", "w")
            hdfile.write(self.mesh,"mesh")
            hdfile.write(self.boundary_markers,"boundaries")
            hdfile.close()
        if hasattr(self,"terrain"):
            np.savetxt(folder_string+"terrain.txt",self.terrain)

        ### TODO: Export boundary names and ID ###

    def BuildBoundaryMarkers(self):
        self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.boundary_markers.set_all(0)
        for i in range(len(self.boundary_subdomains)):
            if self.boundary_subdomains[i] is not None:
                self.boundary_subdomains[i].mark(self.boundary_markers, i+1,check_midpoint=False)

    def BoxRefine(self,x_range,y_range,z_range=[],expand_factor=1):
        refine_start = time.time()

        ### Calculate Expanded Region ###
        x0, x1 = x_range
        y0, y1 = y_range
        ex = (expand_factor-1)*(x1-x0)/2.0
        ey = (expand_factor-1)*(y1-y0)/2.0
        en = min(ex,ey)
        x0, x1 = x0-en, x1+en
        y0, y1 = y0-en, y1+en
        if self.dim == 3:
            z0, z1 = z_range
            ez = (expand_factor-1)*(z1-z0)/2.0
            z0, z1 = z0-ez, z1+ez

        ### Print Region
        self.fprint("Starting Box Refinement",special="header")
        self.fprint("X Range: [{: .2f}, {: .2f}]".format(x0/self.xscale,x1/self.xscale))
        self.fprint("Y Range: [{: .2f}, {: .2f}]".format(y0/self.xscale,y1/self.xscale))
        if self.dim == 3:
            self.fprint("Z Range: [{: .2f}, {: .2f}]".format(z0/self.xscale,z1/self.xscale))

        ### Create cell markers ###
        self.fprint("Marking Cells")
        cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
        for cell in cells(self.mesh):
            in_square = between(cell.midpoint()[0],(x0,x1)) and between(cell.midpoint()[1],(y0,y1))
            if self.dim == 3:
                # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                in_z = between(cell.midpoint()[2],(z0,z1+g_z))
                in_region = in_square and in_z
            else:
                in_region = in_square
            if in_region:
                cell_f[cell] = True

        ### Refine Mesh
        self.Refine(cell_f)

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def CylinderRefine(self,center,radius,height=0,expand_factor=1):

        ### Calculate Expanded Region ###
        radius = expand_factor*radius
        if self.dim == 3:
            ez = (expand_factor-1)*height/2.0
            center[2] -= ez
            height = expand_factor*height

        ### Print Region Data ###
        refine_start = time.time()
        self.fprint("Starting Vertical Cylinder Refinement",special="header")
        self.fprint("Cylinder Radius: {:.2f}".format(radius/self.xscale))
        if self.dim == 3:
            self.fprint("Cylinder Height: {:.2f}".format(height/self.xscale))
            self.fprint("Cylinder Center: ({:.2f}, {:.2f}, {:.2f})".format(center[0]/self.xscale,center[1]/self.xscale,center[2]/self.xscale))
        else:
            self.fprint("Cylinder Center: ({:.2f}, {:.2f})".format(center[0]/self.xscale,center[1]/self.xscale))

        cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
        cells_marked = 0

        ### Check if we are refining in a circle ###
        for cell in cells(self.mesh):
            in_circle = (cell.midpoint()[0]-center[0])**2.0+(cell.midpoint()[1]-center[1])**2.0<=radius**2.0
            if self.dim == 3:
                # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
                in_z = between(cell.midpoint()[2],(center[2],center[2]+height+g_z))
                in_region = in_circle and in_z
            else:
                in_region = in_circle
            if in_region:
                cell_f[cell] = True

        ### Refine Mesh
        self.Refine(cell_f)

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def StreamRefine(self,center,radius,length,theta=0,pivot_offset=0,expand_factor=1):
        refine_start = time.time()

        ### Calculate Expanded Region ###
        en = 2*(expand_factor-1)*radius/2.0
        center[0] -= en*np.cos(theta)
        center[1] -= -en*np.sin(theta)
        length = length+2*en
        radius = expand_factor*radius

        ### Output data
        self.fprint("Starting Vertical Cylinder Refinement",special="header")
        self.fprint("Cylinder Radius: {:.2f}".format(radius/self.xscale))
        self.fprint("Cylinder Length: {:.2f}".format(length/self.xscale))
        if self.dim == 3:
            self.fprint("Cylinder Center: ({:.2f}, {:.2f}, {:.2f})".format(center[0]/self.xscale,center[1]/self.xscale,center[2]/self.xscale))
        else:
            self.fprint("Cylinder Center: ({:.2f}, {:.2f})".format(center[0]/self.xscale,center[1]/self.xscale))
        self.fprint("Cylinder Rotation: {:.2f}".format(theta))

        cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
        cells_marked = 0

        ### Check if we are refining in a circle ###
        for cell in cells(self.mesh):

            xs = cell.midpoint()[0]
            ys = cell.midpoint()[1]
            x = np.cos(theta)*(xs-center[0]-pivot_offset)-np.sin(theta)*(ys-center[1])+ center[0]+pivot_offset
            y = np.sin(theta)*(xs-center[0]-pivot_offset)+np.cos(theta)*(ys-center[1])+ center[1]

            in_stream = between(x,(center[0],center[0]+length))
            if self.dim == 3:
                z = cell.midpoint()[2] + self.ground_function(x,y)
                in_circle = (y-center[1])**2.0+(z-center[2])**2.0<=radius**2.0
            else:
                in_circle = between(y,(center[1]-radius,center[1]+radius))
            in_region = in_circle and in_stream

            if in_region:
                cell_f[cell] = True

        ### Refine Mesh
        self.Refine(cell_f)

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def Refine(self, cellmarkers=None):

        # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # print(parameters["refinement_algorithm"])
        # print(cellmarkers)
        # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        old_verts = self.mesh.num_vertices()
        old_cells = self.mesh.num_cells()
        self.fprint("Refining Mesh")
        self.mesh = refine(self.mesh,cellmarkers)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")

        if self.params.num_procs == 1:
            self.boundary_markers = adapt(self.boundary_markers,self.mesh)
        elif self.params.num_procs > 1:
            # This isn't needed for actual boundary marking, but it helps it pass a test later on
            self.BuildBoundaryMarkers()

        self.my_fs = FunctionSpace
        self.my_fn = Function
        # print('After refinement')
        # Q = FunctionSpace(self.mesh, 'P', 1)
        # f = Function(Q)
        # print('finished after refinement')

        self.fprint("Original Mesh Vertices: {:d}".format(old_verts))
        self.fprint("Original Mesh Cells:    {:d}".format(old_cells))
        self.fprint("New Mesh Vertices:      {:d}".format(self.mesh.num_vertices()))
        self.fprint("New Mesh Cells:         {:d}".format(self.mesh.num_cells()))



















    # def Refine(self,num,region=None,region_type=None,cell_markers=None):
    #     """
    #     This function can be used to refine the mesh. If a region is
    #     specified, the refinement is local

    #     Args:
    #         num (int): the number of times to refine

    #     :Keyword Arguments:
    #         * **region** (*list*):
    #                             for square region use: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
    #                             for circle region use: [[radius],[c_x,c_y],[zmin,zmax]]
    #         * **region_type** (*str*): Either "circle" or "square"
    #         * **cell_markers** (*:meth:dolfin.mesh.MeshFunction*): A cell function marking which cells to refine
    #     """
    #     refine_start = time.time()

    #     ### Print some useful stats
    #     self.fprint("Starting Mesh Refinement",special="header")
    #     if cell_markers is not None:
    #         self.fprint("Region Type: {0}".format("cell_markers"))
    #     elif region is not None:
    #         self.fprint("Region Type: {0}".format(region_type))
    #         if "circle" in region_type:
    #             self.fprint("Circle Radius: {:.2f}".format(region[0][0]/self.xscale))
    #             self.fprint("Circle Center: ({:.2f}, {:.2f})".format(region[1][0]/self.xscale,region[1][1]/self.xscale))
    #             if self.dim == 3:
    #                 self.fprint("Z Range:       [{:.2f}, {:.2f}]".format(region[2][0]/self.xscale,region[2][1]/self.xscale))
    #         else:
    #             self.fprint("X Range: [{: .2f}, {: .2f}]".format(region[0][0]/self.xscale,region[0][1]/self.xscale))
    #             self.fprint("Y Range: [{: .2f}, {: .2f}]".format(region[1][0]/self.xscale,region[1][1]/self.xscale))
    #             if self.dim == 3:
    #                 self.fprint("Z Range: [{: .2f}, {: .2f}]".format(region[2][0]/self.xscale,region[2][1]/self.xscale))
    #     else:
    #         self.fprint("Region Type: {0}".format("full"))

    #     ### Mark cells for refinement
    #     for i in range(num):
    #         if num>1:
    #             step_start = time.time()
    #             self.fprint("Refining Mesh Step {:d} of {:d}".format(i+1,num), special="header")
    #         else:
    #             self.fprint("")


    #         ### Check if cell markers were provided ###
    #         if cell_markers is not None:
    #             cell_f = cell_markers
    #             self.fprint("Cells Marked for Refinement: {:d}".format(sum(cell_markers.array())))

    #         ### Check if a region was provided ###
    #         elif region is not None:
    #             self.fprint("Marking Cells")

    #             ### Create an empty cell marker function
    #             cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
    #             cells_marked = 0

    #             ### Check if we are refining in a circle ###
    #             if "circle" in region_type:
    #                 radius=region[0][0]
    #                 for cell in cells(self.mesh):
    #                     in_circle = (cell.midpoint()[0]-region[1][0])**2.0+(cell.midpoint()[1]-region[1][1])**2.0<=radius**2.0
    #                     if self.dim == 3:
    #                         # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
    #                         g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
    #                         in_z = between(cell.midpoint()[2],(region[2][0],region[2][1]+g_z))
    #                         in_region = in_circle and in_z
    #                     else:
    #                         in_region = in_circle
    #                     if in_region:
    #                         cell_f[cell] = True
    #                         cells_marked += 1

    #             ### or a rectangle ###
    #             else:
    #                 cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),False)
    #                 for cell in cells(self.mesh):
    #                     in_square = between(cell.midpoint()[0],tuple(region[0])) and between(cell.midpoint()[1],tuple(region[1]))
    #                     if self.dim == 3:
    #                         # g_z = 0.0#self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
    #                         g_z = self.ground_function(cell.midpoint()[0],cell.midpoint()[1])
    #                         in_z = between(cell.midpoint()[2],(region[2][0],region[2][1]+g_z))
    #                         in_region = in_square and in_z
    #                     else:
    #                         in_region = in_square
    #                     if in_region:
    #                         cell_f[cell] = True
    #                         cells_marked += 1
    #             self.fprint("Cells Marked for Refinement: {:d}".format(cells_marked))


    #         ### If neither a region or cell markers were provided, Refine everwhere ###
    #         else:
    #             cell_f = MeshFunction('bool', self.mesh, self.mesh.geometry().dim(),True)

    #         old_verts = self.mesh.num_vertices()
    #         old_cells = self.mesh.num_cells()
    #         self.fprint("Refining Mesh")
    #         self.mesh = refine(self.mesh,cell_f)
    #         self.bmesh = BoundaryMesh(self.mesh,"exterior")
    #         self.boundary_markers = adapt(self.boundary_markers,self.mesh)

    #         self.fprint("Original Mesh Vertices: {:d}".format(old_verts))
    #         self.fprint("Original Mesh Cells:    {:d}".format(old_cells))
    #         self.fprint("New Mesh Vertices:      {:d}".format(self.mesh.num_vertices()))
    #         self.fprint("New Mesh Cells:         {:d}".format(self.mesh.num_cells()))
    #         if num>1:
    #             step_stop = time.time()
    #             self.fprint("Step {:d} of {:d} Finished: {:1.2f} s".format(i+1,num,step_stop-step_start), special="footer")

    #     refine_stop = time.time()
    #     self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")


    def WarpSplit(self,h,s):
        """
        This function warps the mesh to shift more cells towards the ground.
        is achieved by spliting the domain in two and moving the cells so
        that a percentage of them are below the split.

        Args:
            h (float): the height that split occurs
            s (float): the percent below split in the range [0,1)
        """

        if self.mesh.topology().dim() != 3:
            raise ValueError("Cannot warp a 2D mesh")


        warp_start = time.time()
        self.fprint("Starting Mesh Warping",special="header")
        self.fprint("Height of Split:     {:1.2f} m".format(h))
        self.fprint("Percent below Split: {:1.2f}%".format(s*100.0))

        z0 = self.z_range[0]
        z1 = self.z_range[1]

        # h1 = 60
        # h2 = 160

        z = copy.deepcopy(self.mesh.coordinates()[:,2])
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
        self.mesh.coordinates()[:,2]=z_new
        self.mesh.bounding_box_tree().build(self.mesh)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")

        warp_stop = time.time()
        self.fprint("Warping Finished: {:1.2f} s".format(warp_stop-warp_start),special="footer")


    def WarpSmooth(self,s):
        """
        This function warps the mesh to shift more cells towards the ground.
        The cells are shifted based on the function:

        .. math::

            z_new = z_0 + (z_1-z_0) \\left( \\frac{z_old-z_0}{z_1-z_0} \\right)^{s},

        where :math:`z_0` is the ground and :math:`z_1` is the top of the domain.

        Args:
            s (float): compression strength
        """

        if self.mesh.topology().dim() != 3:
            raise ValueError("Cannot warp a 2D mesh")

        warp_start = time.time()
        self.fprint("Starting Mesh Warping",special="header")
        self.fprint("Compression Strength: {:1.4f}".format(s))
        self.fprint("Moving Nodes")

        z=self.mesh.coordinates()[:,2].copy()
        z0 = self.z_range[0]
        z1 = self.z_range[1]
        z1 = z0 + (z1 - z0)*(abs((z-z0)/(z1-z0)))**s
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
        ### The way dolfin_adjoint overwrites ALE.move is destructive so we have to use a work around
        if self.params.dolfin_adjoint:
            ALE.move(self.mesh,self.bmesh)
        else:
            with stop_annotating():
                ALE.move(self.mesh,self.bmesh)

        move_stop = time.time()
        self.fprint("Mesh Moved: {:1.2f} s".format(move_stop-move_start),special="footer")

    def ComputeGlobalHmin(self):
        '''
        This function computes the global hmin() value, since the hmin method
        isn't parallel/global by default. Note that this assumes uniform spacing
        in the x-, y-, and z- directions
        '''

        hmin = self.mesh.hmin()/np.sqrt(self.mesh.topology().dim())
        hmin_global = np.zeros(self.params.num_procs)
        self.params.comm.Allgather(hmin, hmin_global)
        self.global_hmin = np.amin(hmin_global)


    def Finalize(self):
        DefaultFinalize(self)

    # def ComputeCellRadius(self):
    #     self.mesh_radius = MeshFunction("double", self.mesh, self.mesh.topology().dim())
    #     for cell in cells(self.mesh):
    #         self.mesh_radius.set_value(cell.index(),cell.h())
    #         # self.mesh_radius.set_value(cell.index(),cell.inradius())
    #         # self.mesh_radius.set_value(cell.index(),cell.circumradius())

    def SetupInterpolatedGround(self):
        self.fprint("Ground Type: Interpolated From File")

        ### Import data from Options ###
        if self.path is not None and self.terrain_path is None:
            self.terrain_path = self.path + "terrain.txt"

        ### Copy Files to input folder ###
        shutil.copy(self.terrain_path,self.params.folder+"input_files/")

        self.fprint("Path: {0}".format(self.terrain_path),offset=1)

        ### import ground data
        self.terrain = np.loadtxt(self.terrain_path)
        x_data = self.terrain[1:,0]*self.xscale
        y_data = self.terrain[1:,1]*self.xscale
        z_data = self.terrain[1:,2]*self.xscale

        ### generate interpolating function
        x_data = np.sort(np.unique(x_data))
        y_data = np.sort(np.unique(y_data))
        z_data = np.reshape(z_data,(int(self.terrain[0,0]),int(self.terrain[0,1])))
        self.terrain_interpolated = RectBivariateSpline(x_data,y_data,z_data.T)


        def InterplatedGroundFunction(x,y,dx=0,dy=0):
            if dx == 0 and dy == 0:
                return float(self.terrain_interpolated(x,y)[0]+self.ground_reference)
            else:
                return float(self.terrain_interpolated(x,y,dx=dx,dy=dy)[0])

        self.ground_function = InterplatedGroundFunction

    def SetupAnalyticGround(self):
        if self.analytic == "gaussian":

            sigma_x = self.gaussian["sigma_x"]*self.xscale
            sigma_y = self.gaussian["sigma_y"]*self.xscale
            theta = self.gaussian["theta"]
            amp = self.gaussian["amp"]*self.xscale
            center = np.array(self.gaussian["center"])*self.xscale
            x0 = center[0]
            y0 = center[1]
            self.fprint("")
            self.fprint("Ground Type: Gaussian Hill")
            self.fprint("Hill Center:   ({: .2f}, {: .2f})".format(x0,y0),offset=1)
            self.fprint("Hill Rotation:  {: <7.2f}".format(theta),offset=1)
            self.fprint("Hill Amplitude: {: <7.2f}".format(amp),offset=1)
            self.fprint("Hill sigma_x:   {: <7.2f}".format(sigma_x),offset=1)
            self.fprint("Hill sigma_y:   {: <7.2f}".format(sigma_y),offset=1)
            a = np.cos(theta)**2/(2*sigma_x**2) + np.sin(theta)**2/(2*sigma_y**2)
            b = np.sin(2*theta)/(4*sigma_y**2) - np.sin(2*theta)/(4*sigma_x**2)
            c = np.cos(theta)**2/(2*sigma_y**2) + np.sin(theta)**2/(2*sigma_x**2)

            def GaussianGroundFuncion(x,y,dx=0,dy=0):
                return amp*exp( - (a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2)**2)+self.ground_reference

            self.ground_function = GaussianGroundFuncion

        elif self.analytic == "plane":
            x0 = np.array(self.plane["intercept"])*self.xscale
            mx = self.plane["mx"]*self.xscale
            my = self.plane["my"]*self.xscale
            self.fprint("")
            self.fprint("Ground Type: Plane")
            self.fprint("Intercept: ({: .2f}, {: .2f}, {: .2f})".format(*x0),offset=1)
            self.fprint("X Slope:     {: <7.6f}".format(mx),offset=1)
            self.fprint("Y Slope:     {: <7.6f}".format(my),offset=1)

            def PlaneGroundFuncion(x,y,dx=0,dy=0):
                if dx == 1:
                    val = mx
                elif dy == 1:
                    val = my
                elif abs(dx)+abs(dy) >=2:
                    val = 0
                else:
                    val = (mx*(x-x0[0])+my*(y-x0[1]))+x0[2]+self.ground_reference
                return val

            self.ground_function = PlaneGroundFuncion

        else:
            raise ValueError(self.analytic + "is not an implemented type")

    def ground_function(self,x,y,dx=0,dy=0):
        if dx == 0 and dy == 0:
            return self.ground_reference
        else:
            return 0.0

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

    def SetupPeriodicMapping(self):
        # TODO: account for periodic in refine
        if self.type not in ["box","rectangle"] and (self.streamwise_periodic or self.spanwise_periodic):
            raise TypeError("Periodic boundary conditions are only supported for box and rectangle")

        if self.streamwise_periodic and self.spanwise_periodic:
            self.fprint("Setting up dual periodic boundary mapping")
            self.periodic_mapping = DoublePeriodicBoundary(self.x_range, self.y_range)
            self.bcs_to_remove = ["east","west","north","south"]

        elif self.streamwise_periodic:
            self.fprint("Setting up streamwise periodic boundary mapping")
            self.periodic_mapping = SinglePeriodicBoundary(self.x_range[0], self.x_range[1], axis=0)
            self.bcs_to_remove = ["east","west"]

        elif self.spanwise_periodic:
            self.fprint("Setting up spanwise periodic boundary mapping")
            self.periodic_mapping = SinglePeriodicBoundary(self.y_range[0], self.y_range[1], axis=1)
            self.bcs_to_remove = ["north","south"]

        else:
            self.periodic_mapping = None
            self.bcs_to_remove = None

    def RecomputeBoundaryMarkers(self,inflow_angle):
        raise NotImplementedError("This Domain type does not support nonzero inflow angles.")


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
        self.x_range = self.x_range*self.xscale
        self.y_range = self.y_range*self.xscale
        self.z_range = self.z_range*self.xscale
        self.dim = 3


        ### Print Some stats ###
        self.fprint("X Range: [{: .2f}, {: .2f}]".format(self.x_range[0]/self.xscale,self.x_range[1]/self.xscale))
        self.fprint("Y Range: [{: .2f}, {: .2f}]".format(self.y_range[0]/self.xscale,self.y_range[1]/self.xscale))
        self.fprint("Z Range: [{: .2f}, {: .2f}]".format(self.z_range[0]/self.xscale,self.z_range[1]/self.xscale))

        ### Create mesh ###
        mesh_start = time.time()
        self.fprint("")
        self.fprint("Generating Mesh")

        if self.mesh_type == "gmsh":

            if (self.params.rank == 0):

                # only dependencies for the gmsh case
                import gmsh
                import meshio
                import tempfile

                # start up gmsh model
                gmsh.initialize()
                gmsh.model.add("boxmesh")

                # true extents
                Lx_true = self.x_range[1] - self.x_range[0]
                Ly_true = self.y_range[1] - self.y_range[0]
                Lz_true = self.z_range[1] - self.z_range[0]

                # either do stretching to preserve
                stretch_anisotropic = True # TODO: turn into an input param
                if stretch_anisotropic:
                    Lx_prime = 1.0
                    Ly_prime = self.ny/self.nx
                    Lz_prime = self.nz/self.nx
                else:
                    Lx_prime = Lx_true
                    Ly_prime = Ly_true
                    Lz_prime = Lz_true

                # create points that constitute box's lower surface
                pt_bx000 = gmsh.model.geo.addPoint(
                    self.x_range[0]/Lx_true*Lx_prime,
                    self.y_range[0]/Ly_true*Ly_prime,
                    self.z_range[0]/Lz_true*Lz_prime,
                )
                pt_bx100 = gmsh.model.geo.addPoint(
                    self.x_range[1]/Lx_true*Lx_prime,
                    self.y_range[0]/Ly_true*Ly_prime,
                    self.z_range[0]/Lz_true*Lz_prime,
                )
                pt_bx010 = gmsh.model.geo.addPoint(
                    self.x_range[0]/Lx_true*Lx_prime,
                    self.y_range[1]/Ly_true*Ly_prime,
                    self.z_range[0]/Lz_true*Lz_prime,
                )
                pt_bx110 = gmsh.model.geo.addPoint(
                    self.x_range[1]/Lx_true*Lx_prime,
                    self.y_range[1]/Ly_true*Ly_prime,
                    self.z_range[0]/Lz_true*Lz_prime,
                )

                # create lines of box's lower surface
                ln_bx000 = gmsh.model.geo.addLine(pt_bx000, pt_bx100)
                ln_bx001 = gmsh.model.geo.addLine(pt_bx100, pt_bx110)
                ln_bx010 = gmsh.model.geo.addLine(pt_bx010, pt_bx110)
                ln_bx011 = gmsh.model.geo.addLine(pt_bx000, pt_bx010)

                # create curve loop for lower surface
                cl_bx0 = gmsh.model.geo.addCurveLoop(
                    [ln_bx000, ln_bx001, -ln_bx010, -ln_bx011]
                )
                # create plane surface
                ps_bx0 = gmsh.model.geo.addPlaneSurface([cl_bx0])

                bx_simple = gmsh.model.geo.extrude(
                    [(2, ps_bx0)],
                    0.0,
                    0.0,
                    (self.z_range[1] - self.z_range[0])/Lz_true*Lz_prime
                )

                # synchronize the geometry
                gmsh.model.geo.synchronize()

                # characteristic length based on regular tetrahedron volume formula
                lc_mean = (
                    6*np.sqrt(2)
                    *(Lx_prime/self.nx)
                    *(Ly_prime/self.ny)
                    *(Lz_prime/self.nz)
                )**(1/3)
                # set the mesh size to the characteristic length
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_mean)

                # generate and output
                gmsh.model.mesh.generate(3)

                # create a temporary directory to hold mesh IO files for conversion
                dir_for_meshes = tempfile.TemporaryDirectory()

                # write and finalize gmsh
                gmsh.write(os.path.join(dir_for_meshes.name, "dummy.msh"))
                gmsh.finalize() # don't need it anymore!

                # use meshio to convert from gmsh to dolfin xml file
                mesh_meshio = meshio.read(os.path.join(dir_for_meshes.name, "dummy.msh"))
                meshio.write(
                    os.path.join(dir_for_meshes.name, "dummy.xml"),
                    mesh_meshio,
                    file_format="dolfin-xml",
                )

            # broadcast the temporary directory for shared use
            dirname_for_meshes = self.params.comm.bcast(dir_for_meshes.name, root=0)

            # load the xml into windse
            self.mesh = Mesh(os.path.join(dirname_for_meshes, "dummy.xml"))

            # delete the temporary directory
            if self.params.rank == 0:
                dir_for_meshes.cleanup()

            # stretch the coordinates back to the physical dimensions (now with the
            # specified aspect ratios)
            self.mesh.coordinates()[:,0] = self.mesh.coordinates()[:,0]*Lx_true/Lx_prime
            self.mesh.coordinates()[:,1] = self.mesh.coordinates()[:,1]*Ly_true/Ly_prime
            self.mesh.coordinates()[:,2] = self.mesh.coordinates()[:,2]*Lz_true/Lz_prime

        else:
            start = Point(self.x_range[0], self.y_range[0], self.z_range[0])
            stop  = Point(self.x_range[1], self.y_range[1], self.z_range[1])
            self.mesh = BoxMesh(start, stop, self.nx, self.ny, self.nz)

        # box = Box(start,stop)
        # self.mesh = generate_mesh(box,self.nx)
        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.fprint("Mesh Generated: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        east    = CompiledSubDomain("near(x[0], x1, tol) && on_boundary",x1 = self.x_range[1], tol = 1e-10)
        north   = CompiledSubDomain("near(x[1], y1, tol) && on_boundary",y1 = self.y_range[1], tol = 1e-10)
        west    = CompiledSubDomain("near(x[0], x0, tol) && on_boundary",x0 = self.x_range[0], tol = 1e-10)
        south   = CompiledSubDomain("near(x[1], y0, tol) && on_boundary",y0 = self.y_range[0], tol = 1e-10)
        bottom  = CompiledSubDomain("near(x[2], z0, tol) && on_boundary",z0 = self.z_range[0], tol = 1e-10)
        top     = CompiledSubDomain("near(x[2], z1, tol) && on_boundary",z1 = self.z_range[1], tol = 1e-10)
        self.boundary_subdomains = [east,north,west,south,bottom,top]
        self.boundary_names = {"east":1,"north":2,"west":3,"south":4,"bottom":5,"top":6,"inflow":None,"outflow":None}
        self.boundary_types = {"inflow":    ["west","south","north"],
                               "no_slip":   ["bottom"],
                               "free_slip": ["top"],
                               "no_stress": ["east"]}

        self.SetupPeriodicMapping()

        ### Generate the boundary markers for boundary conditions ###
        if self.params.num_procs == 1:
            self.BuildBoundaryMarkers()

            ### Rotate Boundary
            if not near(self.inflow_angle,0.0):
                self.RecomputeBoundaryMarkers(self.inflow_angle)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))


        self.fprint("Initial Domain Setup",special="footer")

    def RecomputeBoundaryMarkers(self,inflow_angle):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Get ids for new markings ###
        inflow_id  = self.boundary_names[self.boundary_types["inflow"][0]]
        outflow_id = self.boundary_names[self.boundary_types["no_stress"][0]]
        east_id    = self.boundary_names["east"]
        north_id   = self.boundary_names["north"]
        west_id    = self.boundary_names["west"]
        south_id   = self.boundary_names["south"]

        ### Set up the baseline (angle=0) order ###
        cardinal_ids = [east_id,north_id,west_id,south_id]
        diagonal_ids = [outflow_id,outflow_id,inflow_id,inflow_id]

        ### Count the number of pi/2 sections in the new inflow_angle ###
        turns = inflow_angle/(pi/2)

        ### New order ###
        tol = 1e-3
        if turns % 1 <= tol: # we are at a cardinal direction
            new_order = np.roll(cardinal_ids,int(turns))
        else:
            new_order = np.roll(diagonal_ids,int(turns))

        ### Get a list of all wall facets ###
        wall_facets = []
        for i in cardinal_ids:
            wall_facets += self.boundary_markers.where_equal(i)

        ### Iterate through facets remarking as prescribed by new_order ###
        for facet_id in wall_facets:

            ### Get the facet normal
            facet = Facet(self.mesh,facet_id)
            facet_normal = facet.normal().array()
            nx = np.sign(facet_normal[0])*(abs(facet_normal[0])>tol)
            ny = np.sign(facet_normal[1])*(abs(facet_normal[1])>tol)

            ### Map the normal to the order east, north, west, south
            face_id = int(abs(nx+2*ny-1))

            ### remark the boundary ###
            self.boundary_markers.set_value(facet_id,new_order[face_id])

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
        self.radius   = self.radius *self.xscale
        self.center   = self.center *self.xscale
        self.z_range  = self.z_range*self.xscale
        self.x_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.y_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.dim = 3

        ### Calculating the boundary of the shadow ###
        angles = np.linspace(0,2.0*np.pi,self.nt+1)
        self.boundary_line = np.array((self.radius*np.cos(angles)+self.center[0],self.radius*np.sin(angles)+self.center[1]))

        self.fprint("Radius:        {: .2f}".format(self.radius/self.xscale))
        self.fprint("Center:       ({: .2f}, {: .2f})".format(self.center[0]/self.xscale,self.center[1]/self.xscale))
        self.fprint("Z Range:      [{: .2f}, {: .2f}]".format(self.z_range[0]/self.xscale,self.z_range[1]/self.xscale))
        self.fprint("Meshing Type:  {0}".format(self.mesh_type))

        mesh_start = time.time()
        self.fprint("")
        if self.mesh_type == "mshr":
            raise NotImplementedError("Mshr is no longer supported, use gmsh instead.")
            # self.fprint("Generating Mesh Using mshr")

            # ### Create Mesh ###
            # # mshr_circle = Circle(Point(self.center[0],self.center[1]), self.radius, self.nt)
            # # mshr_domain = Extrude2D(mshr_circle,self.z_range[1]-self.z_range[0])
            # top    = Point(self.center[0],self.center[1],self.z_range[1])
            # bottom = Point(self.center[0],self.center[1],self.z_range[0])
            # mshr_domain = Cylinder(top,bottom,self.radius,self.radius,self.nt)
            # self.mesh = generate_mesh(mshr_domain,self.res)
            # # self.mesh = refine(self.mesh)
            # # self.mesh = refine(self.mesh)
            # # self.mesh = refine(self.mesh)

            # # z = self.mesh.coordinates()[:,2]#+self.z_range[0]
            # # self.mesh.coordinates()[:,2] = z
            # # self.mesh.bounding_box_tree().build(self.mesh)

        elif self.mesh_type == "gmsh":
            self.fprint("Generating Mesh Using gmsh")

            if (self.params.rank == 0):

                # only dependencies for the gmsh case
                import gmsh
                import meshio
                import tempfile

                # start up gmsh model
                gmsh.initialize()
                gmsh.model.add("circmesh")
                factory = gmsh.model.occ

                # true extents
                R_true = self.radius
                Lz_true = self.z_range[1] - self.z_range[0]

                stretch_anisotropic = True # TODO: turn into an input param
                if stretch_anisotropic:
                    R_prime = 1.0
                    Lz_prime = 2*np.pi*self.nz/self.nt

                    lc_mean = Lz_prime/self.nz
                else:
                    raise NotImplementedError("logic for ignoring scaling not implemented. -cfrontin")

                # set up the geometry
                x_circ = 0.0
                y_circ = 0.0
                z_circ = 0.0 # self.z_range[0]

                # create a circle in the horizon plane
                ln_circ = factory.addCircle(x_circ, y_circ, z_circ, R_prime)
                cl_circ = factory.addCurveLoop([ln_circ])
                ps_circ = factory.addPlaneSurface([cl_circ])

                geom_cyl = factory.extrude([(2, ps_circ)], 0.0, 0.0, Lz_prime)

                # synchronize the geometry
                factory.synchronize()

                # set the mesh size
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_mean)

                # generate and output
                gmsh.model.mesh.generate(3)
                # gmsh.write("dummy.msh") # DEBUG!!!!!

                # create a temporary directory to hold mesh IO files for conversion
                dir_for_meshes = tempfile.TemporaryDirectory()

                # write and finalize gmsh
                gmsh.write(os.path.join(dir_for_meshes.name, "dummy.msh"))
                gmsh.finalize() # don't need it anymore!

                # use meshio to convert from gmsh to dolfin xml file
                mesh_meshio = meshio.read(os.path.join(dir_for_meshes.name, "dummy.msh"))

                meshio.write(
                    os.path.join(dir_for_meshes.name, "dummy.xml"),
                    mesh_meshio,
                    file_format="dolfin-xml",
                )

            # broadcast the temporary directory for shared use
            dirname_for_meshes = self.params.comm.bcast(dir_for_meshes.name, root=0)

            # load the xml into windse
            self.mesh = Mesh(os.path.join(dirname_for_meshes, "dummy.xml"))

            # delete the temporary directory
            if self.params.rank == 0:
                dir_for_meshes.cleanup()

            # stretch the coordinates back to the physical dimensions (now with the
            # specified aspect ratios)
            self.mesh.coordinates()[:,0] = self.mesh.coordinates()[:,0]*R_true/R_prime + self.center[0]
            self.mesh.coordinates()[:,1] = self.mesh.coordinates()[:,1]*R_true/R_prime + self.center[1]
            self.mesh.coordinates()[:,2] = self.mesh.coordinates()[:,2]*Lz_true/Lz_prime + self.z_range[0]

        else:
            self.fprint("Generating Box Mesh")

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
        nom_x = np.cos(self.inflow_angle)
        nom_y = np.sin(self.inflow_angle)

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        outflow = CompiledSubDomain(
            "on_boundary",
            nx=nom_x, ny=nom_y,
            z0=self.z_range[0], z1=self.z_range[1],
        )
        inflow  = CompiledSubDomain(
            "nx*(x[0]-c0)+ny*(x[1]-c1)<=0 && on_boundary",
            nx=nom_x, ny=nom_y,
            z0=self.z_range[0], z1=self.z_range[1],
            c0=self.center[0], c1=self.center[1],
        )
        top     = CompiledSubDomain(
            "near(x[2], z1, tol) && on_boundary",
            z1=self.z_range[1], tol=1e-10
        )
        bottom  = CompiledSubDomain(
            "near(x[2], z0, tol) && on_boundary",
            z0=self.z_range[0], tol=1e-10,
        )
        # outflow = CompiledSubDomain("x[2] > z0 && x[2] < z1 && on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1])
        # inflow  = CompiledSubDomain("x[2] > z0 && x[2] < z1 && nx*(x[0]-c0)+ny*(x[1]-c1)<=0 && on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1], c0=self.center[0], c1=self.center[1])
        # top     = CompiledSubDomain("near(x[2], z1, tol) && on_boundary",z1 = self.z_range[1],tol = 1e-10)
        # bottom  = CompiledSubDomain("near(x[2], z0, tol) && on_boundary",z0 = self.z_range[0],tol = 1e-10)
        self.boundary_subdomains = [None,None,None,None,outflow,inflow,bottom,top]
        self.boundary_names = {"west":None,"east":None,"south":None,"north":None,"bottom":7,"top":8,"inflow":6,"outflow":5}
        self.boundary_types = {"inflow":          ["inflow"],
                               "no_stress":       ["outflow"],
                               "free_slip":       ["top"],
                               "no_slip":         ["bottom"]}

        ### Generate the boundary markers for boundary conditions ###
        self.BuildBoundaryMarkers()

        ### Rotate Boundary
        if not near(self.inflow_angle,0.0):
            self.RecomputeBoundaryMarkers(self.inflow_angle)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def RecomputeBoundaryMarkers(self,inflow_angle):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Define Plane Normal ###
        nom_x = np.cos(inflow_angle)
        nom_y = np.sin(inflow_angle)

        # ### Define Boundary Subdomains ###
        # outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1])
        # inflow  = CompiledSubDomain("x[2] >= z0 && x[2] <= z1 && nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, z0 = self.z_range[0], z1 = self.z_range[1], c0=self.center[0], c1=self.center[1])
        # top     = CompiledSubDomain("near(x[2], z1) && on_boundary",z1 = self.z_range[1])
        # bottom  = CompiledSubDomain("near(x[2], z0) && on_boundary",z0 = self.z_range[0])
        # self.boundary_subdomains = [None,None,None,None,bottom,top,inflow,outflow]

        # ### Generate the boundary markers for boundary conditions ###
        # self.BuildBoundaryMarkers()

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
        self.radius   = self.radius*self.xscale
        self.center   = self.center*self.xscale
        self.x_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.y_range  = [self.center[0]-self.radius,self.center[1]+self.radius]
        self.dim = 2

        ### Calculating the boundary of the shadow ###
        angles = np.linspace(0,2.0*np.pi,self.nt+1)
        self.boundary_line = np.array((self.radius*np.cos(angles),self.radius*np.sin(angles)))

        self.fprint("Radius:        {: .2f}".format(self.radius))
        self.fprint("Center:       ({: .2f}, {: .2f})".format(self.center[0],self.center[1]))
        self.fprint("Meshing Type:  {0}".format(self.mesh_type))

        mesh_start = time.time()
        self.fprint("")
        if self.mesh_type == "mshr":
            raise NotImplementedError("Mshr is no longer supported, use gmsh instead.")
            # self.fprint("Generating Mesh Using mshr")

            # ### Create Mesh ###
            # mshr_circle = Circle(Point(self.center[0],self.center[1]), self.radius, self.nt)
            # self.mesh = generate_mesh(mshr_circle,self.res)


        elif self.mesh_type == "gmsh":
            self.fprint("Generating Mesh Using gmsh")

            if (self.params.rank == 0):

                # only dependencies for the gmsh case
                import gmsh
                import meshio
                import tempfile

                # start up gmsh model
                gmsh.initialize()
                gmsh.model.add("circmesh")
                factory = gmsh.model.occ

                # use arc division to set characteristic length on the horizon
                lc_mean = 2*np.pi*self.radius/self.nt

                # set up the geometry
                x_circ = self.center[0]
                y_circ = self.center[1]
                z_circ = 0.0

                # create the circle in the horizon plane
                ln_circ = factory.addCircle(x_circ, y_circ, z_circ, self.radius)
                cl_circ = factory.addCurveLoop([ln_circ])
                ps_circ = factory.addPlaneSurface([cl_circ])

                # synchronize the geometry
                factory.synchronize()

                # set the mesh size
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_mean)

                # generate and output
                gmsh.model.mesh.generate(2)

                # create a temporary directory to hold mesh IO files for conversion
                dir_for_meshes = tempfile.TemporaryDirectory()

                # write and finalize gmsh
                gmsh.write(os.path.join(dir_for_meshes.name, "dummy.msh"))
                gmsh.finalize() # don't need it anymore!

                # use meshio to convert from gmsh to dolfin xml file
                mesh_meshio = meshio.read(os.path.join(dir_for_meshes.name, "dummy.msh"))

                # flatten to 2D
                assert np.all(np.isclose(mesh_meshio.points[0,2], mesh_meshio.points[:,2]))
                mesh_meshio.points = mesh_meshio.points[:,0:2]

                # write the dolfin file, and then finally re-load it into dolfin
                meshio.write(
                    os.path.join(dir_for_meshes.name, "dummy.xml"),
                    mesh_meshio,
                    file_format="dolfin-xml",
                )

            # broadcast the temporary directory for shared use
            dirname_for_meshes = self.params.comm.bcast(dir_for_meshes.name, root=0)

            # load the xml into windse
            self.mesh = Mesh(os.path.join(dirname_for_meshes, "dummy.xml"))

            # delete the temporary directory
            if self.params.rank == 0:
                dir_for_meshes.cleanup()

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
        nom_x = np.cos(self.inflow_angle)
        nom_y = np.sin(self.inflow_angle)

        ### Define Boundary Subdomains ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Marking Boundaries")
        outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y)
        inflow  = CompiledSubDomain("nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, c0=self.center[0], c1=self.center[1])
        self.boundary_subdomains = [None,None,None,None,None,None,outflow,inflow]
        self.boundary_names = {"west":None,"east":None,"south":None,"north":None,"bottom":None,"top":None,"inflow":8,"outflow":7}
        self.boundary_types = {"inflow":  ["inflow"],
                               "no_stress": ["outflow"]}

        ### Generate the boundary markers for boundary conditions ###
        self.BuildBoundaryMarkers()

        ### Rotate Boundary
        if not near(self.inflow_angle,0.0):
            self.RecomputeBoundaryMarkers(self.inflow_angle)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def RecomputeBoundaryMarkers(self,inflow_angle):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Define Plane Normal ###
        nom_x = np.cos(inflow_angle)
        nom_y = np.sin(inflow_angle)

        # ### Define Boundary Subdomains ###
        # mark_start = time.time()
        # outflow = CompiledSubDomain("on_boundary", nx=nom_x, ny=nom_y)
        # inflow  = CompiledSubDomain("nx*(x[0]-c0)+ny*(x[1]-c1)<=0  && on_boundary", nx=nom_x, ny=nom_y, c0=self.center[0], c1=self.center[1])
        # self.boundary_subdomains = [None,None,None,None,None,None,outflow,inflow]

        # ### Generate the boundary markers for boundary conditions ###
        # self.BuildBoundaryMarkers()

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
        self.x_range = self.x_range*self.xscale
        self.y_range = self.y_range*self.xscale
        self.fprint("X Range: [{: .2f}, {: .2f}]".format(self.x_range[0]/self.xscale,self.x_range[1]/self.xscale))
        self.fprint("Y Range: [{: .2f}, {: .2f}]".format(self.y_range[0]/self.xscale,self.y_range[1]/self.xscale))
        self.dim = 2

        ### Create mesh ###
        mesh_start = time.time()
        self.fprint("")
        self.fprint("Generating Mesh")

        if self.mesh_type == "gmsh":

            if (self.params.rank == 0):

                # only dependencies for the gmsh case
                import gmsh
                import meshio
                import tempfile

                # start up gmsh model
                gmsh.initialize()
                gmsh.model.add("rectmesh")

                # true extents
                Lx_true = self.x_range[1] - self.x_range[0]
                Ly_true = self.y_range[1] - self.y_range[0]

                # either do stretching to preserve
                stretch_anisotropic = True # TODO: turn into an input param
                if stretch_anisotropic:
                    Lx_prime = 1.0
                    Ly_prime = self.ny/self.nx
                else:
                    Lx_prime = Lx_true
                    Ly_prime = Ly_true

                # create points that constitute box's lower surface
                pt_bx000 = gmsh.model.geo.addPoint(
                    self.x_range[0]/Lx_true*Lx_prime,
                    self.y_range[0]/Ly_true*Ly_prime,
                    0.0,
                )
                pt_bx100 = gmsh.model.geo.addPoint(
                    self.x_range[1]/Lx_true*Lx_prime,
                    self.y_range[0]/Ly_true*Ly_prime,
                    0.0,
                )
                pt_bx010 = gmsh.model.geo.addPoint(
                    self.x_range[0]/Lx_true*Lx_prime,
                    self.y_range[1]/Ly_true*Ly_prime,
                    0.0,
                )
                pt_bx110 = gmsh.model.geo.addPoint(
                    self.x_range[1]/Lx_true*Lx_prime,
                    self.y_range[1]/Ly_true*Ly_prime,
                    0.0,
                )

                # create lines of box's lower surface
                ln_bx000 = gmsh.model.geo.addLine(pt_bx000, pt_bx100)
                ln_bx001 = gmsh.model.geo.addLine(pt_bx100, pt_bx110)
                ln_bx010 = gmsh.model.geo.addLine(pt_bx010, pt_bx110)
                ln_bx011 = gmsh.model.geo.addLine(pt_bx000, pt_bx010)

                # create curve loop for lower surface
                cl_bx0 = gmsh.model.geo.addCurveLoop(
                    [ln_bx000, ln_bx001, -ln_bx010, -ln_bx011]
                )
                # create plane surface
                ps_bx0 = gmsh.model.geo.addPlaneSurface([cl_bx0])

                # synchronize the geometry
                gmsh.model.geo.synchronize()

                # characteristic length based on regular tetrahedron volume formula
                lc_mean = np.sqrt(4/np.sqrt(3)*(Lx_prime/self.nx)*(Ly_prime/self.ny))
                # set the mesh size to the characteristic length
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_mean)

                # generate and output
                gmsh.model.mesh.generate(2)

                # create a temporary directory to hold mesh IO files for conversion
                dir_for_meshes = tempfile.TemporaryDirectory()

                # write and finalize gmsh
                gmsh.write(os.path.join(dir_for_meshes.name, "dummy.msh"))
                gmsh.finalize() # don't need it anymore!

                # use meshio to convert from gmsh to dolfin xml file
                mesh_meshio = meshio.read(os.path.join(dir_for_meshes.name, "dummy.msh"))

                # flatten to 2D
                assert np.all(np.isclose(mesh_meshio.points[0,2], mesh_meshio.points[:,2]))
                mesh_meshio.points = mesh_meshio.points[:,0:2]

                meshio.write(
                    os.path.join(dir_for_meshes.name, "dummy.xml"),
                    mesh_meshio,
                    file_format="dolfin-xml",
                )

            # broadcast the temporary directory for shared use
            dirname_for_meshes = self.params.comm.bcast(dir_for_meshes.name, root=0)

            # load the xml into windse
            self.mesh = Mesh(os.path.join(dirname_for_meshes, "dummy.xml"))

            # delete the temporary directory
            if self.params.rank == 0:
                dir_for_meshes.cleanup()

            # stretch the coordinates back to the physical dimensions (now with the
            # specified aspect ratios)
            self.mesh.coordinates()[:,0] = self.mesh.coordinates()[:,0]*Lx_true/Lx_prime
            self.mesh.coordinates()[:,1] = self.mesh.coordinates()[:,1]*Ly_true/Ly_prime
        else:
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
        east   = CompiledSubDomain("near(x[0], x1) && on_boundary",x1 = self.x_range[1])
        north  = CompiledSubDomain("near(x[1], y1) && on_boundary",y1 = self.y_range[1])
        west   = CompiledSubDomain("near(x[0], x0) && on_boundary",x0 = self.x_range[0])
        south  = CompiledSubDomain("near(x[1], y0) && on_boundary",y0 = self.y_range[0])
        self.boundary_subdomains = [east,north,west,south]
        self.boundary_names = {"east":1,"north":2,"west":3,"south":4,"bottom":None,"top":None,"inflow":None,"outflow":None}
        self.boundary_types = {"inflow":    ["west","south","north"],
                               "no_stress": ["east"]}

        ### Generate the boundary markers for boundary conditions ###
        self.BuildBoundaryMarkers()

        ### Rotate Boundary
        if not near(self.inflow_angle,0.0):
            self.RecomputeBoundaryMarkers(self.inflow_angle)

        mark_stop = time.time()
        self.fprint("Boundaries Marked: {:1.2f} s".format(mark_stop-mark_start))
        self.fprint("Initial Domain Setup",special="footer")

    def RecomputeBoundaryMarkers(self,inflow_angle):
        mark_start = time.time()
        self.fprint("")
        self.fprint("Remarking Boundaries")

        ### Get ids for new markings ###
        inflow_id  = self.boundary_names[self.boundary_types["inflow"][0]]
        outflow_id = self.boundary_names[self.boundary_types["no_stress"][0]]
        east_id    = self.boundary_names["east"]
        north_id   = self.boundary_names["north"]
        west_id    = self.boundary_names["west"]
        south_id   = self.boundary_names["south"]

        ### Set up the baseline (angle=0) order ###
        cardinal_ids = [east_id,north_id,west_id,south_id]
        diagonal_ids = [outflow_id,outflow_id,inflow_id,inflow_id]

        ### Count the number of pi/2 sections in the new inflow_angle ###
        turns = inflow_angle/(pi/2)

        ### New order ###
        tol = 1e-3
        if turns % 1 <= tol: # we are at a cardinal direction
            new_order = np.roll(cardinal_ids,int(turns))
        else:
            new_order = np.roll(diagonal_ids,int(turns))

        ### Get a list of all wall facets ###
        wall_facets = []
        for i in cardinal_ids:
            wall_facets += self.boundary_markers.where_equal(i)

        ### Iterate through facets remarking as prescribed by new_order ###
        for facet_id in wall_facets:

            ### Get the facet normal
            facet = Facet(self.mesh,facet_id)
            facet_normal = facet.normal().array()
            nx = np.sign(facet_normal[0])*(abs(facet_normal[0])>tol)
            ny = np.sign(facet_normal[1])*(abs(facet_normal[1])>tol)

            ### Map the normal to the order east, north, west, south
            face_id = int(abs(nx+2*ny-1))

            ### remark the boundary ###
            self.boundary_markers.set_value(facet_id,new_order[face_id])

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
        # raise NotImplementedError("Imported Domains need to be updated. Please use an Interpolated domain for now.")
        super(ImportedDomain, self).__init__()

        self.fprint("Importing Domain",special="header")

        ### Import data from Options ###
        if self.path is not None:
            if self.filetype == "xml.gz":
                self.mesh_path  = self.path + "mesh.xml.gz"
                self.boundary_path = self.path + "boundaries.xml.gz"
            elif self.filetype == "h5":
                self.mesh_path = self.path + "mesh_data.h5"
            if self.interpolated and self.terrain_path is not None:
                self.terrain_path = self.path + "terrain.txt"

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
        self.mesh.coordinates()[:] *= self.xscale


        self.bmesh = BoundaryMesh(self.mesh,"exterior")
        mesh_stop = time.time()
        self.dim = self.mesh.topology().dim()

        if self.dim != 3:
            raise ValueError("Currently, only 3D meshes can be imported.")

        self.fprint("Mesh Imported: {:1.2f} s".format(mesh_stop-mesh_start))

        ### Calculate the range of the domain and push to options ###
        self.x_range = [min(self.mesh.coordinates()[:,0]),max(self.mesh.coordinates()[:,0])]
        self.y_range = [min(self.mesh.coordinates()[:,1]),max(self.mesh.coordinates()[:,1])]
        self.z_range = [min(self.mesh.coordinates()[:,2]),max(self.mesh.coordinates()[:,2])]
        self.params["domain"]["x_range"] = [min(self.mesh.coordinates()[:,0])/self.xscale,max(self.mesh.coordinates()[:,0])/self.xscale]
        self.params["domain"]["y_range"] = [min(self.mesh.coordinates()[:,1])/self.xscale,max(self.mesh.coordinates()[:,1])/self.xscale]
        self.params["domain"]["z_range"] = [min(self.mesh.coordinates()[:,2])/self.xscale,max(self.mesh.coordinates()[:,2])/self.xscale]

        ### Load the boundary markers ###
        mark_start = time.time()
        self.fprint("")
        self.fprint("Importing Boundary Markers")
        if self.filetype == "h5":
            self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.geometry().dim()-1)
            hdf5.read(self.boundary_markers, "/boundaries")
        elif self.filetype == "xml.gz":
            self.boundary_markers = MeshFunction("size_t", self.mesh, self.boundary_path)
        self.boundary_names = {"east":1,"north":2,"west":3,"south":4,"bottom":5,"top":6,"inflow":None,"outflow":None}
        self.boundary_types = {"inflow":    ["west","south","north"],
                               "no_slip":   ["bottom"],
                               "free_slip": ["top"],
                               "no_stress": ["east"]}
        mark_stop = time.time()
        self.fprint("Boundary Markers Imported: {:1.2f} s".format(mark_stop-mark_start))

        ### Create the interpolation function for the ground ###
        interp_start = time.time()
        self.fprint("")
        self.fprint("Building Interpolating Function")

        if self.interpolated:
            self.SetupInterpolatedGround()

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

        if self.analytic:
            self.SetupAnalyticGround()
        else:
            self.SetupInterpolatedGround()

        interp_stop = time.time()
        self.fprint("Interpolating Function Built: {:1.2f} s".format(interp_stop-interp_start),special="footer")

    def Finalize(self):
        self.Move(self.ground_function)
        DefaultFinalize(self)

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

        if self.analytic:
            self.SetupAnalyticGround()
        else:
            self.SetupInterpolatedGround()

        interp_stop = time.time()
        self.fprint("Ground Function Built: {:1.2f} s".format(interp_stop-interp_start),special="footer")

    def Finalize(self):
        self.Move(self.ground_function)
        DefaultFinalize(self)


class PeriodicDomain(BoxDomain):
    def __init__(self):
        super(PeriodicDomain, self).__init__()


        # self.extra_forcing_term

