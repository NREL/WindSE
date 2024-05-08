.. _param_tips:

Additional Notes
================

This file is intended to supplement the :doc:`param_api` to give more 
explanation on some of the more involved inputs. 

.. contents:: :local:



.. _imported_domain:

Formatting for Importing a Domain
---------------------------------

When using ``domain:type:imported``, three files are required: 

* mesh.xml.gz - this contains the mesh in a format dolfin can handle
* boundaries.xml.gz - this contains the facet markers that define where the boundaries are
* terrain.txt - this contains the data for the terrain as a function of x and y. 

These files can be specified individually/manually using: ``mesh_path``, ``bound_path``, and ``terrain_path``,
or have WindSE find this files automatically by only setting the ``path`` option, which expects
the three files above with the exact same filenames and types.

The topology file assumes that the coordinates are from a uniform mesh.
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

Note: If using "h5" file format, the mesh and boundary will be in one file.



.. _imported_windfarm:

Formatting for Importing a Wind Farm
------------------------------------

To import a wind farm, set the path to a .csv file containing the per 
turbine information. In the .csv file, each column specifies a turbine
property and each row is a unique turbine. At minimum, the locations for
each turbine must be specified. Here is a small two turbine example::

         x,      y
    200.00, 0.0000
    800.00, 0.0000

Additional turbine properties can be set by adding a column with a header
equal to the yaml parameter found in the "Turbine Options" section. Here
is an example of a two turbine farm with additional properties set::

      x,       y,      HH,      yaw,      RD, thickness,   axial
    0.0,  -325.0,   110.0,   0.5236,   130.0,      13.0,    0.33
    0.0,   325.0,   110.0,  -0.5236,   130.0,      13.0,    0.33

The columns can be in any order and white space is ignored. If a property
is set in both the yaml and the imported .csv, the value in the .csv will
be used and a warning will be displayed.



.. _turbine_representation:

Supported Turbine Representation
--------------------------------

.. todo::
    Describe in details the difference between ``disks``, ``numpy_disks``, ``hybrid_disks``, ``lines``, ``dolfin_lines``, and ``disabled``

There are three types of supported turbine representations:
    
    1. ``disks``: Standard actuator disks
    2. ``hybrid_disks``: Actuator disk build from a series of nested ring allow for more control over the forcing function and include rotation.  
    3. ``lines``: Standard actuator lines to be used with unsteady simulations. 
    4. ``disabled``: These still have turbine location data which is useful for mesh refinement but the forcing function is set to zero.


.. _mesh_refinement:

Mesh Refinement Options
-----------------------

The domain can be refined in special
ways to maximize the efficiency of the number DOFs. None of these options
are required. There are three types of mesh manipulation: warp, farm refine,
turbine refine. Warp shifts more cell towards the ground, refining the farm
refines within the farm extents, and refining the turbines refines within
the rotor diameter of a turbine. When choosing to warp, a "smooth" warp will 
shift the cells smoothly towards the ground based on the strength. A "split"
warp will attempt to create two regions, a high density region near the 
ground and a low density region near the top



.. _refine_custom:

Custom Mesh Refinement
----------------------

The ``refine_custom`` option can be used to supply WindSE with a series of
refinement step that will be performed before the standard "farm" and "turbine"
refinements. This options is used by supplying a list of refinements in order 
with the corresponding refinement parameters. Example::

    refine_custom:
        box:
            x_range: [-500,500]
            y_range: [-500,500]
            z_range: [0,150]
        cylinder:
            center: [0, 0, 0]
            radius: 1000.0
            height: 250.0
        simple:
            radius: 100

The example up above will result in three refinements:

    1. A box refinement bounded by: [[-500,500],[-500,500],[0,150]]
    2. A cylinder centered at origin with radius 750 m and a height of 150 m
    3. A simple turbine refinement with radius 100 m 

If you want to repeat a specific type of refinement (such as two box refinements),
you need to start the entry with a number to help differentiate the two refinements. 
Example::

    refine_custom:
        1:
            type: box
            x_range: [-850, 1200]
            y_range: [-150, 150]
            z_range: [0, 200]
        2:
            type: box
            x_range: [-800, 1200]
            y_range: [-100, 100]
            z_range: [25, 175]

Currently there are 7 types of refinement:
    
    1. ``full``: perform a uniform refinement over the full domaint
    2. ``box``: refine within a box defined by the x, y, z ranges
    3. ``cylinder``: refine within a cylinder oriented parallel with the z-axis
    4. ``stream``: refine within a cylinder oriented parallel with the x-axis and rotated by theta around the pivot offset
    5. ``simple``: for each turbine, refine in a cylinder oriented parallel with the z-axis centered on that turbine
    6. ``tear``: for each turbine, refine in a tear drop shape mimicking the way that can be rotated by theta centered on that turbine
    7. ``wake``: for each turbine, refine in a cylinder oriented parallel with the x-axis and rotated by theta centered on that turbine

The syntax for each refinement type is::

    refine_custom:
        full: {}
        box:
            x_range: [x_min,x_max]
            y_range: [y_min,y_max]
            z_range: [z_min,z_max]
            expand_factor: <float>
        cylinder:
            center: [c_x,c_y,c_z]
            radius: <float>
            height: <float>
            expand_factor: <float>
        stream: 
            center: [c_x,c_y,c_z]
            radius: <float>
            height: <float>
            theta: <float>
            pivot_offset: <float>
            expand_factor: <float>
        simple:
            radius: <float>
            expand_factor: <float>
        tear: 
            radius: <float>
            theta: <float>
            expand_factor: <float>
        wake: 
            radius: <float>
            length: <float>
            theta: <float>
            expand_factor: <float>

.. note::
    * For ``cylinder``, the center is the base of the cylinder
    * For ``stream``, the center is the start of the vertical base and offset indicates the rotation offset
    * For ``stream``, wake, length is the distance center to the downstream end of the cylinder
    * For ``stream``, tear, wake, theta rotates the shape around the center



.. _inflow_angle:

Defining Inflow Angle
---------------------

The inflow angle is defined in the meteorological sense with units of degrees. This means an inflow angle of 0 will result in wind blowing from the north to the south. Increasing the inflow angle will rotate the wind clockwise (i.e. 90 degrees corresponds to wind blowing from the east and moving west). The default value of 270 degrees, which corresponds to wind blowing from the west to the east. 

.. _custom_boundaries:

Customizing Boundary Conditions
-------------------------------

This section describes the boundary condition options. There are four types
of boundary conditions: inflow, no slip, free slip, no stress. By default, inflow is 
prescribed on boundary facing into the wind, no slip on the ground and 
no stress on all other faces. These options describe the inflow boundary
velocity profile. 

If you are importing a mesh or want more control over boundary conditions, you can specify the boundary markers using ``names`` and ``types``.
The default for these two are

Rectangular Mesh::

    boundary_condition:
        boundary_names: 
            east:  1
            north: 2
            west:  3
            south: 4
        boundary_types: 
            inflow:    ["west","north","south"]
            no_stress: ["east"]

Box Mesh::

    boundary_condition:
        boundary_names: 
            east:   1
            north:  2
            west:   3
            south:  4
            bottom: 5
            top:    6
        boundary_types: 
            inflow:    ["west","north","south"]
            free_slip: ["top"]
            no_slip:   ["bottom"]
            no_stress: ["east"]

Circle Mesh::

    boundary_condition:
        boundary_names: 
            outflow: 7
            inflow:  8
        boundary_types: 
            inflow:    ["inflow"]
            no_stress: ["outflow"]

Cylinder Mesh::

    boundary_condition:
        boundary_names: 
            outflow: 5
            inflow:  6
            bottom:  7
            top:     8
        boundary_types: 
            inflow:    ["inflow"]
            free_slip: ["top"]
            no_slip:   ["bottom"]
            no_stress: ["outflow"]

These defaults correspond to an inflow wind direction from West to East.

When marking a rectangular/box domains, from a top-down perspective, start from 
the boundary in the positive x direction and go counter clockwise, the boundary 
names are: "easy", "north", "west", "south". Additionally, in 3D there are also
"top" and "bottom". For a circular/cylinder domains, the boundary names are
"inflow" and "outflow". Likewise, in 3D there are also "top" and "bottom". 
Additionally, you can change the ``boundary_types`` if using one of the built 
in domain types. This way you can customize the boundary conditions without 
importing a whole new mesh.


.. _multi_objectives:

Defining Multiple Objective Functions
-------------------------------------

.. todo::
    automatically compile list of objective function and their keyword arguments.
    For now, just reference :meth:`windse.objective_functions`.

The ``objective_type`` can be defined in two ways. First as a single string such as::

    optimization:
        objective_type: alm_power 

If the object chosen in this way has any keyword arguments, the defaults will 
automatically chosen. If you want to use multiple objective, each one need its
own entry, which allow for setting keyword arguments::

    optimization:
        objective_type: 
            1:
                type: point_blockage
                location: [-390.0,0.0,110.0]
            2:
                type: point_blockage
                location: [0.0,0.0,240.0]
            3:
                type: box_blockage
                p0: [-990, -430, 0]
                p1: [-860, 430,  200]
            4:
                type: box_blockage
                p0: [-730, -430, 240]
                p1: [730,  430,  370]
            5:
                type: cyld_kernel
                orientation: upstream
            6:
                type: cyld_kernel
                orientation: above

For a full list of objective function visit: :meth:`windse.objective_functions`. Notice that we can have multiple version of the same objective starting the dictionary with an integer. This allows us to evaluate objectives of the same type with different keyword arguments. If, instead, each object was unique, then they could be defined as::

    optimization:
        objective_type: 
            point_blockage:
                location: [-390.0,0.0,110.0]
            box_blockage:
                p0: [-990, -430, 0]
                p1: [-860, 430,  200]
            cyld_kernel:
                orientation: upstream

Regardless of the number of objective types listed, currently, only the first one will be used for an optimization. 



.. _multi_constraints:

Defining Multiple/Custom Constraints
------------------------------------

The ``constraint_types`` option is defined in a similar way to the ``objective_type``. By default the minimum distance between turbines is setup::

    constraint_types:
        min_dist:       
            target: 2   
            scale:  1   

This constraint will only be used if the ``control_types`` contains "layout". Additional constraints can be added using the same objective functions from :meth:`windse.objective_functions` by setting::

    constraint_types:
        min_dist:       
            target: 2   
            scale:  1 
        plane_blockage:
            target: 8.0
            scale: -1
            kwargs:
                axis: 2
                thickness: 130
                center: 240.0

This will still enforce the layout constraint but will additionally enforce a "plane_blockage" type constraint. By default, the constrains are setup like:

.. math::

    s * \left( c(m)-t \right) \geq 0

where :math:`c` is the constraint function, :math:`t` is the target, :math:`s` is the scale, and :math:`m` are the controls. In this configuration, we are enforcing that the result of the constraint function is greater than or equal to the target. However, we can set the scale to -1 to flip the inequality. Just like the ``objective_type``, multiple constraints of the same type can be use by appending "_#" followed by a number to the end of the name with the exception of the "min_dist" type. 