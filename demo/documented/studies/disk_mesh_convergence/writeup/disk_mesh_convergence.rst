
.. _disk_mesh_convergence:

Actuator Disk Mesh Convergence
==============================

This study is designed to develop intuition on how refined the actuator disk model needs to be to produce converged power. 



Keywords:
---------

mesh, actuator disk, power



Input files and code version:
-----------------------------

This study was ran using this parameter file and code version: 
    
    * Parameter File: :download:`../../../demo/documented/studies/disk_mesh_convergence/simulation/power.yaml`
    * Code Version: `WindSE 2021.08.01 <https://github.com/NREL/WindSE/releases/tag/2021.08.01>`_



Setup:
------

.. figure:: ../../../demo/documented/studies/disk_mesh_convergence/writeup/wind_farm.png
    :width: 600

    Figure 1: The wind farm layout

This simulation starts with a 3x3 grid arrangement of turbines with 3 rotor diameters padding for the inflow, outflow and sides as seen in Figure 1. The initial mesh has 16x16x10 cells in the x, y, and z directions, respectively. The mesh is then refined up to 3 times to get the mesh seen in Figure 2. Each refinement is local in a cylinder centered on each turbine with a radius of 1.25 time the rotor diameter and extending the full height of the turbine. For each level of refinement, a steady RANS simulation is performed with log layer inflow with hub height inflow speed of 8 m/s with the wind blowing from west to east. The number of refinement was controlled using the command line override parameter::

     windse run power.yaml -p general:name:n=N -p refine:turbine_num:N

where N is the number of refinements. 


.. figure:: ../../../demo/documented/studies/disk_mesh_convergence/writeup/3x3_mesh.png
    :width: 400

    Figure 2: The mesh after 3 turbine refinements


Results:
--------

The full compiled power output for each refinement level can be found in Table 1. The "mesh spacing" column is calculated by taking the full width of the mesh (1512 m) and dividing it by the initial number of cells in the x direction (16), which results in 94.5 m. This a measurement of the distance between mesh nodes. If we divide the rotor diameter by the mesh spacing, we get an approximation for the number of mesh nodes that span an actuator disk. This number is useful for determining how well resolved the disks are with a given resolution. For example, after 3 turbine refinements a disk is represented by about 10 nodes in the mesh. The goal of this study is to determine how many nodes per turbine is necessary to produced converged power calculations. 

.. csv-table:: Table 1: Power data for each turbine and refinement level
   :file: compiled_data.csv
   :widths: 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70
   :header-rows: 1

.. note::
    
    The magnitude of the power produced is not part of this study and has not been calibrated. This study is exclusively looking at mesh convergence. 


First let's look a the total power produced in Figure 3. Looking at this farm scale metric implies that convergence is essentially reached after one level of refinement. Looking at the "nodes/turbine" column in Table 1, this corresponds to needing only ~3 nodes per turbine, which seems exceptionally low. 

.. figure:: ../../../demo/documented/studies/disk_mesh_convergence/writeup/total_power.png
    :width: 600

    Figure 3: Total power with respect to mesh resolution

The story doesn't end there though. We can also look at the power produced by the leading turbines and fully waked turbines. Because the wind is blowing west to east the leading turbine are numbers 0, 3, 6 and the turbines we are calling "fully waked" are numbers 2, 5, 8. The leading turbine power shown in Figure 4 shows that convergence is a bit slower than the full farm's power taking an additional refinement. It is also interesting to note that all three turbines converge to the same power. This is expected because the inflow profile is not turbulent and uniform in the y direction so each of the leading turbines should experience the exact same forces resulting in identical powers. 

.. figure:: ../../../demo/documented/studies/disk_mesh_convergence/writeup/front_row.png
    :width: 600

    Figure 4: Power of the leading edge of turbines 

Finally, looking at the fully waked power production in Figure 5, we see a completely different trend. Now it is possible that these turbines are not yet fully converged. It appears that the power is converging but might require an additional refinement for a total of 4. Currently all of these simulation are running on a laptop, which does not have enough memory to run the 4 refinement simulation. This implies that if the wakes are exceptionally important to the simulation, more refinement is required. That said, after only 3 refinement, all three fully waked turbines produce the same power, just like the leading turbines. Since this is also expected, this could indicate 3 refinements or about 10 nodes per turbine is sufficient. 


.. figure:: ../../../demo/documented/studies/disk_mesh_convergence/writeup/back_row.png
    :width: 600

    Figure 5: Power of the fully wake turbines



Conclusions:
------------

Based on the information presented in this study, we conclude that when performing a steady simulation with actuator disk, aim for around 10 mesh nodes per turbine to get the best computational performance to accuracy. Some future studies that would be useful to refine this recommendation would include investigating mesh convergence of power with respect to:

    * turbulent inflow
    * increasing number of waked turbines
    * refining waked turbines more than leading turbines

