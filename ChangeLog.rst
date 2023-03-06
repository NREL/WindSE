ChangeLog:
==========

2021.08.01:
-----------

Versions:
~~~~~~~~~

* FEniCS: 2019.1.0, build: py38_9
* Dolfin-Adjoint: 2019.1.0

New Features:
~~~~~~~~~~~~~

* Refactored how wind farms and turbines work
* Added DELs for the actuator line
* Added small yaw optimization regression test

Bug Fixes
~~~~~~~~~

* Change boundary velocity reference height to be based on mean(turbines.z)
* Reduced the size of 3D_Skew_Steady.yaml problem
* No longer annotate debug functions
* Fixed a bug where the forward problem was not recomputed correctly during optimization 

2021.08.00:
-----------

Versions:
~~~~~~~~~

* FEniCS: 2019.1.0, build: py38_9
* Dolfin-Adjoint: 2019.1.0

New Features:
~~~~~~~~~~~~~

* Added Regression Tests including Coveralls
* Switched to Github Actions
* Added output for angle of attack and rotor plane forces
* Added Empty and disabled wind farm options
* Linked OpenMDAO and SNOPT
* Parallelized Unsteady solver
* Reorganized Objective functions
* Added ability to save multiple objectives at once
* Added blockage based objective functions
* Grid farm can now be defined using spacing and skew
* Random farm now respect a minimum distance between turbines
* Added the option to use a corrective force for momentum loss
* Added Iterative Steady solvers
* Added example scripts for small jobs 
* Added ability to maximize or minimize


Bug Fixes
~~~~~~~~~

* Fixed mpi optimizations
* Normalized blockage objective functions to report in m/s



2020.5.0:
---------

Versions:
~~~~~~~~~

* FEniCS: 2019.1.0
* Dolfin-Adjoint: 2019.1.0
* Requires: tsfc `install reference: <https://fenics.readthedocs.io/projects/ffc/en/latest/installation.html>`_

New Features:
~~~~~~~~~~~~~

* Added Actuator Line Model
* Added a Parameter checker
* Added Unit Tests
* Added 2.5D options
* Added Power Functional optimized for 2D
* Added Numpy Representation for Actuator Disks
* Driver functions can be used externally

Bug Fixes
~~~~~~~~~

* Fixed import for Jupyter
* Fixed boundary label typo
* Fixed free slip boundary assignment