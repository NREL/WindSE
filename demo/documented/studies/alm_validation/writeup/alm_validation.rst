
.. _alm_validation:

Actuator Line Method Validation
===============================

In this study, we compare four key along-blade forces calculated by WindSE's actuator line method---namely lift, drag, angle of attack, and axial velocity---to a benchmark 2018 study by Martínez-Tossas et al. [1]_


Keywords:
---------

actuator line method, ALM, lift, drag, angle of attack, velocity, validation



Input files and code version:
-----------------------------

The complete details of the domain, mesh, and turbine can be found in Martínez-Tossas et al. [1]_ but are reproduced as a WindSE study in the following YAML file. 

    * Parameter File: :download:`../../../demo/documented/studies/alm_validation/input_files/alm_validation.yaml`
    * Code Version: `WindSE s_2022.10.09 <https://github.com/NREL/WindSE/releases/tag/s_2022.10.09>`_

.. Note: not all studies need a full release, but they should at least link to a tag/commit.



Setup:
------

The setup used for this test and enumerated in the YAML file reproduces almost exactly the setup of Martínez-Tossas et al. [1]_ The turbine used is the NREL 5-MW reference turbine which has a rotor diameter, :math:`D`, of 126 m operating with rotational speed 9.155 RPM. The computational domain spans :math:`-3D \leq x \leq 21D`, :math:`-3D \leq y \leq 3D`, :math:`-3D \leq z \leq 3D` where the rotor is centered at :math:`(0, 0, 0)` and oriented normal to the :math:`x^+` flow direction. The grid size surrounding the rotor is 1.96875 m, which implies 64 mesh cells across the rotor diameter. 64 actuator nodes are used along the length of each blade with a fixed Gaussian size of 10 m. The fluid is assigned a density :math:`\rho=1.0` kg/m\ :sup:`3`, the inflow velocity is a uniform 8 m/s applied at the :math:`x^-` wall, slip boundary conditions (no flow through) are set at all lateral walls, and a 0-pressure outlet condition applied at the :math:`x^+` wall. The Smagorinsky eddy viscosity model is used with :math:`C_s=0.16`.

Note that because we are only comparing along-blade quantities in this study and not features in the far wake, we make a slight deviation from the paper and perform *local* mesh refinements to acheive the specified 1.96875 m grid size in the region surrounding the turbine. This reproduces the mesh size and resolution around the rotor as specified in the paper but leaves relatively coarse mesh cells in the downstream region to reduce the computational workload. By this same token, we only run the simulations long enough for the along-blade quantities to converge, rather than the much longer time needed to satisfy repeated flow throughs, which was experimentally found to be 100 s.
    


Results:
--------

The results generated by the WindSE ALM implementation are plotted alongside 4 other codes, including NREL's SOWFA, in Figure 1:

.. figure:: ../../../demo/documented/studies/alm_validation/writeup/windse_alm_verification_3.png
   :scale: 80 %

   Figure 1: Clockwise from top left, the angle of attack, relative axial velocity, lift, and drag as a function of position along the blade calculated using two different discretizations: the 64-actuator version outlined above and a 12-acuator version testing a greatly reduced number of control points that would be more suitable for an optimization problem (opt scale) featuring actuator-based controls.

Both the 64-node and 12-node, optimization scale, cases recover the along-blade forces well. The axial velocity is normalized by the reference velocity, :math:`U_{ref} = 8` m/s, and so represents the percentage of freestream while the lift and drag are non-dimensionalized by the actuator length, :math:`w`, rotor diameter, density, and reference velocity.  


Conclusions:
------------

We find that along-blade forces are recovered very well using WindSE's actuator line implementation. The largest deviations away from the comparison codes seem to occur in regions of sharp change (see the first 20% of blade span drag profile) where differences in the airfoil section resolution and the interpolation used between airfoil sections may have a large effect. 

With respect to sizing ALM simulations, we are able to produce very good agreement with both the benchmark 64-node setup and the 12-node, optimization scale, setup. Although it is important to keep in mind that we are unlikely to produce such good agreement in the downstream wake using this study, for the purposes of optimizing objectives confined to the rotor plane (e.g., a single turbine's power) with respect to along-the-blade quantities like chord or twist angle, using ~10 actuator nodes sized to be roughly ~1/10 of the rotor diameter seems to be a good starting point. For a more detailed convergence study of actuator node resolution and size, see the 2017 paper of Martínez-Tossas et al. [2]_


References:
-----------

.. [1]  Luis A. Martínez-Tossas, Matthew J. Churchfield, Ali Emre Yilmaz, Hamid Sarlak, Perry L. Johnson, Jens N. Sørensen, Johan Meyers, and Charles Meneveau, "Comparison of four large-eddy simulation research codes and effects of model coefficient and inflow turbulence in actuator-line-based wind turbine modeling", Journal of Renewable and Sustainable Energy 10, 033301 (2018) https://doi.org/10.1063/1.5004710

.. [2]  Martínez-Tossas, L. A., Churchfield, M. J., and Meneveau, C. (2017) Optimal smoothing length scale for actuator line models of wind turbine blades based on Gaussian body force distribution. Wind Energ., 20: 1083– 1096. doi: 10.1002/we.2081.
