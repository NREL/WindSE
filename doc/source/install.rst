Installation
============

General Installation:
---------------------

In order to use WindSE, the python version fenics 2018.1.0 or later must be installed along with a compatible version of dolfin-adjoint. WindSE can be installed by downloading the source files from the GitHub (https://github.com/NREL/WindSE) and running the command::

    pip install -e .

in the root folder. 

Conda Installation (Automatic):
-------------------------------

The best way to use WindSE is via a conda environment. To install conda check out this link: `Conda Installation. <https://conda.io/projects/conda/en/latest/user-guide/install/>`_ After conda is installed, you can automatically setup the WindSE environment using::

    conda create --name <enviroment_name>
    conda activate <enviroment_name>
    conda install -c conda-forge windse

Conda Installation (Source):
----------------------------

If you want to use the latest version or just want to setup the environment manually, follow these steps. After conda is installed, create a new environment using::

    conda create --name <enviroment_name>

You can replace the name <enviroment_name> with a different name for the environment if you want. Next we activate the environment using::

    conda activate <enviroment_name>

or whatever you named your environment. Now we need to install the dependent packages using::

    conda install -c conda-forge fenics=2018.1.0 dolfin-adjoint mshr matplotlib scipy pyyaml

Finally, download/clone the WindSE repo and run::

    pip install -e .

in the root folder. 

