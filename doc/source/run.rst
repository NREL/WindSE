Running WindSE
==============

To run WindSE, first create a parameters file (as described in :doc:'params' and 
demonstrated in :doc:'demo'). Then activate the conda environment using::

    source activate fenics_windse

then run::

    windse_driver run <params file>

Where, <params files> is the path of the parameters file you wish to run.
By default windse searches for "params.txt" in the current directory if 
no file is supplied.   

Sit back and let the magic happen.