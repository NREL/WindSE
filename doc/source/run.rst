Running WindSE
==============

To run WindSE, first create a parameters file (as described in :doc:'params' and 
demonstrated in :doc:'demo'). Then activate the conda environment using::

    source activate <enviroment_name>

where the <enviroment_name> is what was defined in the install process
then run::

    windse run <params file>

Where, <params files> is the path of the parameters file you wish to run.
By default windse searches for "params.txt" in the current directory if 
no file is supplied.   

Sit back and let the magic happen. Additionally, you can run::

    windse run <params file> -p group:option:value

where ``group:option:value`` is a single unbroken string and the group is
the group in the params file, the option is the specific option in that
group and the value is the new value. This allows for overriding parameters
in the yaml file via the terminal. For example: ``wind_farm:HH:140`` will 
change the hub height of all turbines in a "grid" or "random" farm to 140 m
regardless of what was defined in the params.yaml file. 
