"""
The PostprocessingManager submodule contains the various classes used for 
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
    import os
    from sys import platform
    import numpy as np
    from pyadjoint.tape import stop_annotating
    import yaml

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *
        
    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt


def write_dict_to_yaml(yaml_filename, dict_to_write):
    # This custom representer ensures that numpy ndarray data is written 
    # without the extra binary/metadata notation that numpy would normally include
    def ndarray_representer(dumper: yaml.Dumper, array: np.ndarray) -> yaml.Node:
        return dumper.represent_list(array.tolist())

    yaml.add_representer(np.ndarray, ndarray_representer)

    with open(yaml_filename, "w") as fp:
        print(f"Writing file {yaml_filename}")
        yaml.dump(dict_to_write, fp, sort_keys=False)


def write_to_floris(data_to_write, solver):
    """ Write the desired outputs as a FLORIS input-file.

    This function creates a dictionary according to FLORIS formats
    (e.g., "farm":"layout_x": [-10, 0, 10]) and then stores that dictionary as
    a yaml file. This enables a FLORIS simulation to be launched using the
    outputs from a previously-computed WindSE solve.
    
    Args:
        data_to_write (list): A list of the WindSE variables that should be
        written into a FLORIS input file solver

        (windse.solver): A WindSE solver object
    """
    floris_dict = {}

    if "layout" in data_to_write:
        if "farm" not in floris_dict:
            floris_dict["farm"] = {}

        layout = solver.problem.farm.get_hub_locations()

        floris_dict["farm"]["layout_x"] = layout[:, 0]
        floris_dict["farm"]["layout_y"] = layout[:, 1]

    if "yaw" in data_to_write:
        if "farm" not in floris_dict:
            floris_dict["farm"] = {}

        yaws = solver.problem.farm.get_yaw_angles()

        floris_dict["farm"]["yaw"] = np.degrees(yaws)

    # Generate the full path to the FLORIS file
    path_to_floris_file = os.path.join(solver.params.folder,
        solver.params["postprocessing"]["write_floris_filename"])

    # Finally, write the yaml file
    write_dict_to_yaml(path_to_floris_file, floris_dict)

