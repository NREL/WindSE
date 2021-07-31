'''
The objective function live in the windse/objective_functions folder.
These functions should be called using the dictionary
objective_funcs[<name>](solver, \*args, \*\*kwargs), where <name> is the
function name. 
'''


from os.path import dirname, basename, isfile, join
import glob
import importlib
from pyadjoint import stop_annotating, annotate_tape
import __main__
import copy

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = basename(__main__.__file__)
else:
    main_file = "ipython"
    
### Make sure only the objective dictionary is imported ###
__all__ = ["_annotated_objective", "objective_functions", "objective_kwargs"]

### Get all files in folder ###
files = glob.glob(join(dirname(__file__), "*.py"))

### Create a function that will turn off annotation where needed ###
def _annotated_objective(objective, *args, **kwargs):
    '''
    This is a wrapper function that allows dolfin_adjoint to record the objective functions to the tape
    '''
    annotate = annotate_tape(kwargs)
    if annotate:
        out = objective(*args, **kwargs)
    else:
        with stop_annotating():
            out = float(objective(*args, **kwargs))
    return out

### Add the names to a dictionary ###
objective_functions = {}
objective_kwargs = {}
for file in files:

    ### Filter and define which should be imported ###
    if isfile(file) and not basename(file).startswith('_'):
        objective_file = "windse.objective_functions."+basename(file)[:-3]

        ### Try to import the two required elements ###
        try:
            imported_objective = importlib.import_module(objective_file)
            name = imported_objective.name
            _objective = imported_objective.objective
            default_kwargs = imported_objective.keyword_defaults
        except:
            raise ValueError("Objective File '"+file+"' is missing name or objective(), check _template_.py for details")

        ### Check to make sure we don't have conflicting names ###
        if name in objective_functions.keys():
            raise ValueError("Two objectives named: "+name+". Please rename one")

        ### If all is good, add the objective to the file
        objective_functions[name] = _objective
        objective_kwargs[name] = default_kwargs

### Take a couple extra steps so sphinx will detect the objective functions
if main_file == "sphinx-build":
    for name, _func in objective_functions.items():
        vars()[name] = _func
        vars()[name].__module__ = __name__
