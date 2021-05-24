from os.path import dirname, basename, isfile, join
import glob
import importlib

### Make sure only the objective dictionary is imported ###
__all__ = ["objectives_dict"]

### Get all files in folder ###
files = glob.glob(join(dirname(__file__), "*.py"))

### Add the names to a dictionary ###
objectives_dict = {}
for file in files:

    ### Filter and define which should be imported ###
    if isfile(file) and not basename(file).startswith('_'):
        objective_file = "windse.objective_functions."+basename(file)[:-3]

        ### Try to import the two required elements ###
        try:
            imported_objective = importlib.import_module(objective_file)
            name = imported_objective.name
            objective = imported_objective.objective
        except:
            raise ValueError("Objective File '"+file+"' is missing name or objective(), check _template_.py for details")

        ### Check to make sure we don't have conflicting names ###
        if name in objectives_dict.keys():
            raise ValueError("Two objectives named: "+name+". Please rename one")

        ### If all is good, add the objective to the file
        objectives_dict[name] = objective
