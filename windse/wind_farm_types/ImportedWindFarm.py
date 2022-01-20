from . import GenericWindFarm
import numpy as np
import shutil
import pandas

class ImportedWindFarm(GenericWindFarm):
    """
    A ImportedWindFarm produces turbines located based on a text file.
    The params.yaml file determines how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                imported: true
                path: "inputs/wind_farm.csv"

        The "wind_farm.csv" needs at least these columns with the header::

                 x,      y
            200.00, 0.0000
            800.00, 0.0000

        Each row defines a different turbine. The columns can be in any order and the
        header identifies which additional turbine property is being defined. For 
        additional turbine properties, the header is the same as the parameter defined
        in the "turbines" second of the yaml file.

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):

        self.name = "Imported Farm"
        super(ImportedWindFarm, self).__init__(dom)

        ### special stuff here ###

    def load_parameters(self):

        ### Get wind farm file location ###
        self.path = self.params["wind_farm"]["path"]

    def compute_parameters(self):
        ### Copy Files to input folder ###
        shutil.copy(self.path,self.params.folder+"input_files/")

        ### Clean Header (if needed) to remove leading "#" ###
        f = open(self.path, 'r')
        header = f.readline()
        header = header.replace("#", " ")
        header = header.replace(",", " ")

        ### set the separator for legacy wind farm files ###
        file_type = self.path.split(".")[-1]
        if file_type == "csv":
            sep = ","
        elif file_type == "txt":
            sep = "\s+"
        else:
            raise ValueError(f"Unknown imported wind farm file type: {file_type}")

        ### Import the data from path ###
        self.imported_params = pandas.read_csv(f, sep=sep, names=header.split(), dtype=float)

        ### Set the number of turbines ###
        self.numturbs = len(self.imported_params.index)

    def initialize_turbine_locations(self):

        ### Extract x,y position ###
        x = self.imported_params["x"].to_numpy()
        y = self.imported_params["y"].to_numpy()

        ### Calculate farm extents ###
        self.ex_x = [np.min(x), np.max(x)]
        self.ex_y = [np.min(y), np.max(y)]

        return np.array([x,y]).T        

    def setup_turbines(self):
        """
        Using the parameters and initial locations this function will populate the list of turbines
        """
        from windse.turbine_types import turbine_dict
        turbine_method = turbine_dict[self.turbine_type]
        for i,(x,y) in enumerate(self.initial_turbine_locations):
            df_row = self.imported_params.iloc[i]
            df_row = df_row.drop(labels=["x","y"])
            self.turbines.append(turbine_method(i,x,y,self.dom,imported_params=df_row))