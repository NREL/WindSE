"""
The ParameterManager controls the handles importing 
the parameters from the params.yaml file. These
functions don't need to be accessed by the end user.
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    import yaml
    import datetime
    import numpy as np
    from math import ceil
    import shutil
    from dolfin import *
    import sys
    import ast

    # set_log_level(LogLevel.CRITICAL)

######################################################
### Collect all options and define general options ###
######################################################


### THis is a special class that allows prints to go to file and terminal
class Logger(object):
    def __init__(self,filename):
        self.__dict__ = sys.stdout.__dict__.copy() 
        self.terminal = sys.stdout
        self.log = open(filename, "a")
        self.log.seek(0)
        self.log.truncate()

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        self.terminal.flush()
        self.log.flush()
        pass   

    def isatty(self):
        return self.terminal.isatty() 

class Parameters(dict):
    """
    Parameters is a subclass of pythons *dict* that adds
    function specific to windse.
    """
    def __init__(self):
        super(Parameters, self).__init__()
        self.current_tab = 0

    def NestedUpdate(self,dic,keys,value):
        if len(keys) > 1:
            next_dic = dic.setdefault(keys[0],{})
            self.NestedUpdate(next_dic,keys[1:],value)
        elif len(keys) == 1:
            current_value = dic.get(keys[0],"")
            if isinstance(current_value,int):
                dic[keys[0]] = int(value)
            elif isinstance(current_value,float):
                dic[keys[0]] = float(value)
            elif isinstance(current_value,str):
                dic[keys[0]] = value
            elif isinstance(current_value,list):
                dic[keys[0]] = ast.literal_eval(value)


    def Load(self, loc,updated_parameters=[]):
        """
        This function loads the parameters from the .yaml file. 
        It should only be assessed once from the :meth:`windse.initialize` function.

        Args:
            loc (str): This string is the location of the .yaml parameters file.

        """

        ### Load the yaml file (requires PyYaml)
        yaml_file = yaml.load(open(loc),Loader=yaml.SafeLoader)

        ### update any parameters if supplied ###
        for p in updated_parameters:
            keys_list = p.split(":")
            self.NestedUpdate(yaml_file,keys_list[:-1],keys_list[-1])

        ### Set the parameters
        self.update(yaml_file)

        ### Create Instances of the general options ###
        self.name = self["general"].get("name", "Test")
        self.preappend_datetime = self["general"].get("preappend_datetime", False)
        self.output_type = self["general"].get("output_type", "pvd")
        self.dolfin_adjoint = self["general"].get("dolfin_adjoint", False)
        self.output = self["general"].get("output", ["solution"])

        ### Print some stats ###

        ### Set up the folder Structure ###
        timestamp=datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        fancytimestamp=datetime.datetime.today().strftime('%Y/%m/%d_%H:%M:%S')
        if self.preappend_datetime:
            self.name = timestamp+"-"+self.name
            self["general"]["name"]=self.name
        self.folder = "output/"+self.name+"/"
        self["general"]["folder"] = self.folder

        ### Make sure folder exists ###
        if not os.path.exists(self.folder): os.makedirs(self.folder)
        if not os.path.exists(self.folder+"input_files/"): os.makedirs(self.folder+"input_files/")
        
        ### Setup the logger ###
        self.log = self.folder+"log.txt"
        sys.stdout = Logger(self.log)

        ### Copy params file to output folder ###
        shutil.copy(loc,self.folder+"input_files/")

        ### Create checkpoint if required ###
        # if self.save_file_type == "hdf5":
        #     self.Hdf=HDF5File(MPI.mpi_comm(), self.folder+"checkpoint/checkpoint.h5", "w")

        ### Print some more stuff
        self.fprint("General Parameter Information", special="header")
        self.fprint("Run Name: {0}".format(self.name))
        self.fprint("Run Time Stamp: {0}".format(fancytimestamp))
        self.fprint("Output Folder: {0}".format(self.folder))
        if updated_parameters:
            self.fprint("Updated Parameter:")
            for i,p in enumerate(updated_parameters):
                self.fprint("{:d}: {:}".format(i,p),offset=1)
        self.fprint("Parameters Setup", special="footer")

    def Read(self):
        """
        This function reads the current state of the parameters object 
        and prints it in a easy to read way.
        """
        for group in self:
            print(group)
            max_length = 0
            for key in self[group]:
                max_length = max(max_length,len(key))
            max_length = max_length
            for key in self[group]:
                print("    "+key+":  "+" "*(max_length-len(key))+repr(self[group][key]))

    def Save(self, func, filename, subfolder="",val=0,file=None,filetype="default"):
        """
        This function is used to save the various dolfin.Functions created
        by windse. It should only be accessed internally.

        Args:
            func (dolfin.Function): The Function to be saved
            filename (str): the name of the function

        :Keyword Arguments:
            * **subfolder** (*str*): where to save the files within the output folder
            * **n** (*float*): used for saving a series of output. Use n=0 for the first save.

        """
        self.fprint("Saving: {0}".format(filename))

        # if not isinstance(init_func,Function):
        #     func = Function(func)
        # else:
        #     func = init_func

        ### Name the function in the meta data, This should probably be done at creation
        old_filename = func.name()
        func.rename(filename,filename)

        if filetype == "default":
            filetype = self.output_type

        if file is None:
            ### Make sure the folder exists
            if not os.path.exists(self.folder+subfolder): os.makedirs(self.folder+subfolder)

            if filetype == "pvd":
                file_string = self.folder+subfolder+filename+".pvd"
                out = File(file_string)
                out << (func,val)
            elif filetype == "xdmf":
                file_string = self.folder+subfolder+filename+".xdmf"
                out = XDMFFile(file_string)
                out.write(func,val)

            func.rename(old_filename,old_filename)
            return out

        else:
            if filetype == "pvd" or isinstance(func,type(Mesh)):
                file << (func,val)
            elif filetype == "xdmf":
                file.write(func,val)

            func.rename(old_filename,old_filename)
            return file

    def fprint(self,string,tab=None,offset=0,special=None):
        """
        This is just a fancy print function that will tab according to where
        we are in the solve

        Args:
            string (str): the string for printing

        :Keyword Arguments:
            * **tab** (*int*): the tab level

        """
        ### Check Processor ###
        rank = 0
        if rank == 0:
            ### Check if tab length has been overridden
            if tab is None:
                tab = self.current_tab
            
            ### Check if we are starting or ending a section
            if special=="header":
                self.current_tab += 1
                self.fprint("",tab=tab)
            elif special =="footer":
                self.current_tab -= 1
                tab -= 1
                self.fprint("",tab=tab+1)

            ### Apply Offset if provided ###
            tab += offset

            ### Create Tabbed string ###
            tabbed = "|    "*tab

            ### Apply Tabbed string ###
            if isinstance(string,str):
                string = tabbed+string
            else:
                string = tabbed+repr(string)

            ### Print ###
            # print(string, flush=True)
            print(string)
            sys.stdout.flush()

            if special=="header":
                self.fprint("",tab=tab+1)

windse_parameters = Parameters()