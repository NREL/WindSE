"""
The ParameterManager controls the handles importing 
the parameters from the params.yaml file. These
functions don't need to be accessed by the end user.
"""

import __main__
import os
import yaml
import warnings
import copy


### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    import datetime
    import numpy as np
    from math import ceil
    import shutil
    import dolfin
    import sys
    import ast
    import difflib
    import inspect 

    # set_log_level(LogLevel.CRITICAL)

######################################################
### Collect all options and define general options ###
######################################################


### THis is a special class that allows prints to go to file and terminal
class Logger(object):
    def __init__(self, filename, std, rank):
        self.__dict__ = std.__dict__.copy() 
        self.terminal = std
        self.rank = rank
        if self.rank == 0:
            self.log = open(filename, "a")
            self.log.seek(0)
            self.log.truncate()

    def write(self, message):
        self.terminal.write(message)

        if self.rank == 0:
            self.log.write(message)  

    def flush(self):
        self.terminal.flush()

        if self.rank == 0:
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
        self.tagged_output = {}
        self.windse_path = os.path.dirname(os.path.realpath(__file__))
        self.defaults = yaml.load(open(self.windse_path+"/default_parameters.yaml"),Loader=yaml.SafeLoader)

        ### Update self with the defaults ###
        defaults_bak = copy.deepcopy(self.defaults)
        self.update(defaults_bak)

        ### Include the defaults from all the objectives ###
        import windse.objective_functions as obj_funcs
        self.obj_names = obj_funcs.objective_functions.keys()
        self["optimization"]["objective_type"] = obj_funcs.objective_kwargs

        ### Setup the defaults for the constraints based on objective functions###
        constraint_bak_types = self["optimization"]["constraint_types"]
        for key, value in obj_funcs.objective_kwargs.items():
            constraint_bak_types[key] = {
                "target": None,
                "scale": 1.0,
                "kwargs": value
            }

        # print(dir(obj_funcs))
        # print(obj_funcs.alm_power())
        # exit()

    def TerminalUpdate(self,dic,keys,value):
        if len(keys) > 1:
            next_dic = dic.setdefault(keys[0],{})
            self.TerminalUpdate(next_dic,keys[1:],value)
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

    def CheckParameters(self,updates,defaults,out_string=""):
        default_keys = defaults.keys()
        for key in updates.keys():
            split_key = key.split("_#")[0]
            if split_key not in default_keys:
                suggestion = difflib.get_close_matches(key, default_keys, n=1)
                if suggestion:
                    raise KeyError(out_string + key + " is not a valid parameter, did you mean: "+suggestion[0])
                else:
                    raise KeyError(out_string + key + " is not a valid parameter")
            elif isinstance(updates[key],dict):
                in_string =out_string + key + ":"
                self.CheckParameters(updates[key],defaults[split_key],out_string=in_string)

    def RecordUserSupplied(self,yaml_file,defaults):
        user_supplied = {}

        if not isinstance(defaults,dict):
            return True

        for key in defaults.keys():
            user_supplied[key] = False


        for key in yaml_file.keys():
            split_key = key.split("_#")[0]
            if isinstance(yaml_file[key],dict):
                sub_supplied = self.RecordUserSupplied(yaml_file[key],defaults[split_key])
                user_supplied[split_key] = sub_supplied
            else:
                user_supplied[split_key] = True
    
        return user_supplied

    def NestedUpdate(self,dic,subdic=None):
        if subdic is None:
            target_dic = self
        else:
            target_dic = subdic

        for key, value in dic.items():
            if isinstance(value,dict):
                target_dic[key] = self.NestedUpdate(value,subdic=target_dic[key])
            else:
                target_dic[key] = value
        return target_dic

    def Load(self, loc,updated_parameters=[]):
        """
        This function loads the parameters from the .yaml file. 
        It should only be assessed once from the :meth:`windse.initialize` function.

        Args:
            loc (str): This string is the location of the .yaml parameters file.

        """

        # Create an MPI communicator and initialize rank and num_procs 
        self.comm = dolfin.MPI.comm_world
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        ### Load the yaml file (requires PyYaml)
        if isinstance(loc,dict):
            self.fprint("Loading from dictionary")
            yaml_file = loc
        else:
            self.fprint("Loading: "+loc)
            yaml_file = yaml.load(open(loc),Loader=yaml.SafeLoader)

        ### update any parameters if supplied ###
        for p in updated_parameters:
            keys_list = p.split(":")
            self.TerminalUpdate(yaml_file,keys_list[:-1],keys_list[-1])

        ### Check for incorrect parameters ###
        self.CheckParameters(yaml_file,self)
        self.fprint("Parameter Check Passed")

        ### record which parameters were set by the user ###
        self.user_supplied = self.RecordUserSupplied(yaml_file,self.defaults)

        ### Setup objective functions if needed ###
        yaml_op = yaml_file.get("optimization",{})
        objective_type = yaml_op.pop("objective_type", None)
        constraint_types = yaml_op.pop("constraint_types", None)

        ### Load in the defaults objective dictionaries 
        import windse.objective_functions as obj_funcs

        ### Replace the dictionary defaults with the real default
        if objective_type is None:
            objective_type = self.defaults["optimization"]["objective_type"]
        if constraint_types is None:
            constraint_types = self.defaults["optimization"]["constraint_types"]

        ### Process the objective keyword arguments
        if isinstance(objective_type,str):
            objective_type = {objective_type: obj_funcs.objective_kwargs[objective_type]}
        elif isinstance(objective_type,list):
            new_objective_type = {}
            for obj in objective_type:
                new_objective_type[obj] = obj_funcs.objective_kwargs[obj]
            objective_type = new_objective_type
        elif isinstance(objective_type,dict):
            ### make sure to add in any default values the user may not have set for the objectives 
            for key, value in objective_type.items():
                objective_split = key.split("_#")[0]
                obj_default = obj_funcs.objective_kwargs[objective_split]
                for k, v in obj_default.items():
                    if k not in value.keys():
                        value[k] = v

        ### Process the constraints dictionary ###
        for key, value in constraint_types.items():
            if "target" not in value.keys():
                raise ValueError(f"A target needs to be defined for the {key} constraint")
            if "scale" not in value.keys():
                value["scale"] = 1.0

            ### check if the objective function keywords were supplied
            if key != "min_dist":
                constraints_split = key.split("_#")[0]
                constraints_kw_default = obj_funcs.objective_kwargs[constraints_split]
                constraints_kw = value.get('kwargs',{})
                for k, v in constraints_kw_default.items():
                    if k not in constraints_kw.keys():
                        constraints_kw[k] = v
                value["kwargs"] = constraints_kw

        ### Set the parameters ###
        self.update(self.NestedUpdate(yaml_file))
        self["optimization"]["objective_type"] = objective_type
        self["optimization"]["constraint_types"] = constraint_types

        ### Create Instances of the general options ###
        for key, value in self["general"].items():
            setattr(self,key,value)

        ### Check if dolfin_adjoint is unnecessary or required ###
        opt_gradient = yaml_file.get("optimization",{}).get("gradient",False)
        opt_taylor   = yaml_file.get("optimization",{}).get("taylor_test",False)
        opt_optimize = yaml_file.get("optimization",{}).get("optimize",False)
        self.performing_opt_calc = opt_gradient or opt_taylor or opt_optimize
        #if self.performing_opt_calc and not self.dolfin_adjoint:
        #    raise ValueError("Asked to perform gradient, Taylor test, or optimization but general:dolfin_adjoint is set to False. These operations will not work without dolfin_adjoint.")
        #elif not self.performing_opt_calc and self.dolfin_adjoint: 
        #    warnings.warn("general:dolfin_adjoint is set to True but no optimization parameters provided. This will cause unneeded overhead.")

        # print(self.dolfin_adjoint)
        # for module in sys.modules:
        #     if "adjoint" in module:
        #         print(module)

        ### set default name ###
        if self.name is None:
            _, yaml_name = os.path.split(loc)
            self.name = yaml_name.split(".")[0]

        ### Set up the folder Structure ###
        timestamp=datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        fancytimestamp=datetime.datetime.today().strftime('%Y/%m/%d_%H:%M:%S')
        if self.preappend_datetime:
            self.name = timestamp+"-"+self.name
            self["general"]["name"]=self.name
        self.folder = self.output_folder+self.name+"/"
        self["general"]["folder"] = self.folder

        # Create all needed directories ahead of time
        if self.rank == 0:
            # Try to create the parent folder
            os.makedirs(self.folder, exist_ok=True)

            # Try to create all sub folders within this parent
            subfolder_list = ['data',
                              'functions',
                              'input_files',
                              'mesh',
                              'plots',
                              'timeSeries',
                              'profiling']

            for sub in subfolder_list:
                os.makedirs('%s/%s' % (self.folder, sub), exist_ok=True)

            os.makedirs('%s/data/alm' % (self.folder), exist_ok=True)
            os.makedirs('%s/data/alm/rotor_force' % (self.folder), exist_ok=True)
            os.makedirs('%s/data/alm/angle_of_attack' % (self.folder), exist_ok=True)

        # Wait until rank 0 has created the directory structure
        self.comm.barrier()

        # ### Make sure folder exists ###
        # comm = MPI.comm_world
        # rank = comm.Get_rank()
        # if not os.path.exists(self.folder) and rank == 0: os.makedirs(self.folder)
        # if not os.path.exists(self.folder+"input_files/") and rank == 0: os.makedirs(self.folder+"input_files/")
        
        ### Setup the logger ###
        self.log = self.folder+"log.txt"
        sys.stdout = Logger(self.log, sys.stdout, self.rank)
        sys.stderr = Logger(self.log, sys.stderr, self.rank)

        ### Copy params file to output folder ###
        if isinstance(loc,str):
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
            if not os.path.exists(self.folder+subfolder) and self.rank == 0: os.makedirs(self.folder+subfolder)

            if filetype == "pvd":
                file_string = self.folder+subfolder+filename+".pvd"
                out = dolfin.File(file_string)
                out << (func,val)
            elif filetype == "xdmf":
                file_string = self.folder+subfolder+filename+".xdmf"
                out = dolfin.XDMFFile(file_string)
                out.write(func,val)

            func.rename(old_filename,old_filename)
            return out

        else:
            if filetype == "pvd" or isinstance(func,type(dolfin.Mesh)):
                file << (func,val)
            elif filetype == "xdmf":
                file.write(func,val)

            func.rename(old_filename,old_filename)
            return file

    def save_csv(self, filename, data=None, subfolder="", header=None, mode='w'):
        ### Check Processor ###
        if self.rank == 0:

            ### Set the output folder ###
            out_f = subfolder

            ### Check if folder exists ###
            if not os.path.exists(out_f): os.makedirs(out_f, exist_ok=True)

            ### Open the file ###
            f = open(out_f+filename+".csv",mode)

            ### Write the header if present ###
            if header is not None:
                f.write(header)
                f.write("\n")

            ### Write the data ###
            if data is not None:
                np.savetxt(f,data, delimiter=', ')

            ### Close the file ###
            f.close()
        dolfin.MPI.comm_world.barrier()

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
        if self.rank == 0:
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

    def pprint(self):
        '''
        This function print the contents of the params dict in a more readable format
        '''
        dummy = {}
        for key, value in self.items():
            dummy[key] = value

        if self.rank == 0:
            print(yaml.dump(dummy))

    def tag_output(self, key, value, collective_output=None):

        # self.fprint("Tagging Debug item: "+key)
        ### Process value ###
        if not isinstance(value,int):
            value = float(value)

        if self.num_procs > 1:
            send_data = np.float64(value)
            mpi_buff = np.zeros(self.num_procs, dtype=np.float64)
            self.comm.Gather(send_data, mpi_buff, root=0)

            differing_opinions = False

            if not np.all(np.isclose(mpi_buff, mpi_buff[0])):
                differing_opinions = True

            if differing_opinions or collective_output is not None:
                if collective_output == 'sum' or 'sum' in key:
                    value = np.sum(mpi_buff)
                elif collective_output == 'avg' or 'avg' in key:
                    value = np.mean(mpi_buff)
                elif collective_output == 'max' or 'max' in key:
                    value = np.amax(mpi_buff)
                elif collective_output == 'min' or 'min' in key:
                    value = np.amin(mpi_buff)
                else:
                    print('WARNING: tagging %s in parallel may result in disagreement between processors.' % (key))

                value = float(value)

        if self.rank == 0:
            ### Grab the name of the module that called this function ###
            stack = inspect.stack()[1][0]
            mod = inspect.getmodule(stack)
            the_module = mod.__name__.split(".")[-1]
            the_class = stack.f_locals["self"].__class__.__name__
            the_method = stack.f_code.co_name

            ### This will tell exactly where this function was called from ###
            # print("I was called by {}:{}.{}()".format(the_module, the_class, the_method))

            ### Check if that module has called before and add the dictionary entries ###
            if the_module in self.tagged_output.keys():
                self.tagged_output[the_module].update({key: value})
            else:
                self.tagged_output.update({the_module: {key: value}})

            ### Update the yaml file ###
            with open(self.folder+"tagged_output.yaml","w") as file:
                yaml.dump(self.tagged_output, file, sort_keys=False)

            ### Print the new dict ###
            # print(self.tagged_output)

windse_parameters = Parameters()
