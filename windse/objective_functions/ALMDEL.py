################## HEADER DO NOT EDIT ##################
import os
import __main__
import fatpack
import pandas as pd
import numpy as np

def DEL(df1):
   slope = 10 #Wohler Exponent
   elapsed =  df1.Time.max() - df1.Time.min()
   ts = df1.Moment
   ranges = fatpack.find_rainflow_ranges(ts)
   Nrf, Srf = fatpack.find_range_count(ranges, 100)
   DELs = Srf ** slope * Nrf / elapsed
   DEL = DELs.sum() ** (1 / slope)
   return(DEL)

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    from dolfin import *
    from dolfin_adjoint import *
########################################################

### Additional import statements ###
import numpy as np
import math
import os

### Declare Unique name
name = "alm_DELs"

### Set default keyword argument values ###
keyword_defaults = {
    "alm_DEL_type": "flapwise",
    'DEL_start_time': 100
    }

### Define objective function
def objective(solver, inflow_angle=0.0, first_call=False, **kwargs):
    '''
    The "alm_DELs" objective function computes the damage equivalen loading (DEL) using actuator lines 
    by applying the rainflow counting algorithm to the moment arm of the turbine. 
    Can be used for multiple turbines.

    Keyword arguments:
        alm_DEL_type: flapwise 
                            flapwise - out-of-plane bending moment
    '''

    # check if we are attempting to optimize the DEL
    if solver.optimize == True:
        raise(Exception("Can't optimize DEL yet!"))


    ### Extract keyword arguments
    alm_DEL_type = kwargs.pop("alm_DEL_type")
    startTime = kwargs.pop("DEL_start_time")
    DELS = []

    # get current time
    time = solver.problem.simTime_list[-1]

    # check if we are at least 10 s past the burn-in time
    if time < startTime + 10:
       return(np.nan)

    # compute DEL for each turbine
    for turb_i in range(solver.problem.farm.numturbs):

        # use fx to compute flapwise bending
        dat = pd.read_csv(solver.problem.force_files[turb_i][0], sep=', ')

        # drop MPI duplicates
        dat = dat[~np.isnan(dat['r0_n000'])]
        dat.drop_duplicates(subset='time', inplace=True, ignore_index=True, keep='last')
        dat.time = dat.time.astype(float)

        # identify number of nodes
        cols = dat.drop('time', 1).columns
        cols = cols[['r0' in col for col in cols]]
        bladePoss = np.linspace(0, 1, len(cols))

        # create moment nodes
        RD = solver.problem.farm.turbines[turb_i].RD
        Radius = RD / 2.
        bladeLocs = Radius *  bladePoss
        moments = np.sum(bladeLocs * dat[cols], axis=1)
        df = pd.DataFrame({'Time': dat.time, 'Moment': moments})

        DELS.append(DEL(df[df.Time >= startTime]))

    J = np.array([np.sum(DELS)])

    #J_list=np.zeros(solver.problem.farm.numturbs+2)
    #J_list[0]=solver.simTime
    #J = 0.0
    #for i in range(solver.problem.farm.numturbs):
    #    if alm_power_type == "flapwise":
    #        J_temp = assemble(1e-6*(2.0*np.pi*solver.problem.rpm/60.0)*inner(-solver.problem.tf_list[i], solver.problem.cyld_expr_list[i])*dx)
    #    else:
    #        raise ValueError("Unknown ALM Power type: " + repr(alm_power_type))
    #    J_list[i+1] = J_temp
    #    J += J_temp
    #J_list[-1] = float(J)


    #if solver.save_DEL or solver.save_objective:

    #    folder_string = solver.params.folder+"data/"
    #    if not os.path.exists(folder_string): os.makedirs(folder_string)

    #    if first_call:
    #        f = open(folder_string+"alm_power_data.txt",'w')
    #        header = str("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
    #        f.write(header)
    #    else:
    #        f = open(folder_string+"alm_power_data.txt",'a')

    #    np.savetxt(f,[J_list])
    #    f.close()

    return(J)
