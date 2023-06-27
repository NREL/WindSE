import time
import numpy as np

def CreateRefinementList(dom, farm, refine_params):  
    farm_num          = refine_params["farm_num"]       
    farm_type         = refine_params["farm_type"]      
    farm_factor       = refine_params["farm_factor"]    
    turbine_num       = refine_params["turbine_num"]    
    turbine_type      = refine_params["turbine_type"]    
    turbine_factor    = refine_params["turbine_factor"] 
    refine_custom     = refine_params["refine_custom"]  ### Need to fix for if the domain is scaled
    refine_power_calc = refine_params["refine_power_calc"]

    refine_list = []

    RDs = farm.get_rotor_diameters()
    zs  = farm.get_hub_locations()[:,2]

    if refine_custom is not None:
        refine_list = refine_list+refine_custom

    if farm_num > 0:
        bbox = farm.calculate_farm_bounding_box()

        for i in range(farm_num,0,-1):
            # expand_factor = 1+(farm_factor-1)*(i)
            expand_factor = (farm_factor)**(i)

            if farm_type == 'box':
                RD = max(RDs)
                bbox[2] = [bbox[2][0]-RD,bbox[2][1]+RD]
                # refine_list.append(["box",[bbox,expand_factor]])
                refine_list.append({"type": "box",
                                    "x_range": bbox[0],
                                    "y_range": bbox[1],
                                    "z_range": bbox[2]})

            elif farm_type == 'cylinder':
                RD = max(RDs)
                x0 = (bbox[0][1]+bbox[0][0])/2.0
                y0 = (bbox[1][1]+bbox[1][0])/2.0
                if dom.dim == 3:
                    z0 = bbox[2][0]-RD
                    center = [x0,y0,z0]
                    height = bbox[2][1]-bbox[2][0]+2*RD
                else:
                    center = [x0,y0]
                    height = 0
                radius = np.sqrt((bbox[0][1]-bbox[0][0])**2+(bbox[1][1]-bbox[1][0])**2)/2.0
                # refine_list.append(["cylinder",[center,radius,height,expand_factor]])
                refine_list.append({"type": "cylinder",
                                    "center": center,
                                    "radius": radius,
                                    "height": height,
                                    "expand_factor": expand_factor})

            elif farm_type == 'stream':
                RD = max(RDs)
                x0 = bbox[0][0]-RD
                y0 = (bbox[1][1]+bbox[1][0])/2.0
                if dom.dim == 3:
                    z0 = (min(zs)+max(zs))/2.0
                    center = [x0,y0,z0]
                    radius = np.sqrt((bbox[1][1]-bbox[1][0])**2+(bbox[2][1]-bbox[2][0])**2)/2.0
                else:
                    center = [x0,y0]
                    radius = (bbox[1][1]-bbox[1][0])/2.0
                length = bbox[0][1]-bbox[0][0]+6*RD
                theta = dom.inflow_angle
                pivot_offset = 3*max(RDs)/2.0
                # refine_list.append(["stream",[center,radius,length,theta,pivot_offset,expand_factor]])
                refine_list.append({"type": "stream",
                                    "center": center,
                                    "radius": radius,
                                    "length": length,
                                    "theta": theta,
                                    "pivot_offset": pivot_offset,
                                    "expand_factor": expand_factor})

    if turbine_num > 0 and farm.numturbs > 0:
        for i in range(turbine_num,0,-1):
            expand_factor = (turbine_factor)**(i)

            if turbine_type == 'simple':
                radius = max(RDs)
                # refine_list.append(["simple",[radius,expand_factor]])
                refine_list.append({"type": "simple",
                                    "radius": radius,
                                    "expand_factor": expand_factor})

            elif turbine_type == 'sphere':
                radius = max(RDs)
                # refine_list.append(["sphere",[radius,expand_factor]])
                refine_list.append({"type": "sphere",
                                    "radius": radius,
                                    "expand_factor": expand_factor})

            elif turbine_type == 'wake':
                radius = max(RDs)
                length = 5*radius
                theta = dom.inflow_angle
                # refine_list.append(["wake",[radius,length,theta,expand_factor]])
                refine_list.append({"type": "wake",
                                    "radius": radius,
                                    "length": length,
                                    "theta": theta,
                                    "expand_factor": expand_factor})

            elif turbine_type == 'tear':
                radius = max(RDs)
                theta = dom.inflow_angle
                # refine_list.append(["tear",[radius,theta,expand_factor]])
                refine_list.append({"type": "tear",
                                    "radius": radius,
                                    "theta": theta,
                                    "expand_factor": expand_factor})

        if refine_power_calc:
            radius = max(RDs)
            length = radius/5.0
            theta = dom.inflow_angle
            centered = True
            # refine_list.append(["wake",[radius,length,theta,expand_factor,centered]])
            refine_list.append({"type": "wake",
                                "radius": radius,
                                "length": length,
                                "theta": theta,
                                "expand_factor": expand_factor,
                                "centered": centered})

    if refine_custom is not None:
        sorted_refine_custom_keys = sorted(list(refine_custom.keys()))

        for key in sorted_refine_custom_keys:
            refine_list.append(refine_custom[key])

    return refine_list

def RefineMesh(dom,farm):

    ### Define the possible operations ###
    refine_dict = {"full": dom.Refine,
                   "box": dom.BoxRefine,
                   "cylinder": dom.CylinderRefine,
                   "stream": dom.StreamRefine,
                   "simple": farm.SimpleRefine,
                   "sphere": farm.SphereRefine,
                   "tear": farm.TearRefine,
                   "wake": farm.WakeRefine
                   }

    ### Alias the print command ###
    fprint = dom.params.fprint

    ### Convert Parameters to a list of refine instructions ###
    refine_params = dom.params["refine"]
    refine_list = CreateRefinementList(dom,farm,refine_params)

    ### Step through refine instructions ###
    num = len(refine_list)
    for i, refinement_dict in enumerate(refine_list):
        fprint("Refining Mesh Step {:d} of {:d}".format(i+1,num), special="header")
        step_start = time.time()
        
        refine_type = refinement_dict["type"]
        refine_args = dict(refinement_dict)
        del refine_args["type"]

        refine_func = refine_dict[refine_type]
        refine_func(**refine_args)

        step_stop = time.time()
        fprint("Step {:d} of {:d} Finished: {:1.2f} s".format(i+1,num,step_stop-step_start), special="footer")

def WarpMesh(dom):

    warp_type      = dom.params["refine"]["warp_type"]      
    warp_strength  = dom.params["refine"]["warp_strength"]  
    warp_height    = dom.params["refine"]["warp_height"]    
    warp_percent   = dom.params["refine"]["warp_percent"]  

    if warp_type == "smooth":
        dom.WarpSmooth(warp_strength)
    elif warp_type == "split":
        dom.WarpSplit(warp_height*dom.xscale,warp_percent)