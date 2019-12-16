from dolfin import Constant
# from windse import windse_parameters

# current_tab = windse_parameters.current_tab

def BaseHeight(x,y,ground,dx=0,dy=0):
    return Constant(ground(float(x),float(y),dx=dx,dy=dx))

def ConstantForce(x,turbines,dfd=None):
    ''' Returns a numpy vector corresponding to a scaler function space'''
    if dfd is None:
        return "Function Evaluation"
    elif dfd == "layout_x":
        return "Derivative with respect to x"
    elif dfd == "layout_y":
        return "Derivative with respect to y"
    elif dfd == "yaw":
        return "Derivative with respect to yaw"
    
    ...

    else:
    raise ValueError("Unknown derivative: "+dfd)

def DiskKernal(x,turbines,dfd=None):
    ''' Returns a numpy vector corresponding to a scaler function space'''
    if dfd is None:
        return "Function Evaluation"
    elif dfd == "layout_x":
        return "Derivative with respect to x"
    elif dfd == "layout_y":
        return "Derivative with respect to y"
    elif dfd == "yaw":
        return "Derivative with respect to yaw"
    
    ...

    else:
    raise ValueError("Unknown derivative: "+dfd)

def ActuatorDisk(x,turbines,force,space_kernal,dfd=None):
    ''' Returns a numpy vector corresponding to a vector function space '''
    if dfd is None:
        val = force(x,turbines)*space_kernal(x,turbines)
    else:
        val = force(x,turbines)*space_kernal(x,turbines,dfd=dfd) + force(x,turbines,dfd=dfd)*space_kernal(x,turbines)

    return convert_scaler_to_vector(val,yaw)

def TurbineForces(x,turbines,actuator,dfd=None):
    ''' Returns a numpy vector corresponding to a vector function space '''

    if dfd is None:
        tf1 = actuator(x,turbines,force,space_kernal)*cos(yaw)**2
        tf2 = actuator(x,turbines,force,space_kernal)*sin(yaw)**2
        tf3 = actuator(x,turbines,force,space_kernal)*cos(yaw)*sin(yaw)
    elif dfd == "yaw":
        tf1 = actuator(x,turbines,force,space_kernal,dfd=dfd)*cos(yaw)**2       + actuator(x,turbines,force,space_kernal)*(-2)*sin(yaw)     
        tf2 = actuator(x,turbines,force,space_kernal,dfd=dfd)*sin(yaw)**2       + actuator(x,turbines,force,space_kernal)*2*cos(yaw)     
        tf3 = actuator(x,turbines,force,space_kernal,dfd=dfd)*cos(yaw)*sin(yaw) + actuator(x,turbines,force,space_kernal)*(cos(2*yaw))
    else:
        tf1 = actuator(x,turbines,force,space_kernal,dfd=dfd)*cos(yaw)**2      
        tf2 = actuator(x,turbines,force,space_kernal,dfd=dfd)*sin(yaw)**2      
        tf3 = actuator(x,turbines,force,space_kernal,dfd=dfd)*cos(yaw)*sin(yaw)

    return [tf1,tf2,tf3]

def CalculateTurbineForces(x,turbines,actuator,force,space_kernal,dfd=None,df_index=None):
    if dfd is None:
        for i in range(turbines.n):
            tf1_temp,tf2_temp,tf3_temp = TurbineForces(x,turbines,ActuatorDisk,ConstantForce,DiskKernal,dfd=dfd)
            tf1 = tf1_temp
            tf2 = tf2_temp
            tf3 = tf3_temp
    else:
        tf1,tf2,tf3 = TurbineForces(x,turbines[df_index],ActuatorDisk,ConstantForce,DiskKernal,dfd=dfd)

    return numpy_to_function(tf1,tf2,tf3)
