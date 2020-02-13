import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *



def CalculatePowerFunctional(solver,inflow_angle = 0.0):
    J = -assemble(dot(solver.problem.tf,solver.u_next)*dx)
    return J

def CalculateWakeCenter(solver,inflow_angle = 0.0):

    if not hasattr(solver,"outflow_markers"):
        y_factor = 2.0
        z_factor = 1.2
        HH = solver.problem.farm.HH[0]
        x0 = solver.problem.dom.x_range[1]
        y0 = solver.problem.farm.y[0]
        R = solver.problem.farm.RD[0]/2.0
        ly=y0-y_factor*R
        uy=y0+y_factor*R
        lz=HH-z_factor*R
        uz=HH+z_factor*R
        outregion  = CompiledSubDomain("near(x[0], x0, tol) && x[1]>=ly && x[1]<=uy  && x[2]>=lz && x[2]<=uz && on_boundary",x0 = x0, ly=ly, uy=uy, lz=lz, uz=uz, tol = 1e-10)
        solver.outflow_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim() - 1)
        solver.outflow_markers.set_all(0)
        outregion.mark(solver.outflow_markers,1)

        # File("test.pvd") << solver.outflow_markers




    ds = Measure('ds', subdomain_data=solver.outflow_markers)
    x = SpatialCoordinate(solver.problem.dom.mesh)

    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    # u_dif_mag = sqrt((u[0]-u_ref[0])**2.0+(u[1]-u_ref[1])**2.0+(u[2]-u_ref[2])**2.0)
    u_dif_mag = sqrt((u[0]-u_ref[0])**2.0)

    M = assemble(u_dif_mag*ds(1))
    Mx = assemble(x[0]*u_dif_mag*ds(1))
    My = assemble(x[1]*u_dif_mag*ds(1))
    Mz = assemble(x[2]*u_dif_mag*ds(1))

    print(Mx/M,My/M,Mz/M)
    return My/M

objectives_dict = {"power":    CalculatePowerFunctional,
                   "wake_deflection": CalculateWakeCenter
                  }