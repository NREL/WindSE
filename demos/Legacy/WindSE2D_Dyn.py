from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
from scipy import integrate


set_log_level(LogLevel.INFO)


T = 500.0            # final time
num_steps = 1000  # number of time steps
dt = T / num_steps # time step size
mu = 16        # dynamic viscosity
rho = 1            # density
save_int = 5

inflowVel=8

# mesh parameters
degree = 2
Lx = 5000
Ly = 5000

nx = 72
ny = 72

RD=126

#WTG parameters
numturbs = 9

#number of inflow direction bins
bins = 1
WTGexp = 6.
radius = RD/2.
thickness = RD/10.
numRefine = 1
A=RD # weird for 2D
HH=80    
initExtent=1.
mlDenom=5

restart = False
randStart = False
gridStart = True
optimize = False
loadnuT = False

mesh = RectangleMesh(Point(-Lx/2., -Ly/2.), Point(Lx/2., Ly/2.), nx, ny)

site_x = 1000
site_y = 1000
refine_x = 1100
refine_y = 1100

def refine_mesh(mesh, refine_x, refine_y):
    #refines the mesh around the site boundaries
    h = mesh.hmin()
    
    cell_markers = MeshFunction('bool',mesh, mesh.topology().dim())
    cell_markers.set_all(False)
    for cell in cells(mesh):
        if (cell.midpoint()[0] > -(refine_x)) and (abs(cell.midpoint()[1]) < refine_y ):
            cell_markers[cell] = True

    mesh = refine(mesh, cell_markers)

    return mesh

def refine_mesh2(mesh, refine_x, refine_y):
    #refines the mesh around the site boundaries
    h = mesh.hmin()
    
    cell_markers = MeshFunction('bool',mesh, mesh.topology().dim())
    cell_markers.set_all(False)
    for cell in cells(mesh):
        if (cell.midpoint()[0]**2 + cell.midpoint()[1]**2 < refine_x**2+refine_y**2 ):
            cell_markers[cell] = True

    mesh = refine(mesh, cell_markers)

    return mesh

for nums in range(numRefine):
    print('refining mesh')
    mesh=refine_mesh2(mesh, refine_x, refine_y)
    h = mesh.hmin()

Re = Lx*8/mu
print(Re)

print(mesh.hmin())
print(inflowVel*dt/mesh.hmin())

alpha = 7*pi/64

# Define function spaces
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

print(V.dim())

def WTGdist(x,y):
    return np.exp(-((x/thickness)**WTGexp + (y/radius)**WTGexp))

def createLayout(numturbs):
    mx=[]
    my=[]
    mz=[]

    if randStart == True:
        for i in range(numturbs):
            mx.append(Constant(np.random.uniform(low=-(site_x - radius),high=(site_x - radius))))
            my.append(Constant(np.random.uniform(low=-(site_y - radius), high=(site_y - radius))))
            mz.append(Constant(HH)) 

    elif gridStart ==True:

        if numturbs == 16:
            rows = 4
            cols = 4
            xpos = np.linspace(-initExtent*(site_x - radius),initExtent*(site_x - radius),cols)
            ypos = np.linspace(-initExtent*(site_y - radius),initExtent*(site_y - radius),rows)
            for i in range(rows):
                for j in range(cols):
                    mx.append(Constant(xpos[j]))
                    my.append(Constant(ypos[i]))
                    # # some starting noise sometimes helps
                    # mx.append(Constant(xpos[j]+5.*np.random.randn()))
                    # my.append(Constant(ypos[i]+5.*np.random.randn()))
                    mz.append(Constant(HH))

        if numturbs == 9:
            rows = 3
            cols = 3
            xpos = np.linspace(-site_x,site_x,cols)
            ypos = np.linspace(-site_y,site_y,rows)
            for i in range(rows):
                for j in range(cols):
                    mx.append(Constant(xpos[j]))
                    my.append(Constant(ypos[i]))
                    # # some starting noise sometimes helps
                    # mx.append(Constant(xpos[j]+5.*np.random.randn()))
                    # my.append(Constant(ypos[i]+5.*np.random.randn()))
                    mz.append(Constant(HH))

        if numturbs == 1:
            mx.append(Constant(-1500))
            my.append(Constant(0))
            mz.append(Constant(HH))

        if numturbs == 2:
            mx.append(Constant(-1500))
            mx.append(Constant(-1500 + 7*RD))
            my.append(Constant(0))
            my.append(Constant(0))
            mz.append(Constant(HH))
            mz.append(Constant(HH))

        if numturbs == 3:
            mx.append(Constant(-1000))
            mx.append(Constant(0))
            mx.append(Constant(1000))
            my.append(Constant(0))
            my.append(Constant(0))
            my.append(Constant(0))
            mz.append(Constant(HH))
            mz.append(Constant(HH))
            mz.append(Constant(HH))

        if numturbs == 4:
            mx.append(Constant(-1200))
            mx.append(Constant(-400))
            mx.append(Constant(400))
            mx.append(Constant(1200))
            my.append(Constant(0))
            my.append(Constant(0))
            my.append(Constant(0))
            my.append(Constant(0))
            mz.append(Constant(HH))
            mz.append(Constant(HH))
            mz.append(Constant(HH))
            mz.append(Constant(HH))

    return mx, my, mz

def createRotatedTurbineForce(mx,my,ma,A,beta,numturbs,alpha,V):
    x=SpatialCoordinate(mesh)
    tf = Function(V)

    for i in range(numturbs):
        WTGbase = project(Expression(("cos(yaw)","-sin(yaw)"),yaw=myaw[i],degree=2),V)
        # WTGbase = project(Expression(("0","1"),yaw=myaw[i],degree=2),V)


        #rotation
        mxrot = cos(alpha)*mx[i] - sin(alpha)*my[i]
        myrot = sin(alpha)*mx[i] + cos(alpha)*my[i]
        # mxrot=mx[i]
        # myrot=my[i]

        x_centered=x[0]-mxrot
        y_centered=x[1]-myrot

        x_centered_rotated = x_centered*cos(myaw[i]) + y_centered*sin(myaw[i])
        y_centered_rotated = -x_centered*sin(myaw[i]) + y_centered*cos(myaw[i])

        # tf = tf+ 0.0001*exp(-(((x[0] - mx[i])/thickness)**WTGexp +(((x[1] - my[i])**2)/radius**2)**WTGexp))*WTGbase
        # tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-((x_centered/thickness)**WTGexp + ((y_centered-radius/2.)/(radius/2.))**WTGexp))*WTGbase
        # tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-((x_centered/thickness)**WTGexp + ((y_centered+radius/2.)/(radius/2.))**WTGexp))*WTGbase
        tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-((x_centered_rotated/thickness)**WTGexp + ((y_centered_rotated-radius/2.)/(radius/2.))**WTGexp))*WTGbase
        tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-((x_centered_rotated/thickness)**WTGexp + ((y_centered_rotated+radius/2.)/(radius/2.))**WTGexp))*WTGbase

        # tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-(((x[0]*cos(myaw[i]) - x - mxrot)/thickness)**WTGexp + ((x[1] - myrot-radius/2.)/(radius/2.))**WTGexp))*WTGbase
        # tf = tf + 0.5*4.*A*ma[i]/(1.-ma[i])/beta*exp(-(((x[0]*cos(myaw[i]) - mxrot)/thickness)**WTGexp + ((x[1] - myrot+radius/2.)/(radius/2.))**WTGexp))*WTGbase

    return tf

#boundary conditions
class walls(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1]**2 - (Ly/2.)**2, 0.) and on_boundary

class inflow(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],-(Lx/2.)) and on_boundary

class outflow(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],Lx/2.) and on_boundary

wavenum=2*pi/(Ly/4.)
wavenum2=2*pi/(Ly/4.)
freq=2*pi/200.

wavenummod=wavenum
wavenum2mod=wavenum2
freqmod=freq

inflowExpr=Expression(("inflowVel + 0.1*sin(freq*t + wavenum*x[1])","0. + 0.1*sin(freq*t + wavenum2*x[1]) "), inflowVel=inflowVel,t=0,wavenum=wavenum,wavenum2=wavenum2,freq=freq,degree=2)
# inflowExpr=Expression(("inflowVel + 0.05*sin(2*pi*t/100. + wavenum*x[1]) + perturbx*0.2*sin(2*pi*t/100. + wavenum2*x[1]+pi/2.)","0. + 0.01*sin(2*pi*t/100. + wavenum*x[1])+ perturby*0.2*sin(2*pi*t/100. + wavenum2*x[1])"), inflowVel=inflowVel,t=0,perturbx=0,perturby=0,wavenum=wavenum,wavenum2=wavenum2,degree=2)

# inflowExpr=Expression(("inflowVel + 0.5*sin(2*pi*t/100. + wavenum*x[1])","0. + 0.25*sin(2*pi*t/100.)"), inflowVel=inflowVel,t=0,wavenum=wavenum,degree=2)
# inflowExpr=Expression(("inflowVel","0."), inflowVel=inflowVel,degree=2)

# lateral BC
bcu_inflow = DirichletBC(V, inflowExpr, inflow())
# bcu_walls = DirichletBC(V, Expression(("0","0."), inflowVel=inflowVel,degree=2), walls())
bcp_outflow = DirichletBC(Q, Constant(0), outflow())

# bc1a = DirichletBC(V.sub(1), Constant(0.0), NoSlipBoundary())

# inflow BC
# bc2 = DirichletBC(V, Constant((inflowVel,0.0)), InflowBoundary())
# bc2a = DirichletBC(VQ.sub(0).sub(0), Constant(8.), InflowBoundary())

# bcp = [DirichletBC(Q, Constant(0), OutflowBoundary())]
bcp=[bcp_outflow]
# bcu = [bcu_inflow,bcu_walls]
bcu = [bcu_inflow]


# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

mx,my,mz = createLayout(numturbs)

ma=[Constant(mm) for mm in 0.33*np.ones(numturbs)]

# right hand rule from above
# myaw=[Constant(pi/8.),Constant(0),Constant(0)]

yaw=0
myaw = [Constant(mm) for mm in (yaw*pi/180.)*np.ones(numturbs)]

beta = integrate.dblquad(WTGdist,-3*radius,3*radius,lambda x: -3*radius,lambda x: 3*radius)

B=beta[0]

f  = createRotatedTurbineForce(mx,my,ma,A,B,numturbs,alpha,V)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx \
   + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   + dot(f*(cos(myaw[0])**2*u_n[0]*u_n[0]+sin(myaw[0])**2*u_n[1]*u_n[1]), v)*dx  # inner?  other form of vel?
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Create XDMF files for visualization output
# ufile = File('output/fields/velocity_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.pvd')
# pfile = File('output/fields/pressure_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.pvd')
xdmffile_u = XDMFFile('output/velocity_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.xdmf')
xdmffile_p = XDMFFile('output/pressure_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.xdmf')
# # xdmffile_tf = XDMFFile('2DDynamic/turbine_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.xdmf')


# # Create time series (for use in reaction_system.py)
# timeseries_u = TimeSeries('output/velocity_series_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.xdmf')
# timeseries_p = TimeSeries('output/pressure_series_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha)+'.xdmf')

# # Save mesh to file (for use in reaction_system.py)
# File('navier_stokes_cylinder/cylinder.xml.gz') << mesh

# Create progress bar
# progress = Progress('Time-stepping')
# set_log_level(PROGRESS)

# ufile = File('output/u_'+str(float(mu))+'.pvd')
# pfile = File('output/p_'+str(float(mu))+'.pvd')

# DoF=len(u_.vector()[:])
# snapshots = np.zeros((DoF,int(num_steps/save_int)))

# uInterp = Function(V)
# uInterp=project(Expression(("x[0]","x[1]"),degree=2),V)
# basePositions=uInterp.vector()[:]
# np.save('output/basePositions_'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha),basePositions)

# Time-stepping
t = 0
count=0

for n in range(num_steps):

    # Update current time
    t += dt

    # bcu_inflow.perturbx=.1*np.random.rand()
    # bcu_inflow.perturby=.1*np.random.rand()
    inflowExpr.t=t
    # wavenummod = wavenummod + .01*np.random.randn()*wavenum
    # wavenum2mod = wavenum2mod+ .01*np.random.randn()*wavenum2
    # freqmod = freqmod+ .01*np.random.randn()*wavenum2
    # inflowExpr.wavenum=wavenummod
    # inflowExpr.wavenum2=wavenum2mod
    # inflowExpr.freq=freqmod


    bcu_inflow = DirichletBC(V, inflowExpr, inflow())
    bcu=[bcu_inflow]
    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)
    if n % save_int ==0:
        # Save solution to file (XDMF/HDF5)
        # ufile << u_
        # pfile << p_
        xdmffile_u.write(u_, t)
        xdmffile_p.write(p_, t)
        # xdmffile_tf.write(project(f,V),t)

        # # Save nodal values to file
        # timeseries_u.store(u_.vector(), t)
        # timeseries_p.store(p_.vector(), t)

        # snapshots[:,count]=u_.vector()[:]

        print(t)
        # print(wavenummod/wavenum)
        # print(wavenum2mod/wavenum2)
        # print(freqmod/freq)

        count+=1

    # # Update progress bar
    # progress.update(t / T)
    # print('u max:', u_.vector().array().max())

# Hold plot
# interactive()

# np.save('output/snapshots'+str(numturbs) + '_' + str(int(np.round(Re))) + '_' + str(yaw) + '_' + str(alpha),snapshots)