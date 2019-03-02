# code solves for crack morphologies in drying complex droplets
# lubrication theory applied to reduce equations to 2D, compression acts as an effective plane strain
# vanilla version of phase-field model for fracture applied.
# gov eqns are: \nabla \cdot ((1-d)^2 \sigma) = 0, and 
#           \eta \partial_t d + \Gamma_c/l (d-l^2 \nabla^2 d) = (1-d)D, D= (\sigma-p):(e-p/2)
# finite element code solved using fenics
# read constants with yaml

from dolfin import *
from mshr import *
import sys
import numpy as np
import math
import os
import sympy

#---------mat parameters
Kap_raw=.1
Ga_raw=Constant(1e-1)
T_raw=1.e-1
K_raw=1.e-2
G_raw=0.1
alp_raw=1.5
phi_raw=0.2

TIM = 1.6
N_t = 10001
#----------------------


# ------------------
# Parameters
# ------------------
set_log_level(ERROR) # log level
parameters.parse()   # read paramaters from command line
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"]=1
parameters["form_compiler"]["cpp_optimize"] = True
# parameters for linear solver - mess around with this # bicgstab is also good to use
#gmres is default used petsc_amg as precon, try ilu for now as well
solver_parameters =  {"linear_solver" : "cg",
                            "symmetric" : True,
                            "preconditioner" : "jacobi",
                            "krylov_solver" : {
                                "report" : False,
                                "monitor_convergence" : False,
                                "relative_tolerance" : 1e-6
                                }
                            }

# next we define the mesh and key parameters for the problem
#----------geometry-----------
rad = 1.0; cell_size=0.01;res=150;tol=1e-4;Or=1e-4
#---------mesh generation----
geom = Circle(Point(0., 0.), rad)
mesh = generate_mesh(geom,res)

#attempt to define normal vector (doesn't work)
n = FacetNormal(mesh)

ndim = 2 # get number of space dimensions

#------damage width-------------
ellv = 0.05; ell = Constant(ellv)

#-----material constants-----
Gc = 0.1;
kappa=4.;

#-------evap parameters--------
Theta=0.5;
K=0.01;
G = 0.1;
alp=1.5
phi_s0=0.4;

# we should be changing this one
sig_c = 1.
k_ell = Constant(1.e-4) # residual stiffness

time_min = 0.
time_max = 1.2
time_steps = 8001
dt = (time_max-time_min)/(time_steps-1.)
time_multipliers = np.linspace(time_min,time_max,time_steps)
p_m=200

eta = 1.0e-5

# Numerical parameters of the alternate minimization
maxiter = 2
toll = 1e-5

# set the save directory
savedir = "circle-theta%.2f-k%.2f-G%.2f-Gc%.2f-kap%.2f-phi%.2f"%(Theta,K,G,Gc,kappa,phi_s0)
#  loading and initialization of vectors to store time datas
energies = np.zeros((len(time_multipliers),6))
metrics = np.zeros((len(time_multipliers),5))
file_d = File(savedir+"/d.pvd","compressed") # use .pvd if .xdmf is not working
file_u = File(savedir+"/u.pvd","compressed") # use .pvd if .xdmf is not working
file_h = File(savedir+"/h.pvd","compressed")

#----------------------------------------------------------------------------
# Define boundary sets for boundary conditions
#----------------------------------------------------------------------------
class Edg(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and sqrt(x[0]**2+x[1]**2) > sqrt((rad - tol )**2)

# parameters ``mu`` and ``lmbda`` to avoid re-generation of C++ code
# when changing model parameters. Note that ``lambda`` is a reserved
# keyword in Python, hence the misspelling ``lmbda``.

#constit fns for elasticity
def eps(v):
    return 0.5*(nabla_grad(v)+nabla_grad(v).T)

def sig(eps):
    return 2*eps #+ lmbda*tr(eps)*Identity(ndim)

# Constitutive functions of the damage model
def w(d):
    return d**2.

def a(d):
    return (1-d)**2.

def diss_energy(d):
    return 0.5*( w(d)+(ell**2.)*dot(grad(d),grad(d)) )

def Mx(a,b):
    return 0.5*((a+b)+abs(a-b))

#-------------Useful definitions-------------------
# zero and unit vectors
zero_v = Constant((0.,)*ndim)
#e1 = [Constant([1.,0.]),Constant((1.,0.,0.))][ndim-2]
ex = Constant([1.,0.])
ey = Constant([0.,1.])

#-------create vector spaces------------
# Create function space for 2D elasticity + Damage
V_u = VectorFunctionSpace(mesh, "Lagrange", 1)
V_d = FunctionSpace(mesh, "Lagrange", 1)
V_h = FunctionSpace(mesh, "Lagrange",1)

# Define the function, test and trial fields
u, du, v = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
d, dd, beta = Function(V_d), TrialFunction(V_d), TestFunction(V_d)
h, g, dh = Function(V_h), TestFunction(V_h),TrialFunction(V_h)
phi = Function(V_h)
phi_s = Function(V_h)

# h problem
rr = Expression('x[0]*x[0]+x[1]*x[1]',degree=1)
ic = Expression('1.-y',y=rr,degree=1)

h_0 = project(ic,V_h)
h_n = project(ic,V_h)
h   = project(ic,V_h)
Bc0 = Expression('0.0',degree = 1)
#-------------boundary condns-----------------
bc = DirichletBC(V_h,Bc0,Edg())
phi = 1. - phi_s0*(h_0*h**-1.)
phi_s = phi_s0*(h_0*h**-1.)
#------gov eqn-------------------------------
F = ( inner((h-h_n)/dt,g) +  (1.+G*a(d)/phi_s0)*inner(h*grad(phi)*(phi_s)**(-alp),grad(g)) + Theta*phi*inner(1./(K+h),g) )*dx

# u problem
u_0 = interpolate(zero_v,V_u)
bc_u = DirichletBC(V_u,u_0,Edg())
f_phi = (phi_s0 - phi_s0*(h_0*h**-1.))/phi_s0
p = 2.*f_phi*Identity(ndim)
elastic_energy = 0.5*inner(h*a(d)*sig(eps(u)-p/2.),eps(u)-p/2.)*dx+0.5*kappa*inner(u,u)*dx
E_u = derivative(elastic_energy,u,v)
E_du = replace(E_u,{u:du})
problem_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), u)#, bc_u) # bc is dependent upon the study
solver_u      = LinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_parameters)

# d problem
a0 = Expression('0.0',degree = 1)
d_0 = interpolate(a0,V_d)
bc_d = DirichletBC(V_d,a0,Edg())
dissipated_energy = (Gc/ell)*diss_energy(d)*dx
#define dmg problems
driv_temp = inner(sig(eps(u)-p/2.),eps(u)-p/2.)
driv = Mx(driv_temp/sig_c-1.,0.)
dissipated_energy = (Gc/ell)*diss_energy(d)*dx
d_energy = (eta/(2.*dt))*inner(d-d_0,d-d_0)*dx + dissipated_energy + 0.5*a(d)*driv*dx
E_d = derivative(d_energy,d,beta)
E_d_d = derivative(E_d,d,dd)
E_dd = replace(E_d,{d:dd})
problem_d=LinearVariationalProblem(lhs(E_dd), rhs(E_dd), d,bc_d)
solver_d   = LinearVariationalSolver(problem_d)
solver_d.parameters.update(solver_parameters)


# post-processing
total_energy = elastic_energy + dissipated_energy
dmg_area = d*dx
drop_mass = h*dx
h_error = Function(V_h)
diff_h = 0.
d_i = d_0
d_error = Function(V_d)

clc=0.
time_print=0.
fin_print=0
p_t=0.
p_c=0.

initial_mass = assemble(drop_mass)
init_fluid_mass = (1-phi_s0)*initial_mass
solid_mass = phi_s0*initial_mass

initial_stats = np.array([initial_mass,init_fluid_mass,solid_mass])
np.savetxt(savedir+'/init.txt',initial_stats)

for (i_t, t) in enumerate(time_multipliers):
     # Initialization
     iter = 1; err_d = 1
     #solve h problem
     solve(F == 0., h, bc)
     h_n.assign(h)
     while err_d>tol and iter<maxiter:
          # solve elastic problem
          solver_u.solve()
          # solve damage problem
          solver_d.solve()
          # test error
          d_error.vector()[:] = d.vector() - d_i.vector()
          err_d = np.linalg.norm(d_error.vector().get_local(), ord = np.Inf)
          print("Iteration:  %2d, Error: %2.8g" %(iter, err_d))
          iter=iter+1
          d_i.assign(d)
     d_0.assign(d)
     u_0.assign(u)
     h_error.vector()[:] = h.vector() - h_n.vector()
     diff_h = np.linalg.norm(h_error.vector().get_local(), ord = np.Inf)
     h_max = h.vector().max()
     print("Time: %.6f, h_max: %.8g, d_max: %.8g, u_max: %.8g" %(t, h_max, d.vector().max(),u.vector().max()))
     print("\nEnd of timestep %d with load multiplier %g" %(i_t, t))
     print("-----------------------------------------")
     # Some post-processing
     # Calculate the energies
     elastic_energy_value = assemble(elastic_energy)
     surface_energy_value = assemble(dissipated_energy)
     d_area = assemble(dmg_area)
     tot_mass = assemble(drop_mass)
     fluid_mass = tot_mass - solid_mass
     energies[i_t] = np.array([t,elastic_energy_value,surface_energy_value,elastic_energy_value+surface_energy_value,tot_mass,d_area])
     metrics[i_t] = np.array([t,tot_mass,fluid_mass,h_max,d.vector().max()])
     # Calculate the axial force resultant
     max_d = d.vector().max()
     np.savetxt(savedir+'/energies.txt', energies)
     np.savetxt(savedir+'/metrics.txt', metrics)
     max_d = d.vector().max()
     if (fluid_mass <= ((7-p_t)/8)*init_fluid_mass+toll and p_t == p_c):
          p_t+=1
          p_c+=1
          file_d << (d,t)
          file_u << (u,t)
          file_h << (h,t)
          phi_p=project(phi_s0*(h_0/h),FunctionSpace(mesh,'CG',1))
          file_phi << (phi_p,t)
          if (p_c == 7):
              p_c=20.
     if (max_d>= 1.-1e-4 and time_print == 0):
          np.savetxt(savedir+'/nuc.txt', [t, elastic_energy_value+surface_energy_value,d_area,tot_mass,h.vector().max()])
          file_d << (d,t)
          file_u << (u,t)
          file_h << (h,t)
          phi_p=project(phi_s0*(h_0/h),FunctionSpace(mesh,'CG',1))
          file_phi << (phi_p,t)
          time_print = 10
     if (diff_h <= 5e-5 and h_max <= phi_s0+1e-4 and fin_print==0):
          break

np.savetxt(savedir+'/fin.txt', [Theta,i_t,t,phi_s0-h_max])
file_h << (h,t)
file_d << (d,t)
file_u << (u,t)
phi_p=project(phi_s0*(h_0/h),FunctionSpace(mesh,'CG',1))
file_phi << (phi_p,t)  
