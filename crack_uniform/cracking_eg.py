# code solves for crack morphologies in drying complex droplets
# lubrication theory applied to reduce equations to 2D, compression acts as an effective plane strain
# vanilla version of phase-field model for fracture applied.
# gov eqns are: \nabla \cdot ((1-d)^2 \sigma) = 0, and 
#           \eta \partial_t d + \Gamma_c/l (d-l^2 \nabla^2 d) = (1-d)D, D= (\sigma-p):(e-p/2)
# finite element code solved using fenics

from dolfin import *
from mshr import *
import numpy as np

#---------------------------------------------------------
#set fenics relevant parameters
set_log_level(ERROR)
parameters.parse()   # read paramaters from command line
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"]=1
parameters["form_compiler"]["cpp_optimize"] = True
solver_parameters =  {"linear_solver" : "cg",
                            "symmetric" : True,
                            "preconditioner" : "jacobi",
                            "krylov_solver" : {
                                "report" : False,
                                "monitor_convergence" : False,
                                "relative_tolerance" : 1e-6
                                }
                            }

#---------------------------------------------------------
# mesh generation
rad = 1.0; res=150;tol=1e-5;
geom = Circle(Point(0., 0.), rad)
mesh = generate_mesh(geom,res)
ndim = 2 # get number of space dimensions

#---------------------------------------------------------
# physical constants for fracture
Gc = Constant(10e-2);
kappa=16.;
ellv = 0.05; ell = Constant(ellv)
#residual stiffness
k_ell = 1e-6

# evap steps for time normalisation
phi_s0=0.45;
load_min = phi_s0 # load multiplier min value
load_max = 1. # load multiplier max value
load_steps = 5501 # number of time steps
p_m= 100 # points at which to print 
dt = (load_max-load_min)/(load_steps-1)
eta = 1.0e-5
# Numerical parameters of the alternate minimization
maxiter = 200

# printing parameters
clc=0.
p_t=0
time_print = 0
#---------------------------------------------------------
# set the save directory
savedir = "circle-Fixslip-rad%.1f-Gc%.2f-phi%.2f-k%.2f"%(rad,Gc,phi_s0,kappa)
# loading and initialization of vectors to store time datas
load_multipliers = np.linspace(load_min,load_max,load_steps)
energies = np.zeros((len(load_multipliers),5))
file_d = File(savedir+"/d.pvd","compressed") # use .pvd if .xdmf is not working
file_u = File(savedir+"/u.pvd","compressed") # use .pvd if .xdmf is not working

#---------------------------------------------------------
# set for boundary conditions
class Edg(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and sqrt(x[0]**2+x[1]**2) > sqrt((rad - tol )**2)

#---------------------------------------------------------
# constit fns for linear elasticity
def eps(v):
    return 0.5*(nabla_grad(v)+nabla_grad(v).T)
def sig(eps):
    return 2*eps
# constit functions of the damage model
def w(d):
    return d**2
def a(d):
    return (1-d)**2
def diss_energy(d):
    return 0.5*( w(d)+(ell**2.)*dot(grad(d),grad(d)) )
def elas_energy(u,d):
    return 0.5*inner((a(d)*(1.-k_ell)+k_ell)*sig(eps(u)-p/2.),eps(u)-p/2.)+0.5*kappa*inner(u,u)

# Create function space for 2D elasticity + Damage
V_u = VectorFunctionSpace(mesh, "Lagrange", 1)
V_d = FunctionSpace(mesh, "Lagrange", 1)
d, d_0  = Function(V_d), Function(V_d)
d_0 = interpolate(Expression('0.',degree=1),V_d)
u = Function(V_u)
#---------------------------------------------------------
# problem defined below here
phi_s = Expression("t", t = 0.0, degree = 1)
f_phi = Expression("(X-Y)/X",X = phi_s0, Y = phi_s,degree=1)
p = 2.*f_phi*Identity(ndim)


#solves for u given d and p
def u_solution(d,p,V_u):
	u_temp, d_u, v = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
	zero_v = Constant((0.,)*ndim)
	u_0 = interpolate(zero_v,V_u)
	bc_u = DirichletBC(V_u,u_0,Edg())
	#define elastic energy for variational deriv
	elastic_energy = elas_energy(u_temp,d)*dx
	E_u = derivative(elastic_energy,u_temp,v)
	E_du = replace(E_u,{u_temp:d_u})
	problem_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), u_temp, bc_u)
	solver_u = LinearVariationalSolver(problem_u)
	solver_u.parameters.update(solver_parameters)
	return solver_u.solve()

def d_solution(u,d_0,p,V_d):
	d, d_d, beta = Function(V_d), TrialFunction(V_d), TestFunction(V_d)
	# set d
	a0 = Expression('0.0',degree = 1)
	bc_d = DirichletBC(V_d,a0,Edg())
	driv = inner(a(d)*sig(eps(u)-p/2.),eps(u)-p/2.)*dx
	dissipated_energy = (Gc/ell)*diss_energy(d)*dx
	d_energy = (eta/(2.*dt))*inner(d-d_0,d-d_0)*dx + dissipated_energy + 0.5*driv
	E_d = derivative(d_energy,d,beta)
	E_d_d = derivative(E_d,d,d_d)
	E_dd = replace(E_d,{d:d_d})
	problem_d = LinearVariationalProblem(lhs(E_dd), rhs(E_dd), d,bc_d)
	# Set up the solvers
	solver_d = LinearVariationalSolver(problem_d)
	solver_d.parameters.update(solver_parameters_d)
	return solver_d.solve()


# defns for post processing
dissipated_energy = (Gc/ell)*diss_energy(d)*dx
elastic_energy = elas_energy(u,d)*dx
total_energy = elastic_energy + dissipated_energy
dmg_area = d*dx

# error function defined here
d_error = Function(V_d)
d_i = d_0

for (i_t, t) in enumerate(load_multipliers):
     phi_s.t=t
     # Initialise 
     iter = 1; err_d = 1
     while err_d>tol and iter<maxiter:
          # solve elastic problem
          u = u_solution(d,p,V_u)
          # solve damage problem
          d=d_solution(u,d_0,p,V_d)
          # test error
          d_error.vector()[:] = d.vector() - d_i.vector()
          err_d = np.linalg.norm(d_error.vector().get_local(), ord = np.Inf)
          print("Iteration:  %2d, Error: %2.8g, d_max: %.8g, u_max: %.8g" %(iter, err_d, d.vector().max(), u.vector().max()))
          iter=iter+1
          d_i.assign(d)
     d_0.assign(d)
     u_0.assign(u)
     print("\nEnd of timestep %d with load multiplier %g"%(i_t, t))
     print("-----------------------------------------")
     # post-processing
     elastic_energy_value = assemble(elastic_energy)
     surface_energy_value = assemble(dissipated_energy)
     d_area = assemble(dmg_area)
     energies[i_t] = np.array([t,elastic_energy_value,surface_energy_value,elastic_energy_value+surface_energy_value,d_area])
     # Calculate the axial force resultant
     np.savetxt(savedir+'/energies.txt', energies)
     if (max_d>=0.99):
          if (p_t % p_m ==0):
                file_d << (d,t)
                file_u << (u,t)
          p_t+=1
     if (max_d>= 1.-1e-4 and time_print == 0):
          np.savetxt(savedir+'/nuc.txt', [t])
          time_print = 10
     if (t == 1.):
          file_d << (d,t)
          file_u << (u,t)
