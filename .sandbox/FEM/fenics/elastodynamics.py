from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Form compiler options
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Define mesh
mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 3, 2, 2)

# Sub domain for clamp at left end
def left(x, on_boundary):
    return near(x[0], 0.) and on_boundary

# Sub domain for rotation at right end
def right(x, on_boundary):
    return near(x[0], 1.) and on_boundary


# Elastic parameters
E  = 1000.0
nu = 0.3
mu    = Constant(E / (2.0*(1.0 + nu)))
lmbda = Constant(E*nu / ((1.0 + nu)*(1.0 - 2.0*nu)))

# Mass density
rho = Constant(1.0)

# Rayleigh damping coefficients
eta_m = Constant(0.)
eta_k = Constant(0.)

# Generalized-alpha method parameters
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)

# Time-stepping parameters
T       = 4.0
Nsteps  = 50

T       = 4.0
Nsteps  = 100

dt = Constant(T/Nsteps)


p0 = 1.
cutoff_Tc = T/5
# Define the loading as an expression depending on t
p = Expression(("0", "t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)



# Define function space for displacement, velocity and acceleration
V = VectorFunctionSpace(mesh, "CG", 1)
# Define function space for stresses
Vsig = TensorFunctionSpace(mesh, "DG", 0)


# Test and trial functions
du = TrialFunction(V)
u_ = TestFunction(V)
# Current (unknown) displacement
u = Function(V, name="Displacement")
# Fields from previous time step (displacement, velocity, acceleration)
u_old = Function(V)
v_old = Function(V)
a_old = Function(V)

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
force_boundary = AutoSubDomain(right)
force_boundary.mark(boundary_subdomains, 3)

# Define measure for boundary condition integral
dss = ds(subdomain_data=boundary_subdomains)

# Set up boundary condition at left end
zero = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, zero, left)

# Stress tensor
def sigma(r):
    return 2.0*mu*sym(grad(r)) + lmbda*tr(sym(grad(r)))*Identity(len(r))

# Mass form
def m(u, u_):
    return rho*inner(u, u_)*dx

# Elastic stiffness form
def k(u, u_):
    return inner(sigma(u), sym(grad(u_)))*dx

# Rayleigh damping form
def c(u, u_):
    return eta_m*m(u, u_) + eta_k*k(u, u_)

# Work of external forces
def Wext(u_):
    return dot(u_, p)*dss(3)

# Update formula for acceleration
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return (u-u_old-dt_*v_old)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a_old

# Update formula for velocity
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u, u_old, v_old, a_old):
    """Update fields at the end of each time step."""

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # use update functions using vector arguments
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # Update (u_old <- u)
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new


def fenics_time_integration():


    # Residual
    a_new = update_a(du, u_old, v_old, a_old, ufl=True)
    v_new = update_v(a_new, u_old, v_old, a_old, ufl=True)
    res = m(avg(a_old, a_new, alpha_m), u_) + c(avg(v_old, v_new, alpha_f), u_) \
           + k(avg(u_old, du, alpha_f), u_) - Wext(u_)
    a_form = lhs(res)
    L_form = rhs(res)


    # Define solver for reusing factorization
    K, res = assemble_system(a_form, L_form, bc)

    solver = LUSolver(K)
    # solver = LUSolver(K, "mumps")
    solver.parameters["symmetric"] = True

    # Time-stepping
    time = np.linspace(0, T, Nsteps+1)
    u_tip = np.zeros((Nsteps+1,))
    energies = np.zeros((Nsteps+1, 4))
    E_damp = 0
    E_ext = 0
    sig = Function(Vsig, name="sigma")
    xdmf_file = XDMFFile("elastodynamics-results.xdmf")
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True
    xdmf_file.parameters["rewrite_function_mesh"] = False



    def local_project(v, V, u=None):
        """Element-wise projection using LocalSolver"""
        dv = TrialFunction(V)
        v_ = TestFunction(V)
        a_proj = inner(dv, v_)*dx
        b_proj = inner(v, v_)*dx
        solver = LocalSolver(a_proj, b_proj)
        solver.factorize()
        if u is None:
            u = Function(V)
            solver.solve_local_rhs(u)
            return u
        else:
            solver.solve_local_rhs(u)
            return

    for (i, dt) in enumerate(np.diff(time)):

        t = time[i+1]
        print("Time: ", t)

        # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
        p.t = t-float(alpha_f*dt)

        # Solve for new displacement
        res = assemble(L_form)
        bc.apply(res)
        solver.solve(K, u.vector(), res)


        # Update old fields with new quantities
        update_fields(u, u_old, v_old, a_old)

        # Save solution to XDMF format
        xdmf_file.write(u, t)


        File('elastodynamics/displacement_{}.pvd'.format(i)) << u


        # Compute stresses and save to file
        local_project(sigma(u), Vsig, sig)
        xdmf_file.write(sig, t)

        File('elastodynamics/stress_{}.pvd'.format(i)) << sig

        p.t = t
        # Record tip displacement and compute energies
        u_tip[i+1] = u(1., 0.05, 0.)[1]
        E_elas = assemble(0.5*k(u_old, u_old))
        E_kin = assemble(0.5*m(v_old, v_old))
        E_damp += dt*assemble(c(v_old, v_old))
        # E_ext += assemble(Wext(u-u_old))
        E_tot = E_elas+E_kin+E_damp #-E_ext
        energies[i+1, :] = np.array([E_elas, E_kin, E_damp, E_tot])




    # Plot tip displacement evolution
    plt.figure()
    plt.plot(time, u_tip)
    plt.xlabel("Time")
    plt.ylabel("Tip displacement")
    plt.ylim(-0.5, 0.5)
    plt.show()

    # Plot energies evolution
    plt.figure()
    plt.plot(time, energies)
    plt.legend(("elastic", "kinetic", "damping", "total"))
    plt.xlabel("Time")
    plt.ylabel("Energies")
    plt.ylim(0, 0.0011)
    plt.show()




#fenics_time_integration()



u1 = TrialFunction(V)
u2 = TestFunction(V)

print(' -- Assemble mass and stiffness matrix in Fenics -- ')


stiffness_mat = assemble(k(u1,u2))
bc.apply(stiffness_mat)
print('stiffness matrix : ', stiffness_mat.array())

mass_mat = assemble(m(u1,u2))
bc.apply(mass_mat)
print('mass matrix : ',mass_mat.array())


stiffness_mat_np = np.array(stiffness_mat.array())
mass_mat_np = np.array(mass_mat.array())

from siconos.kernel import LagrangianLinearTIDS, NewtonImpactNSL,\
    LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem, MoreauJeanOSI,\
    TimeDiscretisation, LCP, TimeStepping

from siconos.kernel import SimpleMatrix, getMatrix, SPARSE #, SPARSE_COORDINATE


print(' -- create mass and stiffness matrix in Siconos -- ')

n_dof = mass_mat_np.shape[0]

print('n_dof=', n_dof)

# M= SimpleMatrix(n_dof,n_dof,SPARSE_COORDINATE,n_dof)
# K= SimpleMatrix(n_dof,n_dof,SPARSE_COORDINATE,n_dof)
#M= SimpleMatrix(n_dof,n_dof,SPARSE,n_dof)
#K= SimpleMatrix(n_dof,n_dof,SPARSE,n_dof)

M= SimpleMatrix(n_dof,n_dof)
K= SimpleMatrix(n_dof,n_dof)

#input()
#nnz=0
# for i in range(n_dof):
    
#     idx = np.where(stiffness_mat_np[i,:]>=1e-14)
#     for j in idx[0]:
#         K.setValue(i,j,stiffness_mat_np[i,j])

        
#     idx = np.where(mass_mat_np[i,:]>=1e-14)
#     nnz = nnz + len(idx[0])
#     for j in idx[0]:
#         M.setValue(i,j,mass_mat_np[i,j])
# nnz=0

for i in range(n_dof):
    for j in range(n_dof):
        K.setValue(i,j,stiffness_mat_np[i,j])
        M.setValue(i,j,mass_mat_np[i,j])

#print('nnz=', nnz)
#M.display()

print(' -- set initial conditions -- ')

q_0 = u_old.vector().get_local()
v_0 = v_old.vector().get_local()

t0 = 0.0

print(' -- build dynamical system -- ')

body = LagrangianLinearTIDS(q_0,v_0,M)
body.setKPtr(K)

applied_force = np.zeros(n_dof)
applied_force[n_dof-1]=1.
body.setFExtPtr(applied_force)

# -------------
# --- Model ---
# -------------
impactingBar = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system in the non smooth dynamical system
impactingBar.insertDynamicalSystem(body);

# link the interaction and the dynamical system
#impactingBar.link(inter,bar);


# ------------------
# --- Simulation ---
# ------------------

theta=0.5
h = 0.01
# -- (1) OneStepIntegrators --
OSI = MoreauJeanOSI(theta,0.5)

# -- (2) Time discretisation --
t = TimeDiscretisation(t0,h)

# -- (3) one step non smooth problem
osnspb = LCP()

s = TimeStepping(impactingBar, t,OSI,osnspb)

k =0

N = int((T-t0)/h)
dataPlot = np.zeros((N+1, 5))

q = body.q()
v = body.velocity()
p = body.p(1)
#lambda_ = inter.lambda_(1)

# time loop
while s.hasNextEvent():
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()
    print('time=', dataPlot[k, 0])
    print('q=',q)
    dataPlot[k, 1] = q[n_dof-1]
    dataPlot[k, 2] = v[n_dof-1]
    dataPlot[k, 3] = p[n_dof-1]/h
    #dataPlot[k, 4] = lambda_[0]

    k += 1
    s.nextStep()

dataPlot.resize(k,5)


fig_size = [14, 14]
plt.rcParams["figure.figsize"] = fig_size

plt.subplot(411)
plt.title('position')
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title('velocity')
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title('reaction')
plt.grid()

plt.show()
