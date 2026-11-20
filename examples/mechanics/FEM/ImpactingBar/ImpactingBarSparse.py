import numpy as np
import siconos.integrators
import siconos.modeling as sm
import siconos.nonsmooth_formulations
import siconos.plot_config as sicoplot
import siconos.simulation
from scipy.sparse import csc_array

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)


# User-defined main parameters
ndof = 100  # degrees of freedom for the beam
t0 = 1e-8  # initial computation time
T = 0.0015  # final computation time
h = 2e-6  # time step
position_init = 0.0000503  # initial position (contact at 251.5 h: not on a time step)
velocity_init = -0.1  # initial velocity
epsilon = 0.0  # 1e-1
theta = 1.0 / 2.0 + epsilon  # theta for MoreauJeanOSI integrator
# theta = 1.0
E = 210e9  # young Modulus
S = 0.000314  # Beam Section 1 cm  for the diameter
# S=0.1
beam_length = 1.0  # length of the  beam
elem_length = beam_length / ndof  # length of an element
rho = 7800.0  # specific mass
# rho=1.0
g = 9.81  # Gravity
g = 0.0


# mass= SiconosMatrix(ndof,ndof,SPARSE,ndof)
# stiffness= SiconosMatrix(ndof,ndof,SPARSE,ndof)
nnz = 2 * (ndof - 1) + ndof
diag = 2.0 / 3.0 * rho * S * elem_length * np.ones(ndof, dtype=np.float64)
diag[0] *= 0.5
diag[-1] *= 0.5
upper = 1.0 / 6.0 * rho * S * elem_length * np.ones(ndof - 1, dtype=np.float64)
mass = csc_array((ndof, ndof), dtype=np.float64)
mass.setdiag(diag)
mass.setdiag(upper, 1)
mass.setdiag(upper, -1)

diag = 2.0 * E * S / elem_length * np.ones(ndof, dtype=np.float64)
diag[0] *= 0.5
diag[-1] *= 0.5
upper = -E * S / elem_length * np.ones(ndof - 1, dtype=np.float64)
stiffness = csc_array((ndof, ndof), dtype=np.float64)
stiffness.setdiag(diag)
stiffness.setdiag(upper, 1)
stiffness.setdiag(upper, -1)

q0 = np.full((ndof), position_init)
v0 = np.full((ndof), velocity_init)

bar = sm.LagrangianSparseLinearTIDS(q0, v0, sm.alias_t)
bar.setConstantMass(mass, sm.alias_t)
bar.setStiffnessMatrix(stiffness, sm.alias_t)
# bar.display()

weight = np.full((ndof), -g * rho * S / elem_length)
bar.setConstantFext(weight, sm.alias_t)

e = 0.0

H = np.zeros((1, ndof))
H[0, 0] = 1.0

nslaw = sm.NewtonImpactNSL(e)
relation = sm.LagrangianLinearTIR(H)
inter = sm.Interaction(nslaw, relation)

# -------------
# --- Model ---
# -------------
impactingBar = sm.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system in the non smooth dynamical system
impactingBar.insertDynamicalSystem(bar)

# link the interaction and the dynamical system
impactingBar.link(inter, bar)


# ------------------
# --- Simulation ---
# ------------------

# -- (1) OneStepIntegrators --
OSI = siconos.integrators.MoreauJeanOSI(theta, 0.0)

# -- (2) Time discretisation --
t = siconos.simulation.TimeDiscretisation(t0, h)

# -- (3) one step non smooth problem
osnspb = siconos.nonsmooth_formulations.LCP()

s = siconos.simulation.TimeStepping(impactingBar, t, OSI, osnspb)

k = 0

N = int((T - t0) / h)
dataPlot = np.zeros((N, 12))
halfpos = int(ndof * 0.5)

q = bar.q()
v = bar.velocity()
p = bar.p(1)
lambda_ = inter.lambda_python(1)
dataPlot[k, 0] = t0
dataPlot[k, 1] = q[0]
dataPlot[k, 2] = v[0]
dataPlot[k, 3] = p[0]
dataPlot[k, 7] = q[ndof - 1]
dataPlot[k, 8] = v[ndof - 1]
dataPlot[k, 9] = q[halfpos]
dataPlot[k, 10] = v[halfpos]
potentialEnergy = q.dot(stiffness @ q)
kineticEnergy = v.dot(mass @ v)
dataPlot[k, 5] = potentialEnergy
dataPlot[k, 6] = kineticEnergy
k = 1
# time loop
while k < N:
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()
    # print("time=", dataPlot[k, 0])
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0] / h
    dataPlot[k, 4] = lambda_[0]
    dataPlot[k, 7] = q[ndof - 1]
    dataPlot[k, 8] = v[ndof - 1]
    dataPlot[k, 9] = q[halfpos]
    dataPlot[k, 10] = v[halfpos]
    potentialEnergy = q.dot(stiffness @ q)
    kineticEnergy = v.dot(mass @ v)
    dataPlot[k, 5] = potentialEnergy
    dataPlot[k, 6] = kineticEnergy
    v_p_theta = theta * v[0] + (1 - theta) * dataPlot[k - 1, 2]
    impactEnergy = v_p_theta * lambda_[0]
    dataPlot[k, 11] = impactEnergy

    k += 1
    s.nextStep()

ref = np.loadtxt("ImpactingBar.ref", skiprows=1)
assert np.allclose(ref, dataPlot)


fig_size = [14, 14]
plt.rcParams["figure.figsize"] = fig_size

plt.subplot(411)
plt.title("position")
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title("velocity")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title("reaction")
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title("lambda")
plt.grid()

if enable_plot:
    plt.show()
