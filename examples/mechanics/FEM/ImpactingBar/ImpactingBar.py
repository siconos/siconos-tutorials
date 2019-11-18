from numpy.linalg import norm
from siconos.kernel import LagrangianLinearTIDS, NewtonImpactNSL,\
    LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem, MoreauJeanOSI,\
    TimeDiscretisation, LCP, TimeStepping

from siconos.kernel import SimpleMatrix, getMatrix, SPARSE
import numpy as np

# User-defined main parameters
nDof = 500  #degrees of freedom for the beam
t0 = 1e-8    #initial computation time
T = 0.0015                  # final computation time
h = 2e-7                # time step
position_init = 0.00005      # initial position
velocity_init =  -.1      # initial velocity
epsilon = 0.5#1e-1
theta = 1/2.0 + epsilon              # theta for MoreauJeanOSI integrator
#theta = 1.0
E = 210e9 # young Modulus
S = 0.000314 #  Beam Section 1 cm  for the diameter
#S=0.1
L = 1.0 # length of the  beam
l = L/nDof # length of an element
rho = 7800.0  # specific mass
#rho=1.0
g = 9.81 # Gravity
g=0.0


M= SimpleMatrix(nDof,nDof,SPARSE,nDof)
K= SimpleMatrix(nDof,nDof,SPARSE,nDof)
K.setValue(0,0, 1.*E*S/l)
K.setValue(0,1,-1.*E*S/l)
M.setValue(0,0, 1/3.*rho*S*l)
M.setValue(0,1, 1/6.*rho*S*l)

for i in range(1,nDof-1):
    K.setValue(i,i,2.*E*S/l)
    K.setValue(i,i-1,-1.*E*S/l)
    K.setValue(i,i+1,-1.*E*S/l)
    M.setValue(i,i,2/3.*rho*S*l)
    M.setValue(i,i-1,1/6.*rho*S*l)
    M.setValue(i,i+1,1/6.*rho*S*l)


K.setValue(nDof-1,nDof-2,-1.*E*S/l)
K.setValue(nDof-1,nDof-1, 1.*E*S/l)
M.setValue(nDof-1,nDof-2,1/6.*rho*S*l)
M.setValue(nDof-1,nDof-1,1/3.*rho*S*l)


q0 = np.full((nDof), position_init)
v0 = np.full((nDof), velocity_init)

bar = LagrangianLinearTIDS(q0,v0,M)
bar.setKPtr(K)
#bar.display()

weight = np.full((nDof),-g*rho*S/l)
bar.setFExtPtr(weight)

e=0.0

H = np.zeros((1,nDof))
H[0,0]=1.

nslaw = NewtonImpactNSL(e)
relation = LagrangianLinearTIR(H)
inter = Interaction(nslaw, relation)

# -------------
# --- Model ---
# -------------
impactingBar = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system in the non smooth dynamical system
impactingBar.insertDynamicalSystem(bar);

# link the interaction and the dynamical system
impactingBar.link(inter,bar);


# ------------------
# --- Simulation ---
# ------------------

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

q = bar.q()
v = bar.velocity()
p = bar.p(1)
lambda_ = inter.lambda_(1)

# time loop
while s.hasNextEvent():
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()
    print('time=', dataPlot[k, 0])
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0]/h
    dataPlot[k, 4] = lambda_[0]

    k += 1
    s.nextStep()

dataPlot.resize(k,5)


import matplotlib.pyplot as plt

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
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title('lambda')
plt.grid()


plt.show()
