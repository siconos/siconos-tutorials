# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

import numpy as np

from siconos.modeling import (
    LagrangianLinearTIDS,
    LagrangianDS,
    BoundaryCondition,
    NewtonImpactNSL,
    LagrangianLinearTIR,
    Interaction,
    NonSmoothDynamicalSystem,
    alias_t)

from siconos.integrators import (
    MoreauJeanOSI)
from siconos.simulation import (
    TimeDiscretisation,
    TimeStepping)


from siconos.nonsmooth_formulations import    LCP


# ==========================================================
# Parameters
# ==========================================================

nDof = 3

t0 = 0.0
T = 10.0
h = 0.005

theta = 0.5

R = 0.1
m = 1.0
g = 9.81

position_init = 1.0
velocity_init = 0.0

# ==========================================================
# Mass matrix
# ==========================================================

M = np.zeros((nDof, nDof), dtype=np.float64, order="F")
M[0, 0] = m
M[1, 1] = m
M[2, 2] = 2.0 / 5.0 * m * R * R

# ==========================================================
# Ball
# ==========================================================

q0 = np.zeros(nDof)
q0[0] = position_init

v0 = np.zeros(nDof)
v0[0] = velocity_init

ball = LagrangianLinearTIDS(q0, v0, M, alias_t)

weight = np.zeros(nDof, dtype=np.float64)
weight[0] = -m * g

ball.setConstantFext(weight, alias_t )

# ==========================================================
# Moving plane
# ==========================================================

q02 = np.zeros(nDof)

v02 = np.zeros(nDof)
v02[0] = -velocity_init

movingplane = LagrangianDS(q02, v02, alias_t)

movingplane.setConstantMass(M, alias_t)
movingplane.setConstantFext(weight, alias_t)

# ==========================================================
# Boundary condition
# ==========================================================

def prescribed_velocity(time, result):
    result[:]=  np.array(
        [2.0 + np.cos(0.5 * np.pi * time)]
    )
    

bc = BoundaryCondition([0])

bc.setComputePrescribedVelocityFunction(prescribed_velocity)

bc.computePrescribedVelocity(2.0)

print(bc.prescribedVelocity)
movingplane.setBoundaryConditions(bc)




# ==========================================================
# Interaction
# ==========================================================

e = 0.9

H = np.zeros((1, 2 * nDof))
H[0, 0] = 1.0
H[0, 3] = -1.0

nslaw = NewtonImpactNSL(e)
relation = LagrangianLinearTIR(H)

inter = Interaction(nslaw, relation)

# ==========================================================
# Model
# ==========================================================

nsds = NonSmoothDynamicalSystem(t0, T)

nsds.insertDynamicalSystem(ball)
nsds.insertDynamicalSystem(movingplane)

nsds.link(inter, ball, movingplane)

# ==========================================================
# Simulation
# ==========================================================

osi = MoreauJeanOSI(theta)

td = TimeDiscretisation(t0, h)

osnspb = LCP()

simu = TimeStepping(nsds, td, osi, osnspb)

# ==========================================================
# Storage
# ==========================================================

N = int(np.ceil((T - t0) / h))

data = np.zeros((N , 12))

q = ball.q()
v = ball.velocity()

qplane = movingplane.q()
vplane = movingplane.velocity()

lam = inter.lambda_python(1)
y = inter.y(0)

k = 0

data[k, 0] = t0
data[k, 1] = q[0]
data[k, 2] = v[0]
data[k, 7] = qplane[0]
data[k, 8] = vplane[0]

k += 1

# ==========================================================
# Time loop
# ==========================================================

while simu.hasNextEvent() and k < N:
    if k%100 == 0 :
        print("step:", k)
    simu.computeOneStep()

    data[k, 0] = simu.nextTime()

    data[k, 1] = q[0]
    data[k, 2] = v[0]

    data[k, 4] = lam[0]

    data[k, 7] = qplane[0]
    data[k, 8] = vplane[0]

    data[k, 11] = movingplane.reactionToBoundaryConditions()[0]
    simu.nextStep()

    k += 1

# ==========================================================
# Save
# ==========================================================

np.savetxt("BallOnMovingPlane_py.dat", data[:k])

import matplotlib.pyplot as plt

#
plt.subplot(411)
plt.title("position")
plt.plot(data[:, 0], data[:, 1], label='position ball')
plt.plot(data[:, 0], data[:, 7], label='position plane')
plt.legend()
plt.grid()
plt.subplot(412)
plt.title("velocity")
plt.plot(data[:, 0], data[:, 2], label='velocity ball')
plt.plot(data[:, 0], data[:, 8], label='velocity plane')
plt.grid()
plt.legend()
plt.subplot(413)
plt.title("lambda")
plt.plot(data[:, 0], data[:, 4], label='lambda')
plt.grid()
plt.legend()
plt.subplot(414)
plt.title("reaction to boundary conditions")
plt.plot(data[:, 0], data[:, 11], label='reaction to boundary conditions')
plt.grid()
plt.legend()
plt.show()
