# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
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
import numpy.linalg as LA
import siconos.geometry as sg
import siconos.integrators
import siconos.modeling as sm
import siconos.nonsmooth_formulations
import siconos.plot_config as sicoplot
import siconos.simulation

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)


t0 = 0.0  # start time
h = 0.001  # time step
N = 10000
T = 4.61  # h * N
theta = 0.5  # theta scheme

#
# dynamical system
#
# initial configuration
initial_position = np.zeros(7, dtype=np.float64)
initial_position[3] = 1.0
initial_twist = np.zeros(6, dtype=np.float64)
inertia = np.zeros((3, 3), dtype=np.float64, order="F")
inertia[0, 0] = 5.0
inertia[1, 1] = 10.0
inertia[2, 2] = 1.0
mass = 1.0


def compute_mext(time, mExt):
    td = 2.0 - h
    mExt[:] = 0
    if 0 <= time < td:
        mExt[0] = 20.0
    elif td <= time <= td + h:
        mExt[1] = 1.0 / (5.0 * h)


unstableRotation = sm.NewtonEulerDS(
    initial_position, initial_twist, mass, inertia, sm.alias_t
)

unstableRotation.setComputeMextFunction(compute_mext)
unstableRotation.setIsMextExpressedInInertialFrame(True)

rotationVector_init = np.zeros(3, dtype=np.float64)
rotationVector_init[0] = 0.3
initial_twist_heavy_top = np.asarray([0, 0, 0, 0, 0, 50], dtype=np.float64)
initial_position_heavy_top = sg.quaternionFromRotationVector(
    np.asarray([0.3, 0, 0])
)  # rotationVector_init)
inertia_heavytop = np.zeros((3, 3), dtype=np.float64, order="F")
inertia_heavytop[0, 0] = 5.0
inertia_heavytop[1, 1] = 5.0
inertia_heavytop[2, 2] = 1.0
mass_heavytop = 1.0
Mg = 20
length = 1.0


def centermass(q):
    r = np.zeros(3)
    E3 = np.zeros(3)
    E3[2] = 1.0
    sg.rotateVectorFromInertialToBodyFrame(q, E3)
    r[0] = E3[0]
    r[1] = E3[1]
    r[2] = E3[2]
    return r


def compute_mint(twist, pos, time, mint):
    r = centermass(pos)
    mint[...] = np.asarray(Mg * length * np.cross(r, [0, 0, 1.0]))
    sg.rotateVectorFromInertialToBodyFrame(pos, mint)
    print("mInt", mint)


heavytop = sm.NewtonEulerDS(
    initial_position_heavy_top,
    initial_twist_heavy_top,
    mass_heavytop,
    inertia_heavytop,
    sm.alias_t,
)
heavytop.setComputeJacobianMintOver_q_byFD(True)


ds = unstableRotation

# ds = heavytop

print(ds)


# Non-Smooth Dynamical System
#
nsds = sm.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
nsds.insertDynamicalSystem(ds)

#
# Simulation
#

# (1) OneStepIntegrators
OSI = siconos.integrators.MoreauJeanOSI(theta)

# (2) Time discretisation --
t = siconos.simulation.TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = siconos.nonsmooth_formulations.LCP()

# (4) Simulation setup with (1) (2) (3)
s = siconos.simulation.TimeStepping(nsds, t, OSI, osnspb)
# s.setDisplayNewtonConvergence(True)
s.setNewtonTolerance(1e-10)
# s.setNewtonMaxIteration(1)

# end of model definition

#
# computation
#

# Get the values to be plotted
# ->saved in a matrix dataPlot
dataPlot = np.empty((N + 1, 25))

#
# numpy pointers on dense Siconos vectors
#
q = ds.q()
v = ds.twist()
p = ds.p(1)

#
# initial data
#
k = 0
dataPlot[k, 1] = q[0]
dataPlot[k, 2] = q[1]
dataPlot[k, 3] = q[2]
dataPlot[k, 4] = q[3]
dataPlot[k, 5] = q[4]
dataPlot[k, 6] = q[5]
dataPlot[k, 7] = q[6]

dataPlot[k, 8] = v[0]
dataPlot[k, 9] = v[1]
dataPlot[k, 10] = v[2]
dataPlot[k, 11] = v[3]
dataPlot[k, 12] = v[4]
dataPlot[k, 13] = v[5]

omega = v[3:6]
inertia = ds.totalInertiaMatrix
angular_momentum = np.dot(inertia[3:6, 3:6], omega)
sg.rotateVectorFromBodyToInertialFrame(q, angular_momentum)

dataPlot[k, 14] = angular_momentum[0]
dataPlot[k, 15] = angular_momentum[1]
dataPlot[k, 16] = angular_momentum[2]
dataPlot[k, 17] = LA.norm(angular_momentum)

rotationVector = np.zeros(3)
rotationVector = sg.rotationVectorFromQuaternion(q[3], q[4], q[5], q[6])
dataPlot[k, 18] = rotationVector[0]
dataPlot[k, 19] = rotationVector[1]
dataPlot[k, 20] = rotationVector[2]


dataPlot[k, 21] = h * omega[0]
dataPlot[k, 22] = h * omega[1]
dataPlot[k, 23] = h * omega[2]
dataPlot[k, 24] = LA.norm(h * omega)


k = 1

# time loop
while s.hasNextEvent() and k < N:
    # print(' ' )
    # print (
    #     '------- k = ',
    #     k,
    #     '-----------------------------------------')
    # print(' ' )
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = q[1]
    dataPlot[k, 3] = q[2]
    dataPlot[k, 4] = q[3]
    dataPlot[k, 5] = q[4]
    dataPlot[k, 6] = q[5]
    dataPlot[k, 7] = q[6]

    dataPlot[k, 8] = v[0]
    dataPlot[k, 9] = v[1]
    dataPlot[k, 10] = v[2]
    dataPlot[k, 11] = v[3]
    dataPlot[k, 12] = v[4]
    dataPlot[k, 13] = v[5]

    omega = v[3:6]
    inertia = ds.totalInertiaMatrix
    angular_momentum = np.dot(inertia[3:6, 3:6], omega)
    sg.rotateVectorFromBodyToInertialFrame(q, angular_momentum)
    a = np.zeros(1)
    a[0] = angular_momentum[0]
    # a[1] = am(1)
    # print "omega", omega
    # print "angular_momentum", angular_momentum,
    # print "q=", q
    # print " norm(a[1:2])", np.linalg.norm(a)
    # raw_input()
    dataPlot[k, 14] = angular_momentum[0]
    dataPlot[k, 15] = angular_momentum[1]
    dataPlot[k, 16] = angular_momentum[2]
    dataPlot[k, 17] = LA.norm(angular_momentum)
    rotationVector = sg.rotationVectorFromQuaternion(q[3], q[4], q[5], q[6])

    dataPlot[k, 18] = rotationVector[0]
    dataPlot[k, 19] = rotationVector[1]
    dataPlot[k, 20] = rotationVector[2]

    dataPlot[k, 21] = h * omega[0]
    dataPlot[k, 22] = h * omega[1]
    dataPlot[k, 23] = h * omega[2]
    dataPlot[k, 24] = LA.norm(h * omega)

    k = k + 1
    s.nextStep()


dataPlot = np.resize(dataPlot, (k - 2, 25))


np.savetxt("result-py.dat", dataPlot)

ref = np.loadtxt("result-py.ref")

assert np.allclose(ref, dataPlot)

plt.figure(num="Moreau Jean Siconos", figsize=(12, 12))
plt.subplot(321)
plt.title("angular velocities Omega")
plt.plot(dataPlot[:, 0], dataPlot[:, 11])
plt.plot(dataPlot[:, 0], dataPlot[:, 12])
# plt.plot(dataPlot[:, 0], dataPlot[:, 13])

plt.subplot(322)
plt.title("rotation vector")
plt.plot(dataPlot[:, 0], dataPlot[:, 18])
plt.plot(dataPlot[:, 0], dataPlot[:, 19])
plt.plot(dataPlot[:, 0], dataPlot[:, 20])

plt.subplot(323)
plt.title("Theta (h Omega)")
plt.plot(dataPlot[:, 0], dataPlot[:, 21])
plt.plot(dataPlot[:, 0], dataPlot[:, 22])
plt.plot(dataPlot[:, 0], dataPlot[:, 23])

plt.subplot(325)
plt.title("norm of Theta")
plt.plot(dataPlot[:, 0], dataPlot[:, 24])

plt.subplot(324)
plt.title("angular momentum (pi[0])")
plt.plot(dataPlot[:, 0], dataPlot[:, 14])
# plt.plot(dataPlot[:, 0], dataPlot[:, 15])
# plt.plot(dataPlot[:, 0], dataPlot[:, 16])

plt.subplot(326)
plt.title("norm of angular momentum  pi")
plt.plot(dataPlot[:, 0], dataPlot[:, 17])

plt.grid()
if enable_plot:
    plt.show()
