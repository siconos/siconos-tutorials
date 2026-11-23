#!/usr/bin/env python

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


import math

import numpy as np
import siconos.integrators as si
import siconos.modeling as sm
import siconos.nonsmooth_formulations as snsf
import siconos.numerics as sn
import siconos.plot_config as sicoplot
import siconos.simulation as ss
from numpy.linalg import norm

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

"""

Equations of motion

m1 x1" +  k1  x1 =  r(t)
Coulomb(x1'(t),r(t))=0

"""


t0 = 0  # start time
T = 100  # end time
h = 1e-02  # time step
theta = 0.5  # theta scheme

m = 1
stiffness = 1
alpha = 1
xinit = 12.0
vinit = 6.0

x0 = np.array([xinit, vinit], dtype=np.float64)  # initial state

A = np.zeros((2, 2), dtype=np.float64, order="F")

A[0, 1] = 1
A[1, 0] = -stiffness / m

#
# dynamical system
#

oscillator = sm.FirstOrderLinearDS(x0, sm.alias_t)
oscillator.setConstantA(A, sm.alias_t)
#
# Interactions
#


B = np.array([[0], [alpha]], dtype=np.float64, order="F")
C = np.array([[0, 1.0]], dtype=np.float64, order="F")

nslaw = sm.RelayNSL(1, -1, 1)

# relation = sm.FirstOrderLinearR(C, B)
relation = sm.FirstOrderLinearR()
relation.setConstantBAlias(B)
relation.setConstantCAlias(C)
inter = sm.Interaction(nslaw, relation)

print(inter)

#
# Model
#
frictionOscillator = sm.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
frictionOscillator.insertDynamicalSystem(oscillator)

# link the interaction and the dynamical system
frictionOscillator.link(inter, oscillator)


#
# Simulation
#

# (1) OneStepIntegrators
OSI = si.EulerMoreauOSI(theta)

# (2) Time discretisation --
t = ss.TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = snsf.Relay()
osnspb.setSolverId(sn.solver_ids.SICONOS_RELAY_LEMKE)

osnspb.numericsSolverOptions().dparam[0] = 1e-08
# (4) Simulation setup with (1) (2) (3)
s = ss.TimeStepping(frictionOscillator, t, OSI, osnspb)

# end of model definition

#
# computation
#
# the number of time steps
N = math.ceil((T - t0) / h)

# # Get the values to be plotted
# # ->saved in a matrix dataPlot

dataPlot = np.zeros((N + 1, 6))

# #
# # numpy pointers on dense Siconos vectors
# #
x = oscillator.x()
# p = ball.p(1)
lambda_ = inter.lambda_python(0)
y = inter.y(0)


# #
# # initial data
# #
dataPlot[0, 0] = t0
dataPlot[0, 1] = x[0]
dataPlot[0, 2] = x[1]
dataPlot[0, 3] = lambda_[0]
dataPlot[0, 4] = y[0]


k = 1

# time loop
while s.hasNextEvent():
    s.computeOneStep()
    # print('.')
    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]
    dataPlot[k, 3] = lambda_[0]
    dataPlot[k, 4] = y[0]
    k += 1
    s.nextStep()

#
# comparison with the reference file (produced by the cpp version)
#
ref = np.loadtxt("FrictionOscillator.ref", skiprows=1)
if norm(dataPlot[0:10000, 0:6] - ref[0:10000, 0:6]) > 1e-12:
    print(
        "Warning. The result is rather different from the reference file.",
        norm(dataPlot[0:10000, 0:6] - ref[0:10000, 0:6]),
    )
else:
    print(
        "Error with the reference file.",
        norm(dataPlot[0:10000, 0:6] - ref[0:10000, 0:6]),
    )

#
# plots
#
plt.subplot(411)
plt.title("position x1")
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title("velocity dot x1")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.figure()
plt.plot(dataPlot[:, 1], dataPlot[:, 2])
plt.title("lambda")


plt.figure()
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title("lambda")


if enable_plot:
    plt.show()
else:
    plt.savefig("bbts.png")
