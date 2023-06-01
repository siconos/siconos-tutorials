#!/usr/bin/env python

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


import  numpy as np

from numpy.linalg import norm
from siconos.kernel import FirstOrderLinearTIDS, FirstOrderLinearDS, RelayNSL, Interaction, FirstOrderLinearR, NonSmoothDynamicalSystem,\
    EulerMoreauOSI, TimeDiscretisation, TimeStepping, Relay




# LagrangianLinearTIDS, NewtonImpactNSL,\
#     LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem, MoreauJeanOSI,\
#     TimeDiscretisation, LCP, TimeStepping
from siconos.kernel import SimpleMatrix, getMatrix

import siconos.numerics as sn


from numpy import eye, empty, float64, zeros

"""

Equations of motion

m1 x1" +  k1  x1 =  r(t)
Coulomb(x1'(t),r(t))=0 

"""


t0 = 0       # start time
T = 100     # end time
h = 1e-02    # time step
theta = 0.5  # theta scheme

m = 1; stiffness = 1;
alpha = 1
xinit=12.
vinit=6.

x0 = [xinit, vinit]    # initial state

A = np.zeros((2,2))

A[0,1] = 1
A[1,0] = -stiffness/m

#
# dynamical system
#

oscillator = FirstOrderLinearTIDS(x0, A)
#oscillator.display()

#oscillator = FirstOrderLinearDS(x0, A)
#oscillator.setComputebFunction("plugins", "TwoDofsOscillatorB")

#
# Interactions
#


B = np.array([[0] , [alpha]])
C = np.array([[0  , 1.]])

nslaw = RelayNSL(1, -1, 1)

relation = FirstOrderLinearR(C, B)
inter = Interaction(nslaw, relation)

inter.display()

#
# Model
#
frictionOscillator = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
frictionOscillator.insertDynamicalSystem(oscillator)

# link the interaction and the dynamical system
frictionOscillator.link(inter, oscillator)



#
# Simulation
#

# (1) OneStepIntegrators
OSI = EulerMoreauOSI(theta)

# (2) Time discretisation --
t = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = Relay()
osnspb.setSolverId(sn.SICONOS_RELAY_LEMKE);

osnspb.numericsSolverOptions().dparam[0] = 1e-08;
# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(frictionOscillator,t, OSI, osnspb)

# end of model definition

#
# computation
#
import math

# the number of time steps
N = math.ceil((T - t0) / h)

# # Get the values to be plotted
# # ->saved in a matrix dataPlot

dataPlot = zeros((N+1, 6))

# #
# # numpy pointers on dense Siconos vectors
# #
x = oscillator.x()
# p = ball.p(1)
lambda_ = inter.lambda_(0)
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
    #print('.')
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
ref = getMatrix(SimpleMatrix("FrictionOscillator.ref"))
if (norm(dataPlot[0:10000,0:6] - ref[0:10000,0:6]) > 1e-12):
    print("Warning. The result is rather different from the reference file.", norm(dataPlot[0:10000,0:6] - ref[0:10000,0:6]) )
else:
    print("Error with the reference file.", norm(dataPlot[0:10000,0:6] - ref[0:10000,0:6]) )

#
# plots
#
import matplotlib,os
havedisplay = "DISPLAY" in os.environ
if not havedisplay:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.subplot(411)
plt.title('position x1')
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title('velocity dot x1')
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.figure()
plt.plot(dataPlot[:, 1], dataPlot[:, 2])
plt.title('lambda')


plt.figure()
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title('lambda')


if havedisplay:
    plt.show()
else:
    plt.savefig("bbts.png")
