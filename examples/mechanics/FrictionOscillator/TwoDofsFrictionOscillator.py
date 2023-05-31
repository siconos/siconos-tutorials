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


""" siconos -v TwoDofsFrictionOscillator.py --build-plugins """


import  numpy as np

from numpy.linalg import norm
from siconos.kernel import FirstOrderLinearTIDS, FirstOrderLinearDS, RelayNSL, Interaction, FirstOrderLinearR, NonSmoothDynamicalSystem,\
    EulerMoreauOSI, TimeDiscretisation, TimeStepping, Relay

import siconos.numerics as sn


"""

Equations of motion

m1 x1" + d1*x1' + (k1+k2) x1 - k2 x2 = f1*cos(ωt)
m2 x2" + d2*x2' + k2 x2 - k2 x1 = f2*cos(ωt)+r(t)
Coulomb(x2'(t),r(t))=0 

With Coulomb(x2',r) = max(μ*N,|r-ρ*x2'|)*r-μ*N*(r-ρ*x2').

numerical values :

m1 = 1; k1 = 1;
m2 = 1; k2 = 1;
μ = 0.5; N = 1;
d1 = 0; d2 = 0;
f1 = 1;
f2 = 0;
ω = 0.3;
ρ = 0.1;
"""


t0 = 0       # start time
T = 500       # end time
h = 5e-02    # time step
theta = 0.5  # theta scheme

m1 = 1; k1 = 1;
m2 = 1; k2 = 1;
μ = 0.5; N = 1;
d1 = 0; d2 = 0;
f1 = 1;
f2 = 0;
ω = 0.3;


x_1_0=1.0
dot_x_1_0=0.0
x_2_0=0.0
dot_x_2_0=0.0


x = [x_1_0, dot_x_1_0, x_2_0, dot_x_2_0]    # initial state

M = np.zeros((4,4))

M[0,0]=1
M[1,1]=m1
M[2,2]=1
M[3,3]=m2



A = np.zeros((4,4))

A[0,1] = 1

A[1,0] = -(k1+k2)
A[1,1] = - d1
A[1,2] = k2


A[2,3] = 1

A[3,0] =  k2
A[3,2] = -k2
A[3,3] = - d2



#
# dynamical system
#

oscillator = FirstOrderLinearDS(x, A)
oscillator.setMPtr(M)

oscillator.setComputebFunction("plugins", "TwoDofsOscillatorB")
#oscillator.display()

#
# Interactions
#


H = np.array([[0, 0, 0, 1]])

nslaw = RelayNSL(1, -μ*N, μ*N)
relation = FirstOrderLinearR(H, H.T)
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


# the number of time steps
N = int((T - t0) / h)

# # Get the values to be plotted
# # ->saved in a matrix dataPlot

dataPlot = np.zeros((N+1, 6))

# #
# # numpy pointers on dense Siconos vectors
# #
x = oscillator.x()
# p = ball.p(1)
lambda_ = inter.lambda_(0)


# #
# # initial data
# #
dataPlot[0, 0] = t0
dataPlot[0, 1] = x[0]
dataPlot[0, 2] = x[1]
dataPlot[0, 3] = x[2]
dataPlot[0, 4] = x[3]
dataPlot[0, 5] = lambda_[0]


k = 1

# time loop
while s.hasNextEvent():
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]
    dataPlot[k, 3] = x[2]
    dataPlot[k, 4] = x[3]
    dataPlot[k, 5] = lambda_[0]
    
    k += 1
    s.nextStep()




    
# #
# # comparison with the reference file
# #
# ref = getMatrix(SimpleMatrix("BouncingBallTS.ref"))

# if (norm(dataPlot - ref) > 1e-12):
#     print("Warning. The result is rather different from the reference file.")


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
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title('position x2')
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title('velocity dot x2')
plt.grid()


plt.figure()
plt.plot(dataPlot[:, 0], dataPlot[:, 5])
plt.title('lambda')


if havedisplay:
    plt.show()
else:
    plt.savefig("bbts.png")
