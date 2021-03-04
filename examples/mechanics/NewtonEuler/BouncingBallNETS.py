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

import siconos.kernel as sk
from matplotlib.pyplot import subplot, title, plot, grid, show
import numpy as np
import math


class BouncingBallR(sk.NewtonEuler1DR):

    def __init__(self, ballRadius):
        self._ballRadius = ballRadius
        #sk.NewtonEuler1DR.__init__(self)
        super().__init__()

        
    def computeh(self, time, q, y):
        
        height = q[0] - self._ballRadius
        y[0] = height

        nnc = [1,0,0]
        self.setnc(nnc)

        ppc1 = [height, q[1], q[2]]
        self.setpc1(ppc1)
 
        ppc2 = [0.0, q[1], q[2]]
        self.setpc2(ppc2)

#
# dynamical system
#
x = [1.0, 0, 0, 1.0, 0, 0, 0]  # initial configuration
v = [2.0, 0, 0, 0, 0, 0]  # initial velocity
inertia = np.eye(3)       # inertia matrix

# the dynamical system
r = 0.1     # ball radius
g = 9.81    # gravity
m = 1      # ball mass
ball = sk.NewtonEulerDS(x, v, m, inertia)

# set external forces
weight = [-m * g, 0, 0]
ball.setFExtPtr(weight)

#
# Interactions
#

# ball-floor
e = 0.9     # restitution coeficient
nslaw = sk.NewtonImpactNSL(e)
relation = BouncingBallR(r)

inter = sk.Interaction(nslaw, relation)

#
# Model
#
t0 = 0      # start time
T = 10.0    # end time
bouncingBall = sk.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBall.insertDynamicalSystem(ball)

# link the interaction and the dynamical system
bouncingBall.link(inter, ball)

#
# Simulation
#
h = 0.005   # time step
theta = 0.5  # theta scheme

# (1) OneStepIntegrators
OSI = sk.MoreauJeanOSI(theta)

# (2) Time discretisation --
t = sk.TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = sk.LCP()

# (4) Simulation setup with (1) (2) (3)
s = sk.TimeStepping(bouncingBall, t, OSI, osnspb)
s.setNewtonTolerance(1e-10);
s.setNewtonMaxIteration(10);

# end of model definition

#
# computation
#

# the number of time steps
N = int((T - t0) // h) + 1

# Get the values to be plotted
# ->saved in a matrix dataPlot
dataPlot = np.empty((N+1, 16))

#
# numpy pointers on dense Siconos vectors
#
q = ball.q()
v = ball.twist()
p = ball.p(1)
lambda_ = inter.lambda_(1)

#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[0]
dataPlot[0, 2] = v[0]
dataPlot[0, 3] = p[0]
dataPlot[0, 4] = lambda_[0]
dataPlot[0, 5] = math.acos(q[3])
dataPlot[0, 6] = np.linalg.norm(relation.contactForce())
dataPlot[0, 7] = q[0]
dataPlot[0, 8] = q[1]
dataPlot[0, 9] = q[2]
dataPlot[0, 10] = q[3]
dataPlot[0, 11] = q[4]
dataPlot[0, 12] = q[5]
dataPlot[0, 13] = q[6]
dataPlot[0, 14] = v[1]
dataPlot[0, 15] = v[2]
k = 1

# # time loop
while(s.hasNextEvent()):
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0]
    dataPlot[k, 4] = lambda_[0]
    dataPlot[k, 5] = math.acos(q[3])
    dataPlot[k, 6] = np.linalg.norm(relation.contactForce())
    dataPlot[k, 7] = q[0]
    dataPlot[k, 8] = q[1]
    dataPlot[k, 9] = q[2]
    dataPlot[k, 10] = q[3]
    dataPlot[k, 11] = q[4]
    dataPlot[k, 12] = q[5]
    dataPlot[k, 13] = q[6]
    dataPlot[k, 14] = v[1]
    dataPlot[k, 15] = v[2]
    k = k + 1
    s.nextStep()

np.savetxt("result-py.dat", dataPlot)

#
# comparison with the reference file
#
sk.compareRefFile(dataPlot, "BouncingBallNETS.ref", 1e-12)

#
# plots
#
subplot(411)
title('position')
plot(dataPlot[:, 0], dataPlot[:, 1])
grid()
subplot(412)
title('velocity')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
subplot(413)
plot(dataPlot[:, 0], dataPlot[:, 3])
title('reaction')
grid()
subplot(414)
plot(dataPlot[:, 0], dataPlot[:, 4])
title('lambda')
grid()
show()
