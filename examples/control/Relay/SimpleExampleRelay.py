#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2018 INRIA.
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

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig, figure
from numpy import array, eye, empty, zeros, savetxt, loadtxt, linalg
from siconos.kernel import FirstOrderLinearDS,  FirstOrderLinearR, RelayNSL,\
    NonSmoothDynamicalSystem, TimeDiscretisation, TimeStepping, EulerMoreauOSI, \
    Interaction, Relay
from math import ceil


# variables
t0 = 0.0   # start time
T = 1.0     # end time
h = 1.0e-3   # time step
Vinit =1.0
theta = 0.5
N = int((T-t0)/h)

# matrices
A = zeros((2,2))
x0 = array([Vinit,-Vinit])
B = 2.0*eye(2)
C = eye(2)
D = zeros((2,2))

# dynamical systems
process = FirstOrderLinearDS(x0, A)

myNslaw = RelayNSL(2)
myProcessRelation=FirstOrderLinearR(C, B);
#myProcessRelation.setDPtr(D)

myProcessInteraction = Interaction(myNslaw,
                                   myProcessRelation)

simplerelay = NonSmoothDynamicalSystem(t0,T)
simplerelay.insertDynamicalSystem(process)
simplerelay.link(myProcessInteraction,process)


#myProcessRelation.computeJachx(0, x0, x0 , x0, C)

td = TimeDiscretisation(t0, h)
s = TimeStepping(simplerelay, td)

myIntegrator = EulerMoreauOSI(theta)

s.insertIntegrator(myIntegrator)

osnspb = Relay()
s.insertNonSmoothProblem(osnspb)
#s.setComputeResiduY(True)
#s.setComputeResiduR(True)

# matrix to save data
dataPlot = empty((N+1,7))
k=0
dataPlot[k, 0] = t0
dataPlot[k, 1:3] = process.x()
dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
dataPlot[k, 4] = myProcessInteraction.lambda_(0)[1]
dataPlot[k, 5] = myProcessInteraction.y(0)[0]
dataPlot[k, 6] = myProcessInteraction.y(0)[1]
# time loop
k = 1
print('start computation')
while(s.hasNextEvent()):
     if not (k%50):
         sys.stdout.write('.')

     s.computeOneStep()
     #osnspb.display()
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     dataPlot[k, 4] = myProcessInteraction.lambda_(0)[1]
     dataPlot[k, 5] = myProcessInteraction.y(0)[0]
     dataPlot[k, 6] = myProcessInteraction.y(0)[1]
     k += 1
     s.nextStep()
     #print s.nextTime()
sys.stdout.write('\n')
# save to disk
savetxt('SimpleExampleRelay_py.dat', dataPlot)

dataRef = loadtxt('SimpleExampleRelay_py.ref')

print('Comparison with reference file -  error = ',linalg.norm(dataPlot-dataRef) )



assert(linalg.norm(dataPlot-dataRef) <= 1e-12)


# plot interesting stuff
subplot(311)
title('x_1')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(312)
title('x_2')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
subplot(313)
plot(dataPlot[:,0], dataPlot[:,3])
title('lambda')
grid()
savefig("SimpleRelay_py.png")
