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
import math

import numpy as np
import siconos.integrators
import siconos.plot_config as sicoplot
from siconos.modeling import (
    FirstOrderLinearDS,
    FirstOrderLinearTIR,
    Interaction,
    NonSmoothDynamicalSystem,
    RelayNSL,
)

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0.0
T = 100  # Total simulation time
h_step = 1.0e-2  # Time step
xinit = 3.0 * math.sqrt(2.0) / (2.0 * math.pi)  # initial voltage
Modeltitle = "RelayOscillator"


#
# dynamical system
#
init_state = [0.0, xinit, 0]

A = [[0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, -3.0, -2.0]]

LSRelayOscillator = FirstOrderLinearDS(init_state, A)

#
# Interactions
#

C = [[1.0, 0.0, 0.0]]

D = [[0.0]]

B = [[0.0], [0.0], [1.0]]

LTIRRelayOscillator = FirstOrderLinearTIR(C, B)
LTIRRelayOscillator.setConstantD(D)

nslaw = RelayNSL(1)
InterRelayOscillator = Interaction(nslaw, LTIRRelayOscillator)


#
# Model
#
relayOscillator = NonSmoothDynamicalSystem(t0, T)
relayOscillator.setTitle(Modeltitle)
#   add the dynamical system in the non smooth dynamical system
relayOscillator.insertDynamicalSystem(LSRelayOscillator)


#   link the interaction and the dynamical system
relayOscillator.link(InterRelayOscillator, LSRelayOscillator)


#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
aOSI = siconos.integrators.EulerMoreauOSI(theta)
# (2) Time discretisation
aTiDisc = siconos.simulation.TimeDiscretisation(t0, h_step)

# (3) Non smooth problem
aRelay = siconos.nonsmooth_formulations.Relay()

# (4) Simulation setup with (1) (2) (3)
aTS = siconos.simulation.TimeStepping(relayOscillator, aTiDisc, aOSI, aRelay)

# end of model definition

#
# computation
#


k = 0
h = aTS.timeStep()
print("Timestep : ", h)
# Number of time steps
N = (int)((T - t0) / h)
print("Number of steps : ", N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

dataPlot = np.empty([N + 1, 8])

x = LSRelayOscillator.x()
print("Initial state : ", x)
y = InterRelayOscillator.y(0)
print("First y : ", y)
lambda_ = InterRelayOscillator.lambda_(0)

while k < N:
    aTS.computeOneStep()
    # aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]
    dataPlot[k, 3] = x[2]
    dataPlot[k, 4] = y[0]
    dataPlot[k, 5] = lambda_[0]

    k += 1
    if k % 1000 == 0:
        print("step =", k, " < ", N)
    aTS.nextStep()

if enable_plot:
    #
    # plots
    #
    plt.subplot(511)
    plt.title("x1")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 1])
    plt.grid()
    plt.subplot(512)
    plt.title("x2")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 2])
    plt.grid()
    plt.subplot(513)
    plt.title("x3")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 3])
    plt.subplot(514)
    plt.title("y")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 4])
    plt.grid()
    plt.subplot(515)
    plt.title("lambda")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 5])
    plt.grid()
    plt.savefig("relay_oscillator.png")
