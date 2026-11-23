#!/usr/bin/env python
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
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
# -----------------------------------------------------------------------
#
#  DiodeBridge  : sample of an electrical circuit involving :
#    - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
#    - a non smooth system (a 1000 Ohm resistor supplied through a 4
# 	diodes bridge) in parallel with the oscillator
#
#  Expected behavior :
#
#  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
#  an initial energy.
#  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
#  The non smooth system is a full wave rectifier :
#  each phase (positive and negative) of the oscillation allows current to flow
#  through the resistor in a constant direction, resulting in an energy loss :
#  the oscillation damps.
#
#  State variables :
#    - the voltage across the capacitor (or inductor)
#    - the current through the inductor
#
#  Since there is only one dynamical system, the interaction is defined by :
#    - complementarity laws between diodes current and
#  voltage. Depending on the diode position in the bridge, y stands
#  for the reverse voltage across the diode or for the diode
#  current (see figure in the template file)
#    - a linear time invariant relation between the state variables and
# 	y and lambda (derived from Kirchhoff laws)
#
# -----------------------------------------------------------------------

import numpy as np
import siconos.integrators
import siconos.modeling as sm
import siconos.nonsmooth_formulations
import siconos.plot_config as sicoplot
import siconos.simulation
from numpy import linalg as LA

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0.0
T = 5.0e-3  # Total simulation time
h_step = 1.0e-6  # Time step
Lvalue = 1e-2  # inductance
Cvalue = 1e-6  # capacitance
Rvalue = 1e3  # resistance
Vinit = 10.0  # initial voltage
Modeltitle = "DiodeBridge"

#
# dynamical system
#
init_state = np.asarray([Vinit, 0], dtype=np.float64)

A = np.asarray([[0, -1.0 / Cvalue], [1.0 / Lvalue, 0]], dtype=np.float64, order="F")

LSDiodeBridge = sm.FirstOrderLinearDS(init_state, sm.alias_t)
LSDiodeBridge.setConstantA(A, sm.alias_t)

#
# Interactions
#

C = np.asarray(
    [[0.0, 0.0], [0, 0.0], [-1.0, 0.0], [1.0, 0.0]], order="F", dtype=np.float64
)

D = np.asarray(
    [
        [1.0 / Rvalue, 1.0 / Rvalue, -1.0, 0.0],
        [1.0 / Rvalue, 1.0 / Rvalue, 0.0, -1.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ],
    order="F",
    dtype=np.float64,
)

B = np.asarray(
    [[0.0, 0.0, -1.0 / Cvalue, 1.0 / Cvalue], [0.0, 0.0, 0.0, 0.0]],
    order="F",
    dtype=np.float64,
)

LTIRDiodeBridge = sm.FirstOrderLinearTIR(C, B)
LTIRDiodeBridge.setConstantD(D)

nslaw = sm.ComplementarityConditionNSL(4)
InterDiodeBridge = sm.Interaction(nslaw, LTIRDiodeBridge)


#
# Model
#
DiodeBridge = sm.NonSmoothDynamicalSystem(t0, T)
DiodeBridge.setTitle(Modeltitle)

#   add the dynamical system in the non smooth dynamical system
DiodeBridge.insertDynamicalSystem(LSDiodeBridge)

#   link the interaction and the dynamical system
DiodeBridge.link(InterDiodeBridge, LSDiodeBridge)

#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
gamma = 0.5
aOSI = siconos.integrators.EulerMoreauOSI(theta, gamma)
aOSI.useGammaForRelation = True

# (2) Time discretisation
aTiDisc = siconos.simulation.TimeDiscretisation(t0, h_step)

# (3) Non smooth problem
aLCP = siconos.nonsmooth_formulations.LCP()

# (4) Simulation setup with (1) (2) (3)
aTS = siconos.simulation.TimeStepping(DiodeBridge, aTiDisc, aOSI, aLCP)

# end of model definition

#
# computation
#

k = 0
h = aTS.timeStep()
print("Timestep : ", h)
# Number of time steps
N = int((T - t0) / h)
print("Number of steps : ", N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

dataPlot = np.zeros([N, 10])

x = LSDiodeBridge.x()
print("Initial state : ", x)
y = InterDiodeBridge.y(0)
print("First y : ", y)
lambda_ = InterDiodeBridge.lambda_python(0)

# For the initial time step:
# time

#  inductor voltage
dataPlot[k, 1] = x[0]

# inductor current
dataPlot[k, 2] = x[1]

# diode R1 current
dataPlot[k, 3] = y[0]

# diode R1 voltage
dataPlot[k, 4] = -lambda_[0]

# diode F2 voltage
dataPlot[k, 5] = -lambda_[1]

# diode F1 current
dataPlot[k, 6] = lambda_[2]

k += 1
while k < N:
    aTS.computeOneStep()
    # aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    # inductor current
    dataPlot[k, 2] = x[1]
    # diode R1 current
    dataPlot[k, 3] = y[0]
    # diode R1 voltage
    dataPlot[k, 4] = -lambda_[0]
    # diode F2 voltage
    dataPlot[k, 5] = -lambda_[1]
    # diode F1 current
    dataPlot[k, 6] = lambda_[2]
    k += 1
    aTS.nextStep()


ref = np.loadtxt("DiodeBridge.ref", skiprows=1)

error = LA.norm(dataPlot[:, 0:6] - ref[:, 0:6])

print("error = ", error)

# assert (error < 1e-09)
withRef = True
if enable_plot:
    #
    # plots
    #
    plt.subplot(411)
    plt.title("inductor voltage")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 1])
    if withRef:
        plt.plot(ref[0 : k - 1, 0], ref[0 : k - 1, 1])
    plt.grid()
    plt.subplot(412)
    plt.title("inductor current")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 2])
    if withRef:
        plt.plot(ref[0 : k - 1, 0], ref[0 : k - 1, 2])
    plt.grid()
    plt.subplot(413)
    plt.title("diode R1 (blue) and F2 (green) voltage")
    plt.plot(dataPlot[0 : k - 1, 0], -dataPlot[0 : k - 1, 4])
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 5])
    if withRef:
        plt.plot(ref[0 : k - 1, 0], -ref[0 : k - 1, 4])
        plt.plot(ref[0 : k - 1, 0], ref[0 : k - 1, 5])
    plt.grid()
    plt.subplot(414)
    plt.title("resistor current")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 3] + dataPlot[0 : k - 1, 6])
    if withRef:
        plt.plot(dataPlot[0 : k - 1, 0], ref[0 : k - 1, 3] + ref[0 : k - 1, 6])
    plt.grid()
    plt.savefig("diode_brige_tgs.png")
