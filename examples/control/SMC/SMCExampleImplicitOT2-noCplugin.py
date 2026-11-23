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

import sys
from math import ceil, sin

import numpy as np
import siconos.plot_config as sicoplot
from numpy import eye, savetxt, zeros
from numpy.linalg import norm
from siconos.control.controller import LinearSMCOT2
from siconos.control.sensor import LinearSensor
from siconos.control.simulation import ControlZOHSimulation
from siconos.kernel import FirstOrderLinearDS, SiconosMatrix, getMatrix

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)


# Derive our own version of FirstOrderLinearDS
class MyFOLDS(FirstOrderLinearDS):
    """derived FirstOrderLinearDS class to show how to override a method"""

    def computeb(self, time):
        t = sin(50 * time)
        u = [t, -t]
        self.setbPtr(u)


# variable declaration
ndof = 2  # Number of degrees of freedom of your system
t0 = 0.0  # start time
T = 1  # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2  # time step for control
Xinit = 1.0  # initial position
theta = 0.5
N = 2 * ceil((T - t0) / h)  # number of time steps
outputSize = 5  # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = [[0, 1.0]]
Brel = [[0], [2]]

# Simple check
if h > hControl:
    print("hControl must be bigger than h")
    sys.exit(1)

# Declaration of the Dynamical System
processDS = MyFOLDS(x0, A)
# XXX b is not automatically created ...
processDS.setbPtr([0, 0])
# Control simulation
sim = ControlZOHSimulation(t0, T, h)
sim.setSaveOnlyMainSimulation(True)
sim.addDynamicalSystem(processDS)
# Actuator, Sensor & ControlManager
sens = LinearSensor(processDS, sensorC)
sim.addSensor(sens, hControl)
act = LinearSMCOT2(sens)
act.setCsurface(Csurface)
act.setB(Brel)
sim.addActuator(act, hControl)

# Initialization
sim.initialize()

# Run simulation
sim.run()

# Get data
dataPlot = sim.data()

# Save to disk
savetxt("SMCExampleImplicitOT2-noCplugin-py.dat", dataPlot)
# Plot interesting data
plt.subplot(411)
plt.title("x1")
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title("x2")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.title("u")
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.grid()
plt.savefig("ismcOT2-noCplugin.png")

# compare with the reference
ref = getMatrix(SiconosMatrix("SMCExampleImplicitOT2-py.ref"))
print(f"{norm(dataPlot - ref):e}")
assert np.allclose(dataPlot, ref, atol=1e-11)
