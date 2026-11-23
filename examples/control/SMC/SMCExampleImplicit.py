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

import distutils.spawn
from math import ceil

import scipy
from matplotlib import rc
from numpy import eye, savetxt, zeros
from scipy import arange
from siconos.control.controller import LinearSMC
from siconos.control.sensor import LinearSensor
from siconos.control.simulation import ControlZOHSimulation
from siconos.kernel import FirstOrderLinearDS

if distutils.spawn.find_executable("latex"):
    rc("text", usetex=True)
import sys

import siconos.plot_config as sicoplot

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

# variable declaration
ndof = 2  # Number of degrees of freedom of your system
t0 = 0.0  # start time
T = 1  # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2  # time step for control
Xinit = 1.0  # initial position
theta = 0.5
N = ceil((T - t0) / h + 10)  # number of time steps
outputSize = 5  # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = [[0, 1]]
Brel = [[0], [2]]
# Drel = [[0, 0]]
# Simple check
if h > hControl:
    print("hControl must be bigger than h")
    sys.exit(1)

# Declaration of the Dynamical System
processDS = FirstOrderLinearDS(x0, A)
processDS.setComputebFunction("RelayPlugin", "computeB")

# Control simulation
sim = ControlZOHSimulation(t0, T, h)
sim.setSaveOnlyMainSimulation(True)
sim.addDynamicalSystem(processDS)
# Actuator, Sensor & ControlManager
sens = LinearSensor(processDS, sensorC)
sim.addSensor(sens, hControl)
act = LinearSMC(sens)
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
savetxt("SMCExampleImplicit-py.dat", dataPlot)
# Plot interesting data

plt.subplot(211)
plt.ylabel(r"$\sigma$")
plt.xlabel(r"t")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(212)
plt.ylabel(r"$\bar{u}^s$")
plt.xlabel(r"t")
plt.ylim(-2.1, 2.1)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.savefig("ismc_sigma_u.png")

plt.subplot(211)
plt.ylabel(r"$\sigma$")
plt.xlabel(r"t")
plt.xlim(xmin=0.49)
plt.ylim(-0.03, 0.03)
plt.plot(dataPlot[4900:, 0], dataPlot[4900:, 2])
plt.grid()
plt.subplot(212)
plt.ylabel(r"$\bar{u}^s$")
plt.xlabel(r"t")
plt.ylim(-2.1, 2.1)
plt.xlim(xmin=0.49)
p1 = plt.plot(dataPlot[4900:, 0], dataPlot[4900:, 3])
# p2 = plt.plot(dataPlot[4900:, 0], np.sin(50*dataPlot[4900:, 0]))
# plt.legend((p1[0], p2[0]), (r'$\bar{u}^s(t)$', r'$-\rho(t)$'), ncol=2)
plt.savefig("ismc_sigma_u_z.png")

u_z = dataPlot[5100:, 3]
n = len(u_z)
Y = scipy.fft(dataPlot[5100:, 3]) / n
k = arange(n)
T = n * h
frq = k / T
frq = frq[list(range(int(n / 2)))]
Y = Y[list(range(int(n / 2)))]
plt.plot(frq, abs(Y), "r")
plt.xlabel(r"freq (Hz)")
plt.title(r"Frequency spectrum of $\bar{u}^s$")
plt.savefig("ismc_u_freq.png")

# TODO
# compare with the reference
# ref = getMatrix(SiconosMatrix("result.ref"))
# if (norm(dataPlot - ref[1:,:]) > 1e-12):
#    print("Warning. The result is rather different from the reference file.")
