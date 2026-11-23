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
# -----------------------------------------------------------------------
#
#  CircuitRLCD  : sample of an electrical circuit involving :
#  - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
#  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
#    with the oscillator
#
#  Expected behavior :
#  The initial state of the oscillator provides an initial energy.
#  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
#  A positive voltage across the capacitor allows current to flow
#  through the resistor-diode branch , resulting in an energy loss :
#  the oscillation damps.
#
#  State variables :
#  - the voltage across the capacitor (or inductor)
#  - the current through the inductor
#
#  Since there is only one dynamical system, the interaction is defined by :
#  - a complementarity law between diode current and voltage where y stands
#    for the reverse voltage across the diode and lambda stands for the
#    the diode current
#  - a linear time invariant relation between the state variables and
#    y and lambda (derived from Kirchhoff laws)
#
# -----------------------------------------------------------------------


import numpy as np
import siconos.integrators as si
import siconos.modeling as sm
import siconos.nonsmooth_formulations as snsf
import siconos.plot_config as sicoplot
import siconos.simulation as ss

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0.0
T = 5.0e-3  # Total simulation time
h_step = 10.0e-6  # Time step
Lvalue = 1e-2  # inductance
Cvalue = 1e-6  # capacitance
Rvalue = 1e3  # resistance
Vinit = 10.0  # initial voltage


#
# dynamical system
#
init_state = np.array([Vinit, 0], dtype=np.float64, order="F")

A = np.zeros((2, 2), dtype=np.float64, order="F")
A.flat[...] = [0.0, -1.0 / Cvalue, 1.0 / Lvalue, 0.0]

LSCircuitRLCD = sm.FirstOrderLinearDS(init_state, sm.alias_t)
LSCircuitRLCD.setConstantA(A, sm.alias_t)

#
# Interactions
#

C = np.array([[-1.0, 0.0]], dtype=np.float64, order="F")

D = np.array([[Rvalue]], dtype=np.float64, order="F")

B = np.array([[-1.0 / Cvalue], [0.0]], dtype=np.float64, order="F")

LTIRCircuitRLCD = sm.FirstOrderLinearTIR(C, B)
LTIRCircuitRLCD.setConstantD(D)

nslaw = sm.ComplementarityConditionNSL(1)
InterCircuitRLCD = sm.Interaction(nslaw, LTIRCircuitRLCD)


#
# Model
#
CircuitRLCD = sm.NonSmoothDynamicalSystem(t0, T)
CircuitRLCD.setTitle("CircuitRLCD")

#   add the dynamical system in the non smooth dynamical system
CircuitRLCD.insertDynamicalSystem(LSCircuitRLCD)

#   link the interaction and the dynamical system
CircuitRLCD.link(InterCircuitRLCD, LSCircuitRLCD)


#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5000000000001
aOSI = si.EulerMoreauOSI(theta)

# (2) Time discretisation
aTiDisc = ss.TimeDiscretisation(t0, h_step)

# (3) Non smooth problem
aLCP = snsf.LCP()

# (4) Simulation setup with (1) (2) (3)
aTS = ss.TimeStepping(CircuitRLCD, aTiDisc, aOSI, aLCP)

# end of model definition

#
# computation
#


h = aTS.timeStep()
print("Timestep : ", h)
# Number of time steps
N = int((T - t0) / h) + 1
print("Number of steps : ", N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

dataPlot = np.zeros([N + 1, 6], dtype=np.float64)

x = LSCircuitRLCD.x()
y = InterCircuitRLCD.y(0)
lambda_ = InterCircuitRLCD.lambda_python(0)
InterCircuitRLCD.computeInput(t0, 0)
InterCircuitRLCD.computeOutput(t0, 0)
# For the initial time step:
# time
k = 0
#  inductor voltage
dataPlot[k, 1] = x[0]

# inductor current
dataPlot[k, 2] = x[1]

# diode voltage
dataPlot[k, 3] = -y[0]

# diode current
dataPlot[k, 4] = lambda_[0]
dataPlot[k, 5] = LSCircuitRLCD.r()[0]

k += 1
while aTS.hasNextEvent():
    aTS.computeOneStep()
    # aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    # inductor current
    dataPlot[k, 2] = x[1]
    # diode  voltage
    dataPlot[k, 3] = -y[0]
    # diode  current
    dataPlot[k, 4] = lambda_[0]
    dataPlot[k, 5] = 0.0
    k += 1
    aTS.nextStep()

# comparison with reference file

ref = np.loadtxt("CircuitRLCD.ref", skiprows=1)

assert np.allclose(ref, dataPlot, atol=1e-8)

if enable_plot:
    #
    # plots
    #
    plt.subplot(411)
    plt.title("inductor voltage")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 1])
    plt.grid()
    plt.subplot(412)
    plt.title("inductor current")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 2])
    # plt.plot(dataPlot[0:k - 1, 0], ref[0:k - 1, 2])
    plt.grid()
    plt.subplot(413)
    plt.title("diode  voltage")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 3])
    plt.subplot(414)
    plt.title("diode current")
    plt.plot(dataPlot[0 : k - 1, 0], dataPlot[0 : k - 1, 4])
    plt.savefig("circuit_rlcd.png")
