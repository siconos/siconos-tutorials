# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2024 INRIA.
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
import numpy as np
import siconos.integrators
import siconos.modeling as sm
import siconos.nonsmooth_formulations
import siconos.plot_config as sicoplot
import siconos.simulation

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0  # start time
T = 10  # end time
h = 0.005  # time step
r = 0.1  # ball radius
g = 9.81  # gravity
m = 1  # ball mass
e = 0.9  # restitution coeficient
theta = 0.5  # theta scheme

#
# dynamical system
#
ndof = 3
initial_position = np.array([1, 0, 0], dtype=np.float64)
initial_velocity = np.array([0, 0, 0], dtype=np.float64)

mass = np.eye(ndof, dtype=np.float64, order="F")
mass[2, 2] = 2.0 / 5 * r * r


ball = sm.LagrangianLinearTIDS(initial_position, initial_velocity, mass, sm.alias_t)
# set external forces
weight_np = np.array([-m * g, 0, 0], dtype=np.float64)

ball.setConstantFext(weight_np, sm.alias_t)

#
# Interaction ball-floor
H = np.array([[1, 0, 0]], dtype=np.float64, order="F")

nslaw = sm.NewtonImpactNSL(e)
relation = sm.LagrangianLinearTIR(H)
inter = sm.Interaction(nslaw, relation)

# NSDS
bouncingBall = sm.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBall.insertDynamicalSystem(ball)

# link the interaction and the dynamical system
bouncingBall.link(inter, ball)

#
# Simulation
#

# (1) OneStepIntegrators
OSI = siconos.integrators.MoreauJeanOSI(theta)

# (2) Time discretisation --
t = siconos.simulation.TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = siconos.nonsmooth_formulations.LCP()

# (4) Simulation setup with (1) (2) (3)
s = siconos.simulation.TimeStepping(bouncingBall, t, OSI, osnspb)

# the number of time steps
N = int((T - t0) / h)

# Get the values to be plotted
# ->saved in a matrix dataPlot

dataPlot = np.zeros((N + 1, 5))

#
# numpy pointers on dense Siconos vectors
#
q = ball.q()
v = ball.velocity()
p = ball.p(1)
lambda_ = inter.lambda_python(1)


#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[0]
dataPlot[0, 2] = v[0]
dataPlot[0, 3] = p[0]
dataPlot[0, 4] = lambda_[0]

k = 1

# time loop
while s.hasNextEvent():
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0]
    dataPlot[k, 4] = lambda_[0]

    k += 1
    s.nextStep()

#
# comparison with the reference file
#
ref = np.loadtxt("BouncingBallTS.ref", skiprows=1)
assert np.allclose(dataPlot, ref, atol=1e-12)

#
# plots
#
plt.subplot(411)
plt.title("position")
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title("velocity")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title("reaction")
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title("lambda")
plt.grid()

if enable_plot:
    plt.show()
else:
    plt.savefig("bbts.png")
