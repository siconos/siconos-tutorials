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

import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "toolbox")
)
import numpy as np
import result_checks as rc
import siconos.numerics as sn
import siconos.plot_config as sicoplot
from numpy import zeros
from numpy.linalg import norm
from siconos.integrators import MoreauJeanOSI
from siconos.mechanics.collision import (
    RigidBodyDS,
    SiconosBox,
    SiconosContactor,
    SiconosContactorSet,
    SiconosPlane,
)
from siconos.mechanics.collision.bullet import SiconosBulletCollisionManager
from siconos.modeling import (
    NewtonImpactFrictionNSL,
    NonSmoothDynamicalSystem,
    alias_t,
    interactions,
)
from siconos.nonsmooth_formulations import FrictionContact
from siconos.simulation import TimeDiscretisation, TimeStepping

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0  # start time
T = 20  # end time
h = 0.005  # time step

g = 9.81  # gravity

theta = 0.5  # theta scheme

#
# dynamical system
#
position_init = 10
velocity_init = 0

initial_position = np.array([0, 0, position_init, 1, 0, 0, 0], dtype=np.float64)
initial_velocity = np.array([0, 0, velocity_init, 0, 0, 0], dtype=np.float64)
inertia = np.eye(3, dtype=np.float64, order="F")

# a box shape
box1 = SiconosBox(1.0, 1.0, 1.0)

# A Bullet Dynamical System : a shape + a mass (1.0) + position and velocity
body = RigidBodyDS(initial_position, initial_velocity, 1.0, inertia)

# Add the shape, wrapped in a SiconosContactor, to the body's
# contactor set.
contactorOffset = np.array([0, 0, 0, 1, 0, 0, 0], dtype=np.float64)
contactor = SiconosContactor(box1, contactorOffset)
body.contactors().append(contactor)

# set external forces
weight = np.array([0, 0, -body.scalarMass * g], dtype=np.float64)
body.setConstantFext(weight, alias_t)

#
# Model
#
bouncingBox = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBox.insertDynamicalSystem(body)

#
# Simulation
#

# (1) OneStepIntegrators
osi = MoreauJeanOSI(theta)

ground = SiconosPlane()
groundOffset = np.array([0, 0, -0.5, 1, 0, 0, 0], dtype=np.float64)

# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = FrictionContact(3)

osnspb.numericsSolverOptions().iparam[0] = 10000  # as in the C++ version
osnspb.numericsSolverOptions().dparam[0] = 1e-5
osnspb.setMaxSize(16384)
osnspb.setMStorageType(sn.params.NM_SPARSE_BLOCK)
osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)


# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0.8, 0.0, 0.0, 3)

# (5) broadphase contact detection
broadphase = SiconosBulletCollisionManager()

# insert a non smooth law for contactors id 0
broadphase.insertNonSmoothLaw(nslaw, 0, 0)

# The ground is a static object
# we give it a group contactor id : 0
scs = SiconosContactorSet()
scs.append(SiconosContactor(ground, contactorOffset))
broadphase.addStaticBody(scs, groundOffset)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = TimeStepping(bouncingBox, timedisc)
simulation.insertInteractionManager(broadphase)

simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)


# Get the values to be plotted
# ->saved in a matrix dataPlot

N = int((T - t0) / h)
dataPlot = zeros((N + 1, 4))

#
# numpy pointers on dense Siconos vectors
#
q = body.q()
v = body.twist()

#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[2]
dataPlot[0, 2] = v[2]

k = 1

# data for the checks of the results (see the end of the script)
E0 = rc.energy(q, v, body.totalInertiaMatrix, g)
energy_increase = rc.MaxIncrease()
energy_increase.push(E0)
min_lowest_z = rc.box_lowest_z(q, 0.5, 0.5, 0.5)
first_impact_time, v_before, v_after = -1.0, 0.0, 0.0

# time loop
while simulation.hasNextEvent():
    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    energy_increase.push(rc.energy(q, v, body.totalInertiaMatrix, g))
    min_lowest_z = min(min_lowest_z, rc.box_lowest_z(q, 0.5, 0.5, 0.5))
    if first_impact_time < 0 and dataPlot[k - 1, 2] < 0 and v[2] > 0:
        first_impact_time, v_before, v_after = dataPlot[k, 0], dataPlot[k - 1, 2], v[2]

    # if (broadphase.collisionWorld().getDispatcher().getNumManifolds() > 0):
    if (
        broadphase.statistics().new_interactions_created
        + broadphase.statistics().existing_interactions_processed
    ) > 0 and bouncingBox.topology().numberOfIndexSet() == 2:
        index1 = interactions(simulation.indexSet(1))
        if len(index1) == 4:
            dataPlot[k, 3] = (
                norm(index1[0].lambda_python(1))
                + norm(index1[1].lambda_python(1))
                + norm(index1[2].lambda_python(1))
                + norm(index1[3].lambda_python(1))
            )

    k += 1
    simulation.nextStep()

#
# checks of the results (as in BulletBouncingBox.cpp): the trajectory depends on the
# contact detection by Bullet, it is not compared with a reference file.
#
e, ground_z = 0.8, -0.5
v_impact = np.sqrt(2 * g * position_init)
chk = rc.Checks("BulletBouncingBox.py")
chk.near(
    "time of the first impact", first_impact_time, np.sqrt(2 * position_init / g), 2 * h
)
chk.near("restitution of the first impact", -v_after / v_before, e, 0.01)
chk.at_least("lowest point of the box", min_lowest_z, ground_z - (h * v_impact + 0.05))
chk.at_most("max energy increase in one step", energy_increase.value, 1e-5 * E0)
chk.at_most(
    "final energy (ground at -0.5 m g)",
    rc.energy(q, v, body.totalInertiaMatrix, g),
    1e-4 * E0,
)
chk.near("final lowest point", rc.box_lowest_z(q, 0.5, 0.5, 0.5), ground_z, 0.05)
chk.at_most("final vertical velocity", abs(v[2]), 1e-2)
chk.at_least("final orientation: lying on a face", rc.vertical_axis_cosine(q), 0.999)


#
# plots
#

if enable_plot:
    plt.subplot(511)
    plt.title("position")
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 1])
    plt.grid()
    plt.subplot(513)
    plt.title("velocity")
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 2])
    plt.grid()
    plt.subplot(515)
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 3])
    plt.title("lambda")
    plt.grid()
    plt.savefig("result.png")
    if enable_plot:
        plt.show()

np.savetxt("BouncingBox-py.dat", dataPlot)
sys.exit(chk.report())
