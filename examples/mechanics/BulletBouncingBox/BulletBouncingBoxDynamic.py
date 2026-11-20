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
    copy_t,
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
mass = 1.0
weight = np.array([0, 0, -mass * g], dtype=np.float64)


def makeBox(pos=initial_position, vel=initial_velocity):
    box = SiconosBox(1.0, 1.0, 1.0)

    # A Bullet Dynamical System : a shape + a mass (1.0) + position and velocity
    body = RigidBodyDS(pos, vel, mass, inertia)

    # set external forces
    body.setConstantFext(weight, copy_t)

    # Add the shape, wrapped in a SiconosContactor, to the body's
    # contactor set.
    body.contactors().append(SiconosContactor(box))

    return body


# Initial box body
body = makeBox()

# set external forces
weight = np.array([0, 0, -body.scalarMass * g], dtype=np.float64)
body.setConstantFext(weight, copy_t)

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
scs.append(SiconosContactor(ground))
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

# time loop
new_position = np.array([0, 0, 3, 1, 0, 0, 0], dtype=np.float64)
new_velocity = np.array([0, 0, 0, 0, 0, 0], dtype=np.float64)

# data for the checks of the results (see the end of the script)
body2 = None  # the box added at step 100


def total_energy():
    E = rc.energy(q, v, body.totalInertiaMatrix, g)
    if body2 is not None:
        E += rc.energy(body2.q(), body2.twist(), body2.totalInertiaMatrix, g)
    return E


E0 = total_energy()
energy_increase = rc.MaxIncrease()
energy_increase.push(E0)
min_lowest_z = rc.box_lowest_z(q, 0.5, 0.5, 0.5)
min_centre_distance = np.inf
add_time, impact2_time, impact1_time = -1.0, -1.0, -1.0
v1_old, v2_old = v[2], 0.0

while simulation.hasNextEvent():
    # Add a second box dynamically to the simulation
    if k == 100:
        ds = makeBox(new_position, new_velocity)
        bouncingBox.insertDynamicalSystem(ds)
        body2 = ds
        add_time = simulation.startingTime()
        energy_increase.restart(total_energy())

    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    energy_increase.push(total_energy())
    min_lowest_z = min(min_lowest_z, rc.box_lowest_z(q, 0.5, 0.5, 0.5))
    if body2 is not None:
        q2, v2 = body2.q(), body2.twist()[2]
        min_lowest_z = min(min_lowest_z, rc.box_lowest_z(q2, 0.5, 0.5, 0.5))
        min_centre_distance = min(min_centre_distance, rc.centre_distance(q, q2))
        if impact2_time < 0 and v2_old < 0 and v2 > 0:
            impact2_time = dataPlot[k, 0]
        if impact1_time < 0 and v1_old < 0 and v[2] > v1_old + 1:
            impact1_time = dataPlot[k, 0]
        v2_old = v2
    v1_old = v[2]

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
# checks of the results (as in BulletBouncingBoxDynamic.cpp): the trajectory depends on
# the contact detection by Bullet, it is not compared with a reference file.
#
e, ground_z = 0.8, -0.5
z2 = new_position[2]
chk = rc.Checks("BulletBouncingBoxDynamic.py")
t2 = add_time + np.sqrt(2 * z2 / g)
chk.near("time of the first impact of box 2", impact2_time, t2, 2 * h)
# first impact of box 1, on box 2 bouncing up (ideal 1D impacts): z1 - z2 = 1
tt = np.arange(t2, T, 1e-6)
gap = (
    (position_init - 0.5 * g * tt**2)
    - (e * g * (t2 - add_time) * (tt - t2) - 0.5 * g * (tt - t2) ** 2)
    - 1
)
chk.near(
    "time of the first impact of box 1 (on box 2)",
    impact1_time,
    tt[np.argmax(gap <= 0)],
    3 * h,
)
margin = h * np.sqrt(2 * g * position_init) + 0.05
chk.at_least("lowest point of the boxes", min_lowest_z, ground_z - margin)
chk.at_least("distance between the centres", min_centre_distance, 1 - margin)
chk.at_most("max energy increase in one step", energy_increase.value, 1e-5 * E0)
for name, b in (("box 1", body), ("box 2", body2)):
    chk.at_most(f"final vertical velocity of {name}", abs(b.twist()[2]), 1e-2)
    chk.at_least(
        f"final orientation of {name}: on a face", rc.vertical_axis_cosine(b.q()), 0.999
    )


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
    plt.savefig("result_dynamic.png")
    if enable_plot:
        plt.show()

np.savetxt("BouncingBoxDynamic-py.dat", dataPlot)
sys.exit(chk.report())
