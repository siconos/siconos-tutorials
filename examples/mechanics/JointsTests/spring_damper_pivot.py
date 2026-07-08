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
import numpy as np

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.modeling as sm

# An example of applying force to the axis of a joint, and applying
# spring and virtual damping by measuring position and velocity along
# the same axis.

# Note: This example is to demonstrate external measurement of joint
# positions and application of forces to dynamical systems attached to
# joints.  In practice it is better to use internal forces (fInt,
# mInt) to model joint spring-dampers, see folder
# JointsTestsWithInternalForces, and extra Relations with associated
# Non-Smooth Laws to model non-linearities such as joint stops and
# friction, see JointsTestsWithContactDetection.

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of two bars connected by a prismatic joint
    io.add_primitive_shape("Bar", "Box", (1, 0.1, 0.1))
    io.add_object(
        "bar1",
        [Contactor("Bar")],
        [0, 0, 2],
        orientation=((0, 0, 1), np.pi / 2),
        mass=1.0,
        velocity=[0, 0, 0, 0, 0, 1],
    )
    io.add_object(
        "bar2",
        [
            Contactor("Bar", relative_translation=np.array([0.0, 0.1, 0.0])),
            Contactor("Bar", relative_translation=np.array([0.0, -0.1, 0.0])),
        ],
        [0, 0, 2],
        orientation=((0, 0, 1), np.pi / 2),
        mass=1.0,
    )
    io.add_joint(
        "joint1",
        "bar1",
        "bar2",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )

    # Definition of the ground
    io.add_primitive_shape("Ground", "Box", (5, 5, 0.1))
    io.add_object("ground", [Contactor("Ground")], [0, 0, -0.05])
    io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0.0)


class Ctrl(object):
    def initialize(self, io):
        self.count = 0
        self.topo = io._nsds.topology()
        self.ds1 = self.topo.getDynamicalSystem("bar1")
        self.ds2 = self.topo.getDynamicalSystem("bar2")
        self.joint1 = self.topo.getInteraction("joint1").relation()
        self.ds1.setIsMextExpressedInInertialFrame(True)
        self.ds2.setIsMextExpressedInInertialFrame(True)

        # Apply initial forces
        self.step()

    def step(self):
        self.count += 1
        torque1 = np.zeros(3)
        torque2 = np.zeros(3)

        # Get the position and use it to project a torque vector
        # onto the DoF (spring torque)
        angle = np.array([0], dtype=np.float64)
        self.joint1.computehDoF(self.ds1.q(), self.ds2.q(), angle, 0)
        setpoint = np.pi / 4
        ang_diff = setpoint - angle[0]
        spring_torque = (
            np.array(self.joint1.normalDoF(self.ds1.q(), axis=0)) * ang_diff * 500.0
        )
        # Get the velocity of each body projected onto the DoF and
        # calculate their difference (damping torque)
        vel1 = self.joint1.projectVectorDoF(
            self.ds1.angularVelocityInBodyFrame(), self.ds1.q(), axis=0
        )
        vel2 = self.joint1.projectVectorDoF(
            self.ds2.angularVelocityInBodyFrame(), self.ds1.q(), axis=0
        )
        vel_diff = vel1 - vel2
        damping_torque = vel_diff * 5.0
        # Calculate total torques for each body
        torque1 += -(spring_torque + damping_torque) / 2
        torque2 += +(spring_torque + damping_torque) / 2
        self.ds1.setConstantMext(torque1, sm.copy_t)
        self.ds2.setConstantMext(torque2, sm.copy_t)


options = sn.SolverOptions(sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

# Load and run the simulation
with MechanicsHdf5Runner(mode="r+") as io:
    io.run(
        t0=0,
        T=2.0,
        h=0.001,
        theta=0.5,
        Newton_max_iter=1,
        controller=Ctrl(),
        solver_options=options,
        output_frequency=1,
    )
