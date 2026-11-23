#!/usr/bin/env python
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
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
# Example of a small driving 4-wheeled car with no steering
#

import numpy as np
import siconos.modeling as sm
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

bumpy_terrain = True  # false = flat ground
use_torque = True  # false = boundary condition

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.add_Newton_impact_friction_nsl("contact", mu=0.6, e=0.0)

    # Definition of a wheel (Radius = 1.5, width = 0.5)
    io.add_primitive_shape("Wheel", "Cylinder", (1.5, 0.5))

    # Definition of the car body
    io.add_primitive_shape("Body", "Box", (6, 5.4, 0.5))

    if not bumpy_terrain:
        # Definition of the ground shape
        io.add_primitive_shape("Ground", "Box", (50, 50, 1))

        # The ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0, 0])

    else:
        # The heightmap is a 100x200 grid with a height of 1.
        heightmap = (
            np.random.uniform(0, 0.3, size=(100, 200))
            + [np.linspace(0, 2, 200)] * 100
            + [np.sin(np.linspace(0, np.pi * 12, 200))] * 100
        ).T

        # It has tangential extents 100 by 50.
        io.add_height_map("Terrain", heightmap, (100, 50))

        # It is also a static object (no mass)
        io.add_object("ground", [Contactor("Terrain")], translation=[40, 0, 0])

    h = 3
    io.add_object("wheel1", [Contactor("Wheel")], translation=[-3, 3, h], mass=1)
    io.add_object("wheel2", [Contactor("Wheel")], translation=[3, 3, h], mass=1)
    io.add_object("wheel3", [Contactor("Wheel")], translation=[-3, -3, h], mass=1)
    io.add_object("wheel4", [Contactor("Wheel")], translation=[3, -3, h], mass=1)
    io.add_object("body", [Contactor("Body")], translation=[0, 0, h], mass=1)

    io.add_joint(
        "joint1",
        "wheel1",
        "body",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint2",
        "wheel2",
        "body",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint3",
        "wheel3",
        "body",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint4",
        "wheel4",
        "body",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )

    # Constant value on Y-axis angular velocity (index 4) to make the wheels spin
    if not use_torque:
        ang_vel = 8.0
        io.add_boundary_condition(
            "spin1", "wheel1", indices=[4], bc_class="BoundaryCondition", v=[ang_vel]
        )
        io.add_boundary_condition(
            "spin2", "wheel2", indices=[4], bc_class="BoundaryCondition", v=[ang_vel]
        )
        io.add_boundary_condition(
            "spin3", "wheel3", indices=[4], bc_class="BoundaryCondition", v=[ang_vel]
        )
        io.add_boundary_condition(
            "spin4", "wheel4", indices=[4], bc_class="BoundaryCondition", v=[ang_vel]
        )


# A controller to give constant torque to the wheels
class roll:
    def initialize(self, io):
        ang_force = 6.0
        self.io = io
        topo = io._nsds.topology()
        self.wheels = [(topo.getDynamicalSystem(f"wheel{i}")) for i in [1, 2, 3, 4]]
        print(dir(self.wheels[0]))
        print(self.wheels[0].fext())
        self.wheel_force = np.array(self.wheels[0].fext())
        self.wheel_torque = np.array([0, ang_force, 0.0])

        # Same force and torque to all wheels
        [w.setConstantFext(self.wheel_force, sm.copy_t) for w in self.wheels]
        [w.setConstantMext(self.wheel_torque, sm.copy_t) for w in self.wheels]

    def step(self):
        pass


controller = None
if use_torque:
    controller = roll()

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0.0
run_options["T"] = 5
run_options["h"] = 0.005
run_options["controller"] = controller
run_options["theta"] = 0.50001
run_options["Newton_max_iter"] = 1


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
