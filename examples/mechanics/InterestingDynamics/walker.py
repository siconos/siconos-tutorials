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
# A simple walking mechanism (in progress, does not work!)
#
import math

import numpy as np
import siconos.modeling as sm
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

pi = math.pi

motor_id = None

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Only touch the ground, ignore contacts between links of the robot
    io.add_Newton_impact_friction_nsl(
        "contact", collision_group1=0, collision_group2=1, mu=0.3, e=0.0
    )

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (100, 100, 2))
    io.add_object(
        "ground", [Contactor("Ground", collision_group=1)], translation=[0, 0, -1]
    )

    # Define shape of a bar-shaped link
    io.add_primitive_shape("Bar1", "Box", (10, 1, 1))
    io.add_primitive_shape("Bar2", "Box", (12.67, 1, 1))

    # Core
    bar1 = io.add_object(
        "bar1",
        [
            Contactor("Bar1"),
            Contactor(
                "Bar1",
                relative_translation=[5.5, 0, 0],
                relative_orientation=[(0, 0, 1), pi / 2],
            ),
            Contactor(
                "Bar1",
                relative_translation=[-5.5, 0, 0],
                relative_orientation=[(0, 0, 1), pi / 2],
            ),
        ],
        translation=[0, 0, 11],
        mass=1,
    )
    bar1_id = bar1.attrs["id"]

    # Legs
    io.add_object(
        "bar2",
        [Contactor("Bar1")],
        translation=[5.5, 5.5, 6],
        orientation=[(0, 1, 0), -pi / 2],
        mass=1,
    )

    io.add_object(
        "bar3",
        [Contactor("Bar1")],
        translation=[5.5, -5.5, 6],
        orientation=[(0, 1, 0), -pi / 2],
        mass=1,
    )

    io.add_object(
        "bar4",
        [Contactor("Bar1")],
        translation=[-5.5, 5.5, 6],
        orientation=[(0, 1, 0), -pi / 2],
        mass=1,
    )

    io.add_object(
        "bar5",
        [Contactor("Bar1")],
        translation=[-5.5, -5.5, 6],
        orientation=[(0, 1, 0), -pi / 2],
        mass=1,
    )

    io.add_joint(
        "joint1",
        "bar2",
        "bar1",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint2",
        "bar3",
        "bar1",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint3",
        "bar4",
        "bar1",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint4",
        "bar5",
        "bar1",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )

    # Stabilizing leg links to body

    # See mech.spbstu.ru/Dzenushko_Dainis:_Walking_mechanisms_survey
    # for a survey of better walking mechanisms

    io.add_object(
        "bar6",
        [Contactor("Bar2")],
        translation=[0, 6.5, 7],
        orientation=[(0, 1, 0), -pi / 8],
        mass=1,
    )
    io.add_object(
        "bar7",
        [Contactor("Bar2")],
        translation=[0, -6.5, 7],
        orientation=[(0, 1, 0), -pi / 8],
        mass=1,
    )
    io.add_object(
        "bar8",
        [Contactor("Bar2")],
        translation=[0, 6.5, 7],
        orientation=[(0, 1, 0), -pi / 8],
        mass=1,
    )
    io.add_object(
        "bar9",
        [Contactor("Bar2")],
        translation=[0, -6.5, 7],
        orientation=[(0, 1, 0), -pi / 8],
        mass=1,
    )
    io.add_joint(
        "joint5",
        "bar6",
        "bar2",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint6",
        "bar7",
        "bar3",
        [[4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint7",
        "bar8",
        "bar4",
        [[-4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint8",
        "bar9",
        "bar5",
        [[-4.5, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )

    io.add_joint(
        "joint9", "bar6", "bar8", None, [[1, 0, 0]], "PrismaticJointR", absolute=False
    )
    io.add_joint(
        "joint10", "bar7", "bar9", None, [[1, 0, 0]], "PrismaticJointR", absolute=False
    )

    io.add_joint(
        "joint11",
        "bar2",
        "bar3",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )
    io.add_joint(
        "joint12",
        "bar4",
        "bar5",
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
    )

    # # Harmonic oscillator on Y-axis angular velocity = a+b*cos(omega*time+phi))
    freq = 2
    amp = 0.2
    io.add_boundary_condition(
        "vibration",
        "bar1",
        indices=[4],  # Y-angular axis is index 4 into ds->v
        bc_class="HarmonicBC",
        a=[0.0],
        b=[amp],
        omega=[pi * freq],
        phi=[0.0],
    )

forces_on_body = {}


def my_forces(body):
    g = 9.81
    weight = np.array([0, 0, -body.scalarMass * g])
    # twist = np.array([0, 0, 0, 0, 300, 0])
    push = np.zeros(3, dtype=np.float64)
    force = weight.copy()
    if body.number() == bar1_id:
        force = weight + push  # + twist
    body.setConstantFext(force, sm.copy_t)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
test = True

if test:
    T = 0.5
else:
    T = 30.0


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 5e-4

run_options["multipoints_iterations"] = True
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = None
run_options["set_external_forces"] = my_forces


with MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
