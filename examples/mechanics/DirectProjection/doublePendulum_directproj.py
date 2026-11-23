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
#
#
# Example of a double pendulum
#
# Req : mechanics component with bullet.

from math import pi

import siconos.integrators
import siconos.numerics as sn
import siconos.simulation
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.mechanics.collision.tools import Contactor

# length of first branch
l1 = 30

# length of second branch
l2 = 10

# mass 1
m1 = 1

# mass 2
m2 = 10

# radius of first ball
r1 = 1

# radisu of second ball
r2 = 1

# gap between ball branch to avoid friction
gap = 0.0

# gap between ground
hgap = 0.1

# size of a brick
bx = 5
by = 2
bz = 2

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    io.add_primitive_shape("Arm1", "Cylinder", (0.3, l1))
    io.add_primitive_shape("Arm2", "Cylinder", (0.3, l2))
    io.add_primitive_shape("Mass1", "Sphere", (r1,))
    io.add_primitive_shape("Mass2", "Sphere", (r2,))

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (100, 100, 0.5))

    # the brick shape
    io.add_primitive_shape("Brick", "Box", (bx, by, bz))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1)

    # first branch + first mass the center of gravity is at the center of the
    # Mass1
    io.add_object(
        "arm1",
        [
            Contactor("Mass1"),
            Contactor("Arm1", relative_translation=[0, r1 + l1 / 2.0, 0]),
        ],
        translation=[0, 0, r2 + gap + r2 + l2 + r1 + hgap],
        orientation=((1, 0, 0), pi / 2),
        mass=m1,
    )

    # second branch + second mass
    io.add_object(
        "arm2",
        [
            Contactor("Mass2"),
            Contactor("Arm2", relative_translation=[0, r2 + l2 / 2.0, 0]),
        ],
        translation=[0, 0, r2 + gap],
        orientation=((1, 0, 0), pi / 2),
        velocity=[0, 20, 0, 0, 0, 0],
        mass=m2,
    )

    io.add_joint(
        "joint1",
        "arm1",
        "arm2",
        points=[[0, 0, -r1]],
        axes=[[1, 0, 0]],
        joint_class="PivotJointR",
        absolute=False,
    )

    io.add_joint(
        "joint2",
        "arm1",
        points=[[0, r2 + gap + r2 + l2 + r1 + hgap + l1, 0]],
        axes=[[1, 0, 0]],
        joint_class="PivotJointR",
        absolute=False,
    )

    # a brick wall
    H = 3  # heigh
    L = 2  # length
    for k in range(H - 1):
        for n in range(L):
            io.add_object(
                f"brick{k + n * H}",
                [Contactor("Brick")],
                translation=[
                    n * bx - L * bx / 2.0 + (k % 2) * bx / 2.0,
                    -5,
                    k * bz + bz / 2.0,
                ],
                mass=2,
            )

    k = H - 1
    for n in range(1, L):
        io.add_object(
            f"brick{k + n * H}",
            [Contactor("Brick")],
            translation=[
                n * bx - L * bx / 2.0 + (k % 2) * bx / 2.0,
                -5,
                k * bz + bz / 2.0,
            ],
            mass=2,
        )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.25])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8
options.update_internal(1, sn.solver_ids.SICONOS_FRICTION_3D_ONECONTACT_NSN)
# options=None


options_pos = sn.SolverOptions(sn.solver_ids.SICONOS_MLCP_ENUM)
options_pos.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options_pos.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12


test = True
h = 0.01
if test:
    T = 3 * h
else:
    T = 20

with MechanicsHdf5Runner(mode="r+") as io:
    io.run(
        h=h,
        T=T,
        solver_options=options,
        solver_options_pos=options_pos,
        time_stepping=siconos.simulation.TimeSteppingDirectProjection,
        osi=siconos.integrators.MoreauJeanDirectProjectionOSI,
        projection_itermax=3,
        projection_tolerance=1e-5,
        projection_tolerance_unilateral=1e-5,
    )
