#!/usr/bin/env python
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

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
import numpy as np
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of the ground shape
    box_width = 5.0
    box_height = 1.0
    io.add_primitive_shape("Ground", "Box", (box_width, box_width, box_height))

    # Definition of a capsule
    R = 0.1
    L = 2.0
    io.add_primitive_shape("Cap", "Capsule", (R, L))

    wall_height = 5.0

    io.add_primitive_shape("Wall_1", "Box", (wall_height, box_width, 0.2))

    io.add_primitive_shape("Wall_2", "Box", (box_width, wall_height, 0.2))

    io.add_primitive_shape("sphere_shape", "Sphere", (R,))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.0)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    mass_test = 1.0
    inertia_test = np.eye(3)

    inertia_test[0, 0] = 0.25 * mass_test * R * R + 1 / 3.0 * mass_test * L * L
    inertia_test[1, 1] = 0.5 * mass_test * R * R
    inertia_test[2, 2] = 0.25 * mass_test * R * R + 1 / 3.0 * mass_test * L * L

    io.add_object(
        "cap_1",
        [Contactor("Cap")],
        translation=[-box_width / 4.0, 0, 0.2 + box_height / 2.0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
        inertia=inertia_test,
    )

    io.add_object(
        "cap_2",
        [Contactor("Cap")],
        translation=[box_width / 4.0, 0, 0.2 + box_height / 2.0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
        inertia=inertia_test,
    )
    import math

    frequency = 1.0
    amplitude = 3.0
    io.add_boundary_condition(
        "vibration",
        "cap_1",
        indices=[0, 1, 2, 3, 4, 5],  # rotation around the  z-axis
        bc_class="HarmonicBC",
        a=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        b=[0.0, 0.0, 0.0, 0.0, 0.0, amplitude * frequency * 2.0 * math.pi],
        omega=[0.0, 0.0, 0.0, 0.0, 0.0, frequency * 2.0 * math.pi],
        phi=[0.0, 0.0, 0.0, 0.0, 0.0, math.pi / 2.0],
    )

    io.add_boundary_condition(
        "fixed_rotation",
        "cap_2",
        indices=[0, 1, 2, 3, 4, 5],  # rotation around the  z-axis
        bc_class="BoundaryCondition",
        v=[0.0, 0.0, 0.0, 0.0, 0.0, amplitude],
    )

    n_row = 10
    n_col = 10
    n_layer = 6

    spacing = box_width / (n_row + 1)

    for i in range(n_col):
        for j in range(n_row):
            for k in range(n_layer):
                x_loc = (i + 1) * spacing - box_width / 2.0
                y_loc = (j + 1) * spacing - box_width / 2.0
                z_loc = (k + 1) * (2 * R + 0.01)
                io.add_object(
                    "sphere_" + str(i) + "_" + str(j) + "_" + str(k),
                    [Contactor("sphere_shape")],
                    translation=[x_loc, y_loc, z_loc + 0.5 + box_height / 2.0],
                    velocity=[0, 0, 0, 0, 0, 0],
                    mass=1,
                )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.0])
    io.add_object(
        "wall_1",
        [Contactor("Wall_1")],
        translation=[-box_width / 2.0, 0, box_height],
        orientation=([0, 1.0, 0.0], math.pi / 2.0),
    )
    io.add_object(
        "wall_2",
        [Contactor("Wall_1")],
        translation=[box_width / 2.0, 0, box_height],
        orientation=([0, 1.0, 0.0], math.pi / 2.0),
    )
    io.add_object(
        "wall_3",
        [Contactor("Wall_2")],
        translation=[0, box_width / 2.0, box_height],
        orientation=([1, 0, 0.0], math.pi / 2.0),
    )
    io.add_object(
        "wall_4",
        [Contactor("Wall_2")],
        translation=[0, -box_width / 2.0, box_height],
        orientation=([1, 0, 0.0], math.pi / 2.0),
    )

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-3
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 100

test = True
if test:
    T = 0.2
    # T =0.01
else:
    T = 20.0

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.001


run_options["solver_options"] = options


run_options["Newton_max_iter"] = 3
run_options["output_frequency"] = 10

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0


with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
