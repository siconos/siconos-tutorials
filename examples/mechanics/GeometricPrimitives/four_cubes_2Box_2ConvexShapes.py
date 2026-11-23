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
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

edge_length = 0.1
plane_length = 2.0

velocity_init = 0.0
angular_velocity_init = 0.0

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # # Definition of a cube as a convex shape
    io.add_convex_shape(
        "CubeCS1",
        [
            (-edge_length, edge_length, -edge_length),
            (-edge_length, -edge_length, -edge_length),
            (-edge_length, -edge_length, edge_length),
            (-edge_length, edge_length, edge_length),
            (edge_length, edge_length, edge_length),
            (edge_length, edge_length, -edge_length),
            (edge_length, -edge_length, -edge_length),
            (edge_length, -edge_length, edge_length),
        ],
    )

    io.add_convex_shape(
        "CubeCS2",
        [
            (-edge_length, edge_length, -edge_length),
            (-edge_length, -edge_length, -edge_length),
            (-edge_length, -edge_length, edge_length),
            (-edge_length, edge_length, edge_length),
            (edge_length, edge_length, edge_length),
            (edge_length, edge_length, -edge_length),
            (edge_length, -edge_length, -edge_length),
            (edge_length, -edge_length, edge_length),
        ],
    )

    # Alternative to the previous convex shape definition.
    io.add_primitive_shape(
        "CubePrim1", "Box", (2 * edge_length, 2 * edge_length, 2 * edge_length)
    )

    io.add_primitive_shape(
        "CubePrim2", "Box", (2 * edge_length, 2 * edge_length, 2 * edge_length)
    )

    # Definition of the ground shape
    io.add_primitive_shape(
        "Ground", "Box", (plane_length, plane_length, plane_length / 10.0)
    )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0.5)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object(
        "cube1",
        [Contactor("CubeCS1")],
        translation=[0, 0, 2],
        velocity=[
            velocity_init,
            0,
            0,
            angular_velocity_init,
            angular_velocity_init,
            angular_velocity_init,
        ],
        mass=1,
    )

    io.add_object(
        "cube2",
        [Contactor("CubeCS2")],
        translation=[0, 0, 2 + 3 * edge_length],
        velocity=[
            velocity_init,
            0,
            0,
            angular_velocity_init,
            angular_velocity_init,
            angular_velocity_init,
        ],
        mass=1,
    )

    io.add_object(
        "cubeP1",
        [Contactor("CubePrim1")],
        translation=[0, 3 * edge_length, 2],
        velocity=[
            velocity_init,
            0,
            0,
            angular_velocity_init,
            angular_velocity_init,
            angular_velocity_init,
        ],
        mass=1,
    )

    io.add_object(
        "cubeP2",
        [Contactor("CubePrim2")],
        translation=[0, 3 * edge_length, 2 + 3 * edge_length],
        velocity=[
            velocity_init,
            0,
            0,
            angular_velocity_init,
            angular_velocity_init,
            angular_velocity_init,
        ],
        mass=1,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, 0])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

step = 1000
hstep = 0.005

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = step * hstep
run_options["h"] = hstep


run_options["solver_options"] = options
run_options["multipoints_iterations"] = True
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 1
run_options["output_frequency"] = 10

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

with MechanicsHdf5Runner(mode="r+", collision_margin=0.05) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    io.run(run_options)
