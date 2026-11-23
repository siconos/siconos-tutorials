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
# Example of delayed object introduction with time_of_birth parameter
#
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a cube as a convex shape
    io.add_convex_shape(
        "Cube",
        [
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (-1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, -1.0, 1.0),
        ],
    )

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (100, 100, 0.5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.3)

    # The cube objects are made with an unique Contactor : the cube shape.
    # As a mass is given, they are dynamic systems involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0

    # A first cube is introduced a the beginning of the simulation
    io.add_object(
        "cube0",
        [Contactor("Cube")],
        translation=[0, 0, 2],
        velocity=[10, 0, 0, 1, 1, 1],
        mass=1,
    )

    # the second cube introduction is delayed. It is crearted in the simulation
    # a time 0.5
    io.add_object(
        "cube1",
        [Contactor("Cube")],
        translation=[0, 0, 2],
        velocity=[10, 0, 0, 1, 1, 1],
        mass=1,
        time_of_birth=0.5,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, 0])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

test = True
if test:
    T = 2.0
else:
    T = 10.0

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.005


run_options["solver_options"] = options
run_options["multipoints_iterations"] = True
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 20
run_options["output_frequency"] = None

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
