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
import os

import numpy as np
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# ctest sets SICONOS_EXAMPLES_SHORT_RUN: shorter simulation for the tests,
# full simulation otherwise.
short_run = "SICONOS_EXAMPLES_SHORT_RUN" in os.environ


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a capsule
    R = 0.1
    L = 2.0
    io.add_primitive_shape("Cap", "Capsule", (R, L))

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (10, 10, 1.0))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.9)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    mass_test = 1.0
    inertia_test = np.eye(3)

    inertia_test[0, 0] = 0.25 * mass_test * R * R + 1 / 3.0 * mass_test * L * L
    inertia_test[1, 1] = 0.5 * mass_test * R * R
    inertia_test[2, 2] = 0.25 * mass_test * R * R + 1 / 3.0 * mass_test * L * L
    print(inertia_test)
    orientation_test = [0.14, 0.7, 0.7, 0]
    io.add_object(
        "cap",
        [Contactor("Cap")],
        translation=[0, 0, 4],
        orientation=orientation_test,
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
        inertia=inertia_test,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.1])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-6

test = short_run
if test:
    T = 0.5
else:
    T = 20.0

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.0005


run_options["solver_options"] = options
# run_options["bullet_options"]=bullet_options
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
    io.run(run_options)
