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

import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Load a mesh.  The example mesh is a low-poly version of the
    # Stanford Bunny by Saf, license: Creative Commons - Attribution.
    # Taken from http://www.thingiverse.com/thing:466857
    io.add_mesh_from_file(
        "Bunny", "bunny.stl", scale=0.01, insideMargin=0.0, outsideMargin=0.0
    )

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (10, 10, 1.0))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0.2)

    # The mesh object made with an unique Contactor : the bunny shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0.  Ideally we would calculate or add
    # the mesh inertia matrix here, but it is not done in this
    # example.
    for i in range(3):
        io.add_object(
            f"bunny{i}",
            [Contactor("Bunny")],
            translation=[0, 0, 3 + i],
            orientation=[1, 0, 0, 0],
            velocity=[0, 0, 0, 1, 0, 0],
            mass=1,
        )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.5])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8
test = True
if test:
    T = 5.0
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-3
else:
    T = 20.0

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0.0
run_options["T"] = T
run_options["h"] = 0.005
run_options["theta"] = 0.50001
run_options["solver_options"] = options
run_options["multipoints_iterations"] = True

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
