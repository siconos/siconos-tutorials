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

# Here we demonstrate the use of two very small cubes simulated with
# geometry at 1000x its scale.  Without scaling, contact detection is
# incorrect for an object of this size.  When contact detection is
# performed with scaled-up geometries, the contact points are
# correctly generated and used by Siconos.
#
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# We need to pass some options to the Bullet backend

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a cube as a convex shape
    io.add_convex_shape(
        "CubeCH",
        [
            (-0.001, 0.001, -0.001),
            (-0.001, -0.001, -0.001),
            (-0.001, -0.001, 0.001),
            (-0.001, 0.001, 0.001),
            (0.001, 0.001, 0.001),
            (0.001, 0.001, -0.001),
            (0.001, -0.001, -0.001),
            (0.001, -0.001, 0.001),
        ],
    )

    # Definition of a cube as a primitive shape
    io.add_primitive_shape("CubeP", "Box", (0.002, 0.002, 0.002))

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (0.1, 0.1, 0.01))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object(
        "cubeCH",
        [Contactor("CubeCH")],
        translation=[0, 0.003, 0.005],
        velocity=[0.1, 0, 0, 1, 1, 1],
        mass=0.1,
    )

    # The primitive cube geometry object
    io.add_object(
        "cubeP",
        [Contactor("CubeP")],
        translation=[0, -0.003, 0.005],
        velocity=[0.1, 0, 0, 1, 1, 1],
        mass=0.1,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.005])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1000.0
bullet_options.contactBreakingThreshold = 0.001


options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 10.0
run_options["h"] = 0.005


run_options["solver_options"] = options
run_options["bullet_options"] = bullet_options
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 4
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
