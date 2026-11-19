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
import os

import siconos.numerics as sn
from siconos.integrators import MoreauJeanDirectProjectionOSI
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.mechanics.collision.tools import Contactor
from siconos.simulation import TimeSteppingDirectProjection

# ctest sets SICONOS_EXAMPLES_SHORT_RUN: shorter simulation for the tests,
# full simulation otherwise.
short_run = "SICONOS_EXAMPLES_SHORT_RUN" in os.environ


# import time

# print("Script démarré, en attente...")
# time.sleep(10)

# A demonstration of how to couple the two free axes of a
# CylindricalJointR in order to construct a screw relation. (Coupled
# rotational and translational motion.)

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Bouncy contact with the ground
    io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0.6)

    # Definition of a bar
    io.add_primitive_shape("Bar", "Box", (0.2, 0.2, 1))
    io.add_object("bar", [Contactor("Bar")], [0, 0, 1], mass=1)

    # Definition of the ground
    io.add_primitive_shape("Ground", "Box", (2, 3, 0.1))
    io.add_object("ground", [Contactor("Ground")], [0, 0, -0.05])

    # Add a cylindrical joint with a coupling between its two degrees
    # of freedom with a ratio of 5.0 (rotation of 5 radians for every
    # translation of 1.0 units)
    io.add_joint(
        name="joint1",
        object1="bar",
        points=[[0, 0, 0]],
        axes=[[0, 0, 1]],
        joint_class="CylindricalJointR",
        coupled=[(0, 1, 5.0)],
        absolute=True,
    )

options = sn.SolverOptions(sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12


# Load and run the simulation
with MechanicsHdf5Runner(mode="r+") as io:
    io.run(
        t0=0,
        T=0.5 if short_run else 3.0,
        h=0.001,
        theta=0.5,
        Newton_max_iter=1,
        solver_options=options,
        projection_itermax=3,
        projection_tolerance=1e-5,
        projection_tolerance_unilateral=1e-5,
        time_stepping=TimeSteppingDirectProjection,
        osi=MoreauJeanDirectProjectionOSI,
    )
