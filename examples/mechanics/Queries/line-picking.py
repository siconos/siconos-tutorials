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
# Example of how to "pick" objects, i.e., identifying objects based on
# a line intersection test.  Can be used to implement mouse
# interaction.  Here we just detect objects that pass through the
# vertical line at (x,y)=(0,0).
#

import numpy as np
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.mechanics.collision.tools import Contactor

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a sphere
    io.add_primitive_shape("Sphere", "Sphere", (0.5,))

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (50, 10, 0.1))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0.1)

    # Define a series of 10 spheres with an initial velocity to throw
    # them towards 0,0 so that we can detect them as they roll by.
    for i in range(10):
        io.add_object(
            f"sphere{i}",
            [Contactor("Sphere")],
            translation=[i * 2, 0, 2],
            velocity=[-5, 0, 0, 0, 0, 0],
            mass=1,
        )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.1])


class LinePicker:
    def initialize(self, io):
        self.io = io

    def step(self):
        start = np.zeros(3)
        start[:] = [0, 0, 10]
        end = np.zeros(3)
        end[:] = [0, 0, -10]
        query = self.io._interman.lineIntersectionQuery(start, end, False, True)
        names = [self.io._nsds.name(q.body) for q in query if q.body is not None]
        if len(names) > 0:
            print(
                "time {:f}, picked {}".format(self.io.current_time(), ", ".join(names))
            )


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-5

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0.0
run_options["T"] = 10
run_options["h"] = 0.005
run_options["theta"] = 0.50001
run_options["controller"] = LinePicker()
run_options["Newton_max_iter"] = 2
run_options["solver_options"] = options
run_options["bullet_options"] = bullet_options

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
