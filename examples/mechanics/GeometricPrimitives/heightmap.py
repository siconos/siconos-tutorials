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
# Example of how to use the heightmap object to interact with static terrain.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
import siconos.numerics as sn

import numpy as np

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Generate a numpy matrix defining a slope + sinusoid + noise
    noise_amp = 1
    sine_amp = 2
    sine_freq = 0.3
    slope_offset = 50
    slope_height = 50

    s = np.array([np.cos(np.linspace(0, 100, 100) * sine_freq) * sine_amp] * 100)
    heightmap = (
        np.random.uniform(0, noise_amp, size=(100, 100))
        + s
        + s.T
        + [np.linspace(slope_offset, slope_offset + slope_height, 100)] * 100
    )

    # The heightmap has tangential extents 50 by 50, and its height is
    # not scaled; actual height must be reflected by values in
    # heightmap data.
    io.add_height_map("MyTerrain", heightmap, (50, 50))

    # Definition of a sphere
    io.add_primitive_shape("Ball", "Sphere", (1,))

    # Definition of a cube
    io.add_primitive_shape("Cube", "Box", (8, 8, 8))

    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.add_Newton_impact_friction_nsl(
        "contact", mu=0.3, e=0.0, collision_group1=0, collision_group2=1
    )

    # Rain down a 2D array of spheres -- examining the initial collision
    # locations will ensure the terrain is where we think it is.
    for i in range(20):
        for j in range(20):
            io.add_object(
                "ball_%d_%d" % (i, j),
                [Contactor("Ball", collision_group=1)],
                translation=[j * 5 - 47.5, i * 5 - 47.5, 110],
                mass=1,
            )

    # Add a nice big cube for effect
    io.add_object(
        "cube", [Contactor("Cube", collision_group=1)], translation=[0, 18, 150], mass=1
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object(
        "ground", [Contactor("MyTerrain", collision_group=0)], translation=[0, 0, 0]
    )

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-5
test = True
if test:
    T = 3.0
else:
    T = 20.0

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.01


run_options["solver_options"] = options
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 1
run_options["output_frequency"] = None

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

with MechanicsHdf5Runner(mode="r+") as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.

    # We run the simulation with relatively low precision because we
    # are interested mostly just in the location of contacts with the
    # height field, but we are not evaluating the performance of
    # individual contacts here.
    io.run(run_options)
