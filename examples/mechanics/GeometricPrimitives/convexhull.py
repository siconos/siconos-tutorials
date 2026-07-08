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

import numpy
import math

import random
from siconos.mechanics.collision.tools import Contactor
import siconos.mechanics.collision.convexhull
from siconos.mechanics.collision.bullet import SiconosBulletOptions

from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)

import siconos.numerics as sn

unscaled_polyhedron_size = 0.5
unscaled_density = 2500

scale = 1.0 / unscaled_polyhedron_size * 1000.0

polyhedron_size = unscaled_polyhedron_size * scale
density = unscaled_density / (scale**3)

box_height = 3.683 * scale
box_width = 3.430 * scale

# create some bodies

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    angle = 0.0  # math.pi/4.0
    translation = [0.0, 0.0, 0.0]
    orientation = [math.cos(angle / 2.0), 0.0, math.sin(angle / 2.0), 0.0]
    orientation = [1.0, 0.0, 0.0, 0.0]
    # orientation = [0.0, 1.0, 0.0, 0.0]
    # --> cos(\theta/2.0) = 0, sin(\theta/2.0) =1 ==> \theta = pi

    # ######## left_up
    v1 = numpy.array([0, 0, 1.0 * box_height])
    v2 = numpy.array([box_width, box_width, 0.0])
    v3 = numpy.array([box_width, -box_width, 0.0])
    v4 = numpy.array([-box_width, -box_width, 0.0])
    v5 = numpy.array([-box_width, box_width, 0.0])

    left_up_vertices = numpy.array([v1, v2, v3, v4, v5])
    print(left_up_vertices)
    io.add_convex_shape("Left_up", left_up_vertices)
    io.add_object("left_up", [Contactor("Left_up")], translation=translation)

    # ######## right_up
    io.add_convex_shape("Right_up", left_up_vertices)
    translation = [0.0, 2.0 * box_width, 0.0]
    io.add_object("right_up", [Contactor("Right_up")], translation=translation)

    n_polyhedron = 1
    n_row = 1
    n_col = 1
    x_shift = 3.0

    angle = math.pi / 2.0
    orientation_polyhedron = [math.cos(angle / 2.0), math.sin(angle / 2.0), 0.0, 0.0]
    orientation_polyhedron = [math.cos(angle / 2.0), 0.0, math.sin(angle / 2.0), 0.0]
    # orientation_polyhedron = [1.0, 0.0, 0.0, 0.0]
    # orientation_polyhedron = [0.0, 1.0, 0.0, 0.0]
    # --> cos(\theta/2.0) = 0, sin(\theta/2.0) =1 ==> \theta = pi

    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_polyhedron):
                polyhedron_size_rand = polyhedron_size * (
                    1.0 + 0.5 * (random.random() - 1.0)
                )
                polyhedron_vertices = [
                    (-polyhedron_size_rand, polyhedron_size_rand, 0.0),
                    (-polyhedron_size_rand, -polyhedron_size_rand, 0.0),
                    (polyhedron_size_rand, -polyhedron_size_rand, 0.0),
                    (polyhedron_size_rand, polyhedron_size_rand, 0.0),
                    (0.0, 0.0, 10.0 * polyhedron_size_rand),
                ]
                print("polyhedron_vertices", polyhedron_vertices)
                ch = siconos.mechanics.collision.convexhull.ConvexHull(
                    polyhedron_vertices
                )
                cm = ch.centroid()
                print("cm", cm)

                # correction of vertices such that o is the centroid
                polyhedron_vertices = numpy.array(polyhedron_vertices)[:] - cm[:]
                print("corrected polyhedron_vertices", polyhedron_vertices)
                ch = siconos.mechanics.collision.convexhull.ConvexHull(
                    polyhedron_vertices
                )
                cm = ch.centroid()
                print("cm", cm)

                # computation of inertia and volume
                inertia, volume = ch.inertia(cm)
                print("inertia,volume", inertia, volume)
                mass = volume * density
                inertia = inertia * density
                print("inertia,mass", inertia, mass)

                # Definition of a polyhedron as a convex shape
                io.add_convex_shape(
                    "PolyhedronCS" + str(n) + "_" + str(i) + "_" + str(j),
                    polyhedron_vertices,
                )

                trans = [
                    box_width / 5.0 + i * x_shift * polyhedron_size,
                    x_shift * (j + 2) * polyhedron_size,
                    box_height + polyhedron_size * x_shift * n,
                ]
                # trans= [0,0,0]
                io.add_object(
                    "polyhedron" + str(n) + "_" + str(i) + "_" + str(j),
                    [Contactor("PolyhedronCS" + str(n) + "_" + str(i) + "_" + str(j))],
                    translation=trans,
                    orientation=orientation_polyhedron,
                    velocity=[0, 0, 0, 0, 0, 0],
                    mass=mass,
                    inertia=inertia,
                )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.3)

step = 10000
hstep = 1e-4

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3


options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = step * hstep
run_options["h"] = hstep

# run_options['bullet_options']=bullet_options
run_options["solver_options"] = options
# run_options['constraint_activation_threshold']=1e-05

run_options["gravity_scale"] = 1.0 / scale

run_options["Newton_max_iter"] = 10
run_options["output_frequency"] = 100

run_options["verbose"] = False
run_options["with_timer"] = False
run_options["violation_verbose"] = True

with MechanicsHdf5Runner(mode="r+", collision_margin=0.01) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
