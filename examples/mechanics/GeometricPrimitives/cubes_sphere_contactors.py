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
import math

import numpy
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.mechanics.collision.tools import Contactor

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.01

planthickness_ = 0.05

density = 2679.1838


def normal_plane(p1, p2, p3):

    x1 = p1[0]
    y1 = p1[1]
    z1 = p1[2]
    x2 = p2[0]
    y2 = p2[1]
    z2 = p2[2]
    x3 = p3[0]
    y3 = p3[1]
    z3 = p3[2]

    vector1 = [x2 - x1, y2 - y1, z2 - z1]
    vector2 = [x3 - x1, y3 - y1, z3 - z1]

    cross_product = [
        vector1[1] * vector2[2] - vector1[2] * vector2[1],
        -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]),
        vector1[0] * vector2[1] - vector1[1] * vector2[0],
    ]

    a = cross_product[0]
    b = cross_product[1]
    c = cross_product[2]
    # d = -(cross_product[0] * x1 + cross_product[1] * y1 + cross_product[2] * z1)

    return numpy.array([a, b, c]) / numpy.linalg.norm([a, b, c])


# create plans

# Creation of the hdf5 file
with MechanicsHdf5Runner(use_compression=True) as io:
    # ######## amont
    v0 = numpy.array([5.00, 1.6131, 1.0751])
    v1 = numpy.array([-2.50, 1.6131, 1.0751])
    v2 = numpy.array([-2.50, 0.9354, 0.3535])
    v3 = numpy.array([5.00, 0.9354, 0.3535])

    amont_normal = normal_plane(v1, v2, v3)
    print("amont_normal=", amont_normal)

    v0_extruded = v0 + numpy.dot(planthickness_, amont_normal)
    v1_extruded = v1 + numpy.dot(planthickness_, amont_normal)
    v2_extruded = v2 + numpy.dot(planthickness_, amont_normal)
    v3_extruded = v3 + numpy.dot(planthickness_, amont_normal)

    amont_vertices = numpy.array(
        [v0, v1, v2, v3, v0_extruded, v1_extruded, v2_extruded, v3_extruded]
    )
    print("amont_vertices", amont_vertices)

    io.add_convex_shape("amont", amont_vertices)
    io.add_object("amont", [Contactor("amont")], translation=[1.50, -1.45, -1.5331])

    # ######## aval
    v4 = numpy.array([-2.50, 0, 0])
    v5 = numpy.array([5.00, 0, 0])

    aval_normal = normal_plane(v2, v4, v3)
    print("aval_normal=", aval_normal)

    v4_extruded = v4 + numpy.dot(planthickness_, aval_normal)
    v5_extruded = v5 + numpy.dot(planthickness_, aval_normal)

    aval_vertices = numpy.array(
        [v2, v3, v4, v5, v2_extruded, v3_extruded, v4_extruded, v5_extruded]
    )
    print("aval_vertices", aval_vertices)

    io.add_convex_shape("aval", aval_vertices)
    io.add_object("aval", [Contactor("aval")], translation=[1.50, -1.45, -1.5331])

    # ######## sol
    v6 = numpy.array([5.00, -5.00, 0])
    v7 = numpy.array([-2.50, -5.00, 0])

    sol_normal = normal_plane(v4, v6, v5)
    print("sol_normal=", sol_normal)

    v6_extruded = (
        v6 - [planthickness_, 0.0, 0.0] + numpy.dot(planthickness_, sol_normal)
    )
    v7_extruded = (
        v7 + [planthickness_, 0.0, 0.0] + numpy.dot(planthickness_, sol_normal)
    )

    sol_vertices = numpy.array(
        [
            v4 - [planthickness_, 0.0, 0.0],
            v5 + [planthickness_, 0.0, 0.0],
            v6 - [planthickness_, 0.0, 0.0],
            v7 + [planthickness_, 0.0, 0.0],
            v4_extruded - [planthickness_, 0.0, 0.0],
            v5_extruded + [planthickness_, 0.0, 0.0],
            v6_extruded,
            v7_extruded,
        ]
    )
    print("sol_vertices", sol_vertices)

    io.add_convex_shape("sol", sol_vertices)
    io.add_object("sol", [Contactor("sol")], translation=[1.50, -1.45, -1.5331])

    n_cube = 1
    n_row = 50
    n_col = 1
    cube_size = 0.0144
    x_shift = 0.030
    x_translate = 0.1
    sphere_count = 0
    spheres = []
    radius = 0.005
    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_cube):
                # Definition of a cube
                vertices = [
                    (-cube_size, cube_size, -cube_size),
                    (-cube_size, -cube_size, -cube_size),
                    (-cube_size, -cube_size, cube_size),
                    (-cube_size, cube_size, cube_size),
                    (cube_size, cube_size, cube_size),
                    (cube_size, cube_size, -cube_size),
                    (cube_size, -cube_size, -cube_size),
                    (cube_size, -cube_size, cube_size),
                ]
                io.add_convex_shape(
                    "CubeCS" + str(n) + "_" + str(i) + "_" + str(j), vertices
                )

                for v in vertices:
                    sphere_count += 1
                    spheres.append(f"Sphere{sphere_count:03d}")
                    io.add_primitive_shape(spheres[-1], "Sphere", [radius])

                # computation of inertia and volume
                ch = ConvexHull(vertices)
                inertia, volume = ch.inertia(ch.centroid())
                # print inertia, volume
                # raw_input()
                # angle_init = math.pi/4.0
                # angle_init = 0.001
                angle_init = 0.0
                # contactor = [Shape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))]
                contactor = []
                for sph, loc in zip(spheres, vertices):
                    contactor.append(
                        Contactor(sph, relative_translation=loc, collision_group=1)
                    )
                # add a "fake" contactor for visualization ...
                contactor.append(
                    Contactor(
                        "CubeCS" + str(n) + "_" + str(i) + "_" + str(j),
                        collision_group=-1,
                    )
                )
                io.add_object(
                    "cube" + str(n) + "_" + str(i) + "_" + str(j),
                    contactor,
                    translation=[
                        i * (x_translate + x_shift * cube_size),
                        x_shift * j * (x_translate + cube_size),
                        (x_translate + cube_size * x_shift) * n,
                    ],
                    velocity=[0, 0, 0, 0, 0, 0],
                    orientation=[
                        math.cos(angle_init / 2.0),
                        0,
                        math.sin(angle_init / 2.0),
                        0,
                    ],
                    mass=volume * density,
                    inertia=inertia * density,
                )

    # Definition of a non smooth law
    io.add_Newton_impact_friction_nsl(
        "contact", e=0.01, mu=0.9, collision_group1=0, collision_group2=1
    )


test = True
if test:
    step = 100
    hstep = 0.0005
else:
    step = 10000
    hstep = 0.0005

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

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


with MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
