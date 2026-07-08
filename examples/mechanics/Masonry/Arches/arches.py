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
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
import siconos.numerics as sn
import numpy


import pickle

brick = list(pickle.load(open("brick.dat", "rb")))
vertices = list(pickle.load(open("vertices.dat", "rb")))


def one_brick(
    io, name, cname, vertices, size, density=1, trans=None, velo=None, tob=None
):

    estimated_size = max(
        numpy.array(vertices).max(axis=0) - numpy.array(vertices).min(axis=0)
    )
    print(estimated_size)
    # scale = size / max(numpy.array(vertices).max(axis=0)
    #                         - numpy.array(vertices).min(axis=0))
    scale = 1.0
    ch = ConvexHull(vertices)
    cm_ori = ch.centroid()
    # print('cm_ori', cm_ori)
    if numpy.linalg.norm(cm_ori) <= 1e-2:
        input()
    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm_ori[:]) * scale

    ch = ConvexHull(vertices)
    # cm = ch.centroid()
    # print('cm', cm)
    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.001 * size)

    # computation of inertia and volume
    inertia, volume = ch.inertia(ch.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)

    # io.add_object(name,
    #              [Contactor(cname, relative_translation = cm_ori)],
    #              translation=-cm_ori,
    #              velocity=velo,
    #              mass=volume*density,
    #              time_of_birth=tob,
    #              inertia=inertia*density)
    io.add_object(
        name,
        [
            Contactor(
                cname,
                relative_orientation=[1.0, 0.0, 0.0, 0.0],
            )
        ],
        translation=cm_ori,
        orientation=[1.0, 0.0, 0.0, 0.0],
        velocity=velo,
        mass=volume * density,
        time_of_birth=tob,
        inertia=inertia * density,
    )


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    k = 0

    vertices_array = {}

    i = 0
    for v in vertices:
        # print(v['number'])
        vertices_array[v["number"]] = i
        i = i + 1
    # print(vertices_array)
    # print(brick)
    for b in brick:
        v = []
        for vb in b:
            # print('vb',vb)
            # print(len(vertices))
            v.append(vertices[vertices_array[vb]]["coord"])
        # print('v',v)
        name = "brick%03d" % k
        cname = "brick_shp%03d" % k
        one_brick(
            io,
            name,
            cname,
            v,
            1.0,
            density=2300,
            trans=[0.0, 0.0, 0.0],
            velo=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            tob=0.0,
        )
        k = k + 1
        # if (k > 20):
        #     break
        # input()

    # Definition of the ground
    io.add_primitive_shape("Ground", "Box", (10, 1, 0.1))
    io.add_object("ground", [Contactor("Ground")], [2.5, 0, -0.1])

    # # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[25,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.6, e=0.0)

T = 3.0
# T = 3e-2
h_step = 5e-3

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 100

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = h_step

run_options["solver_options"] = options
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = True

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

run_options["output_frequency"] = 1


# Load and run the simulation
with MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
