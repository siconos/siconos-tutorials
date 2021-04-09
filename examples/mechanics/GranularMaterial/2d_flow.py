#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import numpy
import random
import math

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = 1
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3

density = 1000.0


def create_grain(io, name, cname, grain_size=0.05, density=1, trans=None,
                 tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi / 2 * random.gauss(0.5, 0.2) for _ in range(16)]

    def vert(id1, id2, a, b):
        return (a * math.cos(rd[id1]) * math.cos(rd[id2]),
                b * math.sin(rd[id1]) * math.cos(rd[id2]))

    vertices = [vert(0, 1, 1, 1),
                vert(2, 3, 1, -1),
                vert(4, 5, -1, 1),
                vert(6, 7, -1, -1)]

    #print('vertices', vertices)

    scale = grain_size / max(
        numpy.array(vertices).max(axis=0) - numpy.array(vertices).min(axis=0))

    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.001 * grain_size)

    # computation of inertia and volume
    inertia, area = ch2d.inertia(ch2d.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)
    io.add_object(name,
                  [Contactor(cname)],
                  translation=trans,
                  #velocity=veloci,
                  mass=area * density,
                  time_of_birth=tob,
                  inertia=inertia * density)


def create_grains(io, n_row=5, n_col=5, x_shift=3.0,
                  grain_size=0.05, top=0, rate=0.01, density=1,
                  distribution=('uniform', 0.1)):

    N = n_row * n_col

    dist, rng = distribution

    if dist == 'uniform':
        sizes = numpy.random.uniform(low=grain_size - rng / 2,
                                     high=grain_size + rng / 2,
                                     size=N)
    elif dist == 'double':
        sizes = numpy.hstack(
            (numpy.random.normal(scale=rng * 0.2,
                                 loc=grain_size - rng / 2,
                                 size=N / 2),
             numpy.random.normal(scale=rng * 0.2,
                                 loc=grain_size + rng / 2,
                                 size=N / 2)))
        numpy.random.shuffle(sizes)
        # Gives a rock size distribution with two sizes of rocks, with
        # the mean between both at grain_size
    elif dist == 'exp':
        # In a normal distribution, 10-
        # and 90-percentiles are loc +/- rng*1.28.
        # Let's arrange the 10- and 90-percentiles of the distribution
        # to be in the desired range.
        sizes = numpy.random.exponential(1, N)
        bottom = numpy.percentile(sizes, 10)
        top = numpy.percentile(sizes, 90)
        scale = (rng * 1.28) / (top - bottom)
        sizes = (sizes - bottom) * scale + grain_size - rng / 2 * 1.28

    k = 0
    print('Creation of the rocks')
    for i in range(n_row):
        for j in range(n_col):
            # initial translation
            if (k % 100 == 0):
                print('.', end='', flush=True)
            trans = [(i - n_row / 2.0) * x_shift * grain_size,
                     (j) * x_shift * grain_size]
            name = 'rock' + '_' + str(i) + '_' + str(j)
            cname = 'RockCS' + '_' + str(i) + '_' + str(j)
            create_grain(io, name, cname, sizes[k], density, trans,
                         tob=random.random() * rate)
            k += 1


test = True
if test:
    n_row = 5
    n_col = 3
    T = 1.0
    hstep = 1e-3
else:
    n_row = 50
    n_col = 100
    T = 5.
    hstep = 1e-3


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    grains = create_grains(io, n_row=n_row, n_col=n_col,
                           x_shift=2.0, grain_size=0.1, top=3,
                           rate=0.2, density=density)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box2d', (25, 1),
                           insideMargin=0.0, outsideMargin=0.0)

    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, -0.5 - 0.2])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.5)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8
fclocal = sk.solver_options_get_internal_solver(options, 0)
fclocal.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100

with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(verbose=True,
           with_timer=False,
           bullet_options=bullet_options,
           face_class=None,
           edge_class=None,
           t0=0,
           T=T,
           h=hstep,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=True,
           output_frequency=100,
           explode_Newton_solve=False,
           display_Newton_convergence=False,
           )
