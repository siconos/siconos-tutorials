#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk


import numpy as np


sphere_count = 0
cluster_count = 0


def add_sphere_cluster(io, radius, n_spheres=8, dispersion=None,
                       translation=None, orientation=None,
                       mass=1, tob=-1):
    global sphere_count
    global cluster_count
    if dispersion is None:
        dispersion = radius
    spheres = []
    locations = np.random.normal(0, dispersion, (n_spheres, 3))
    for n in range(n_spheres):
        sphere_count += 1
        spheres.append('Sphere%03d' % sphere_count)
        io.add_primitive_shape(spheres[-1], 'Sphere', [radius])

    cluster_count += 1
    io.add_object('cluster%d' % cluster_count,
                  [Contactor(sph, relative_translation=loc)
                   for sph, loc in zip(spheres, locations)],
                  translation=translation, orientation=orientation,
                  mass=mass, time_of_birth=tob)


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cube
    io.add_primitive_shape('Cube1', 'Box', (2, 2, 2),
                           insideMargin=0.04, outsideMargin=0.0)
    io.add_primitive_shape('Cube2', 'Box', (2, 2, 2),
                           insideMargin=0.04, outsideMargin=0.0)
    io.add_primitive_shape('Cube3', 'Box', (2, 2, 2),
                           insideMargin=0.04, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (20, 20, 0.1),
                           insideMargin=0.04, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.4)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('cube', [Contactor('Cube1',
                                     relative_translation=[1, 0.5, 1]),
                           Contactor('Cube2',
                                     relative_translation=[-1, -0.5, -1])],
                  translation=[0, 0, 4],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1)

    io.add_object('cube2', [Contactor('Cube3')],
                  translation=[0, 0, 6],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, time_of_birth=3)

    for i in range(10):
        add_sphere_cluster(io, radius=np.max((0.1,
                                              np.random.normal(0.3, 0.1))),
                           translation=np.random.normal(0, 2, 3)+[0, 0, 10],
                           orientation=[1, 0, 0, 0])

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, -0.1])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

test=True
if test:
    T=3.0
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
else:
    T=10.0


with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           bullet_options=bullet_options,
           t0=0,
           T=T,
           h=0.005,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
