#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

import math

import random

diameter =0.02
depth= 0.1

density = 1300

volume = math.pi*(diameter/2.0)**2*depth

mass = volume*density

inertia = 1/4.0 * mass * (diameter/2.0)**2
margin_ratio=1e-05

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Disk', 'Disk',
                           (diameter/2.0,),
                           insideMargin=diameter*margin_ratio,
                           outsideMargin=diameter*margin_ratio)

    # Definition of the ground shape

    box_x_scale=30*diameter
    box_y_scale=2*diameter



    # We use a convex hull rather than a box primitive. With the box primitive, the outside margin
    # are not taken into account. With multipoints_iterations=False, i.e, only one contact point
    # the spheres are slowly penetrating the ground. This is a well-known issue with bullet collision
    # engine between Sphere and Box

    # io.add_primitive_shape('Ground', 'Box', (50*diameter,50*diameter , diameter))


    vertices = [(0*box_x_scale,0*box_y_scale),
                (0*box_x_scale,1*box_y_scale),
                (1*box_x_scale,0*box_y_scale),
                (1*box_x_scale,1*box_y_scale)]

    io.add_convex_shape('ConvexHull', vertices, outsideMargin=diameter*margin_ratio)

    test=True
    if test==True:
        n_disks = 10
    else:
        n_disks = 250
        
    delta_tob = math.sqrt(2.0*diameter/9.81)
    print('delta_tob', delta_tob)
    T = (n_disks*delta_tob)*1.1
    print('T', T)
    for s in range(n_disks):
        tob = s*delta_tob
        trans =  [0.5*diameter*random.random(), 12*diameter]
        io.add_object('disk_'+str(s), [Contactor('Disk')],
                      translation=trans,
                      velocity=[0, 0, 0],
                      mass=mass,
                      inertia = inertia,
                      time_of_birth=tob)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('ConvexHull')],
                  translation=[-box_x_scale/2.0, -box_y_scale/2.0])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    #io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl('contact_rolling', e= 0.0, mu=0.3, mu_r=1e-03)
    #io.add_Newton_impact_friction_nsl('contact', e= 0.9, mu=0.3)

# Create solver options
options = sk.solver_options_create(sn.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4
# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1000.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = 1


with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           multipoints_iterations=False,
           bullet_options=bullet_options,
           gravity_scale=1,
           t0=0,
           T=T,
           h=1e-3,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver_options=options,
           violation_verbose=False,
           numerics_verbose=False,
           output_frequency=10,
           constraint_activation_threshold=1e-08)