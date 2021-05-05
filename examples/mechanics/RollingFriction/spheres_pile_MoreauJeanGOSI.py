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

density = 1300

volume = 2./3. * math.pi * (diameter/2.0)**3

mass = volume*density

margin_ratio=1e-05

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere',
                           (diameter/2.0,),
                           insideMargin=diameter*margin_ratio,
                           outsideMargin=diameter*margin_ratio)

    # Definition of the ground shape

    box_x_scale=30*diameter
    box_y_scale=30*diameter
    box_z_scale=2*diameter


    # We use a convex hull rather than a box primitive. With the box primitive, the outside margin
    # are not taken into account. With multipoints_iterations=False, i.e, only one contact point
    # the spheres are slowly penetrating the ground. This is a well-known issue with bullet collision
    # engine between Sphere and Box

    # io.add_primitive_shape('Ground', 'Box', (50*diameter,50*diameter , diameter))


    vertices = [(0*box_x_scale,0*box_y_scale,0*box_z_scale),
                (0*box_x_scale,1*box_y_scale,0*box_z_scale),
                (1*box_x_scale,0*box_y_scale,0*box_z_scale),
                (1*box_x_scale,1*box_y_scale,0*box_z_scale),
                (0*box_x_scale,0*box_y_scale,1*box_z_scale),
                (0*box_x_scale,1*box_y_scale,1*box_z_scale),
                (1*box_x_scale,0*box_y_scale,1*box_z_scale),
                (1*box_x_scale,1*box_y_scale,1*box_z_scale)]

    io.add_convex_shape('Ground', vertices, outsideMargin=diameter*margin_ratio)

    test=True
    if test==True:
        n_spheres = 10
    else:
        n_spheres = 250
        
    delta_tob = math.sqrt(2.0*diameter/9.81)
    print('delta_tob', delta_tob)
    T = (n_spheres*delta_tob)*1.1
    print('T', T)
    for s in range(n_spheres):
        tob = s*delta_tob
        trans =  [0.5*diameter*random.random(), 0.5*diameter*random.random(), 7*diameter]
        io.add_object('sphere_'+str(s), [Contactor('Sphere')],
                      translation=trans,
                      velocity=[0, 0, -0, 0, 0, 0],
                      mass=mass,
                      time_of_birth=tob)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[-box_x_scale/2.0, -box_y_scale/2.0, 0])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    #io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl('contact_rolling', e= 0.0, mu=0.3, mu_r=1e-03)
    #io.add_Newton_impact_friction_nsl('contact', e= 0.0, mu=0.3)

# Create solver options
options = sk.solver_options_create(sn.SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR)
#options = sk.solver_options_create(sn.SICONOS_GLOBAL_FRICTION_3D_ADMM)
#options = sk.solver_options_create(sn.SICONOS_GLOBAL_FRICTION_3D_NSGS_WR)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-04
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 20
# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           multipoints_iterations=False,
           gravity_scale=1,
           t0=0,
           T=T,
           h=1e-3,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           osi=sk.MoreauJeanGOSI,
           solver_options=options,
           violation_verbose=False,
           numerics_verbose=False,
           numerics_verbose_level=0,
           output_frequency=50,
           constraint_activation_threshold=1e-08)
