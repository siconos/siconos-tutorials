#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics

import math

import random

diameter =0.02

density = 1300

volume = 2./3. * math.pi * (diameter/2.0)**3

mass = volume*density
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere',
                           (diameter/2.0,),
                           insideMargin=0.0, #diameter*0.1,
                           outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (50*diameter,50*diameter , diameter))

    test=True
    if test==True:
        n_spheres = 5
    else:
        n_spheres = 1000
        
    delta_tob = math.sqrt(2.0*diameter/9.81)
    print('delta_tob', delta_tob)
    T = (n_spheres*delta_tob)*1.1
    print('T', T)
    for s in range(n_spheres):
        tob = s*delta_tob
        trans =  [0.5*diameter*random.random(), 0.5*diameter*random.random(), 15*diameter]
        print(trans)
        io.add_object('sphere_'+str(s), [Contactor('Sphere')],
                      translation=trans,
                      velocity=[0, 0, -0, 0, 0, 0],
                      mass=mass,
                      time_of_birth=tob)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[-0, 0, 0])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    #io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl('contact_rolling', e= 0.0, mu=0.3, mu_r=.5)
    #io.add_Newton_impact_friction_nsl('contact', e= 0.9, mu=0.3)

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
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=None)
