#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import siconos.numerics as sn
import siconos.kernel as sk
import numpy as np
import math


density=1000.0

disk_radius =1.

box_height =1.
box_width =1.


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:


    # Definition of a sphere
    io.add_primitive_shape('Disk', 'Disk', (disk_radius,),
                           insideMargin=0.0, outsideMargin=0.0)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0


    mass=density*disk_radius**2 * math.pi
    inertia = 1/2. * mass * disk_radius**2
    io.add_object('disk', [Contactor('Disk')],
                  translation=[0, 5.],
                  velocity=[0, 0, 0.5],
                  mass=mass , inertia = inertia)

    # Definition of a box
    io.add_primitive_shape('Box', 'Box2d', (box_height,box_width),
                           insideMargin=0.0, outsideMargin=0.0)

    mass=density*box_height*box_width
    inertia = 1/12. * mass * ( box_height**2 + box_width**2)
    io.add_object('box', [Contactor('Box')],
                  translation=[0, 7.],
                  velocity=[0, 0, 0.5],
                  mass=mass, inertia = inertia)


    
    #  # Definition of a convex hull
    vertices = np.array([[0., 0.], [2., 0.], [0,2.]])

    # computation of the centroid
    ch2d= ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # move the vertices to center the center of mass at 0.0
    vertices = np.array(vertices)[:]-cm[:]
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()
    
    
    # computation of inertia and volume
    inertia,area=ch2d.inertia(cm)
    mass = area*density
    inertia = inertia*density
    print ('inertia,mass',inertia,mass)

    io.add_convex_shape('ConvexHull',vertices )
    io.add_object('convexhull', [Contactor('ConvexHull')],
                  translation=[0, 1.],
                  velocity=[0, 0, 0.0],
                  mass=mass, inertia = inertia)
    
    
    # Definition of a second convex hull
    vertices = [[0., 0.], [2., 0.], [0,2.], [1,3.], [3, 2.]]  

    # computation of the centroid
    ch2d= ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # move the vertices to center the center of mass at 0.0
    vertices = np.array(vertices)[:]-cm[:]
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    
    # computation of inertia and volume
    inertia,area=ch2d.inertia(cm)
    mass = area*density
    inertia = inertia*density
    print ('inertia,mass',inertia,mass)

    io.add_convex_shape('ConvexHull2',vertices )
    io.add_object('convexhull2', [Contactor('ConvexHull2')],
                  translation=[0, 9.],
                  velocity=[0, 0, 0.0],
                  mass=mass, inertia = inertia)
    

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box2d', (15, 1),
                           insideMargin=0.0, outsideMargin=0.0)
    
    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, -.5])
    
    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.5)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = 1
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3


options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

test= True
if test:
    T=2.0
else:
    T=10.0


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
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=True,
           output_frequency=1)
