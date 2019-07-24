#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor
from mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import numpy as np

options = SiconosBulletOptions()
options.worldScale = 1.0
options.contactBreakingThreshold = 0.04
options.dimension = 1

density=1000.0


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    
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
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(verbose=True,
        with_timer=False,
           options=options,
           face_class=None,
           edge_class=None,
           t0=0,
           T=10.,
           h=0.001,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,  
           itermax=1000,
           tolerance=1e-8,
           numerics_verbose=True,
           output_frequency=1)
