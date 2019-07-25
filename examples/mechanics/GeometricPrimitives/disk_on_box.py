#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics
from siconos.mechanics.collision.bullet import SiconosBulletOptions

options = SiconosBulletOptions()
options.worldScale = 1.0
options.contactBreakingThreshold = 0.04
options.dimension = 1

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Disk', 'Disk', (2,),
                           insideMargin=0.0, outsideMargin=0.0)
    
    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box2d', (10, 1),
                           insideMargin=0.0, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.5)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('disk', [Contactor('Disk')],
                  translation=[0, 5.],
                  velocity=[0, 0, 0.5],
                  mass=1., inertia =2.0)
    
    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, -.5])

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
           T=8,
           h=0.001,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8,
           numerics_verbose=True,
           output_frequency=None)
