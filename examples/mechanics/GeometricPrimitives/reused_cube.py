#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cube
    io.add_primitive_shape('Cube', 'Box', (2, 2, 2),
                           insideMargin=0.04, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 0.1),
                           insideMargin=0.04, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.01, e=0.4)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('twocubes', [Contactor('Cube',
                                         relative_translation=[1, 0.5, 1]),
                               Contactor('Cube',
                                         relative_translation=[-1, -0.5, -1])],
                  translation=[0, 0, 4],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1)

    io.add_object('cube2', [Contactor('Cube')],
                  translation=[0, 0, 6],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, time_of_birth=3)

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
bullet_options.contactBreakingThreshold = 0.01


options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-5


test= True
if test:
    T=4.0
else:
    T=10.0


with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    #import pydoc
    #print(pydoc.render_doc(io.run, "Help on %s"))

    io.run(with_timer=False,
           t0=0,
           T=T,
           h=0.005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
