#!/usr/bin/env python

#
# Example of delayed object introduction with time_of_birth parameter
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cube as a convex shape
    io.add_convex_shape('Cube', [
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0)])

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (100, 100, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # The cube objects are made with an unique Contactor : the cube shape.
    # As a mass is given, they are dynamic systems involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0

    # A first cube is introduced a the beginning of the simulation
    io.add_object('cube0', [Contactor('Cube')], translation=[0, 0, 2],
                  velocity=[10, 0, 0, 1, 1, 1],
                  mass=1)

    # the second cube introduction is delayed. It is crearted in the simulation
    # a time 0.5
    io.add_object('cube1', [Contactor('Cube')], translation=[0, 0, 2],
                  velocity=[10, 0, 0, 1, 1, 1],
                  mass=1, time_of_birth=0.5)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, 0])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=10,
           h=0.0005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
