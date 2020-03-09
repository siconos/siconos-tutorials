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
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cone
    R = 0.1
    L = 2.0
    io.add_primitive_shape('Con', 'Cone', (R, L))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 1.0))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.9)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    mass_test = 1.0
    inertia_test = np.eye(3)

    inertia_test[0, 0] = 0.25*mass_test*R*R + 1/3.0*mass_test*L*L
    inertia_test[1, 1] = 0.5*mass_test*R*R
    inertia_test[2, 2] = 0.25*mass_test*R*R + 1/3.0*mass_test*L*L
    print(inertia_test)
    orientation_test = [0.14, 0.7, 0.7, 0]
    io.add_object('con', [Contactor('Con')],
                  translation=[0, 0, 4],
                  orientation=orientation_test,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, -0.1])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-6

with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           face_class=None,
           edge_class=None,
           t0=0,
           T=6,
           h=0.0005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
