#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere', (2,),
                           insideMargin=0.2, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (100, 100, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    #io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl('contact_rolling', e= 0.9, mu=0.3, mu_r=.5)
    #io.add_Newton_impact_friction_nsl('contact', e= 0.9, mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('sphere', [Contactor('Sphere')], translation=[0, 0, 3],
                  velocity=[10, 0, 0, 0, 1, 0],
                  mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 49, 0])


# Create solver options
options = sk.solver_options_create(sn.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8
test=True
if test:
    T=1.0
    h=5e-3
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4
else:
    T=10.0
    h=5e-3
    
# # Run the simulation from the inputs previously defined and add
# # results to the hdf5 file. The visualisation of the output may be done
# # with the vview command.
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
           h=h,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
