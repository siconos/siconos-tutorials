
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk
# A demonstration of how to couple the two free axes of a
# CylindricalJointR in order to construct a screw relation. (Coupled
# rotational and translational motion.)

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Bouncy contact with the ground
    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.6)

    # Definition of a bar
    io.add_primitive_shape('Bar', 'Box', (0.2, 0.2, 1))
    io.add_object('bar', [Contactor('Bar')], [0,0,1], mass=1)

    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (2, 3, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0,0,-0.05])

    # Add a cylindrical joint with a coupling between its two degrees
    # of freedom with a ratio of 5.0 (rotation of 5 radians for every
    # translation of 1.0 units)
    io.add_joint('joint1', 'bar', None, [[0,0,0]], [[0,0,1]], 'CylindricalJointR',
                coupled=[(0,1,5.0)], absolute=True)
            
options = sk.solver_options_create(sn.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-12


# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=3,
           h=0.001,
           theta=0.5,
           Newton_max_iter=1,
           solver_options=options,
           projection_itermax=3,
           projection_tolerance=1e-5,
           projection_tolerance_unilateral=1e-5,
           time_stepping=sk.TimeSteppingDirectProjection,
           osi=sk.MoreauJeanDirectProjectionOSI,
    )
