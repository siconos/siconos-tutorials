#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


import os


filename = 'n_cubes_scene.hdf5'
filename_py = 'n_cubes_scene.py'

if not os.path.isfile(filename):
    print('creation of the scene that is not existing')
    import random
    from siconos.mechanics.collision.tools import Contactor
    n_cube=4
    n_row=4
    n_col=4
    # Creation of the hdf5 file for input/output
    with MechanicsHdf5Runner() as io:
        for i in range(n_row):
            for j in range(n_col):
                for n in range(n_cube):
                # Definition of a cube as a convex shape
                    io.add_convex_shape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j), [ (-1.0, 1.0, -1.0),
                                                                               (-1.0, -1.0, -1.0),
                                                                               (-1.0, -1.0, 1.0),
                                                                               (-1.0, 1.0, 1.0),
                                                                               (1.0, 1.0, 1.0),
                                                                               (1.0, 1.0, -1.0),
                                                                               (1.0, -1.0, -1.0),
                                                                               (1.0, -1.0, 1.0)])



        # Alternative to the previous convex shape definition.
        #io.add_primitive_shape('CubePrim', 'Box', (2, 2, 2))

        # Definition of the ground shape
        io.add_primitive_shape('Ground', 'Box', (200, 200, .5), insideMargin=0.04)

        # Definition of the left shape
        # io.add_primitive_shape('Left', 'Box', (100, 0.5, 50.))

        # Definition of the right shape
        #io.add_primitive_shape('Right', 'Box', (100, 0.5, 50.))

        # Definition of the rear shape
        #io.add_primitive_shape('Rear0', 'Box', (0.5, 100., 50.))

        # Definition of the front shape
        #io.add_primitive_shape('Front', 'Box', (100, 0.5, 50.))

        # Definition of a non smooth law. As no group ids are specified it
        # is between contactors of group id 0.
        io.add_Newton_impact_friction_nsl('contact', mu=0.3)

        # The cube object made with an unique Contactor : the cube shape.
        # As a mass is given, it is a dynamic system involved in contact
        # detection and in the simulation.  With no group id specified the
        # Contactor belongs to group 0
        for i in range(n_row):
            for j in range(n_col):
                for n in range(n_cube):
                    io.add_object('cubeCS'+str(n)+'_'+str(i)+'_'+str(j), [Contactor('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))],
                                 translation=[3.0*i, 3.0*j, 2.05*(n+1)],
                                 velocity=[10*(1.0+2.0*(random.random()-1.0)/2.0), 10*(1.0+2.0*(random.random()-1.0)/2.0), 0, 1, 1, 1],
                                 mass=1)

        # io.add_object('cube2', [Contactor('CubePrim')], translation=[0, 3, 2],
        #              velocity=[10, 0, 0, 1, 1, 1],
        #              mass=1)

        # the ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        io.add_object('ground', [Contactor('Ground')],
                     translation=[50, 50, 0])
        # io.add_object('left', [Contactor('Left')],
        #              translation=[0, 50., 25.])
        # io.add_object('right', [Contactor('Right')],
        #              translation=[0, -50., 25.])
        # io.add_object('rear00', [Contactor('Rear0')],
        #              translation=[25., 0., 250.])


    filename = 'n_cubes_simulation.hdf5'

test=True
if test:
    nstep=100
    step=0.0005
else:
    nstep=20000
    step=0.0005

# Create solver options
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.iparam[sn.SICONOS_DPARAM_TOL] = 1e-4
with MechanicsHdf5Runner(mode='r+',io_filename=filename) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=nstep*step,
           h=step,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=100)
