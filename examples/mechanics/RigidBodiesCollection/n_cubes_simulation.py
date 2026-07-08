#!/usr/bin/env python
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
import siconos.numerics as sn
import siconos.modeling as sm

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
    n_cube = 4
    n_row = 4
    n_col = 4
    # Creation of the hdf5 file for input/output
    with MechanicsHdf5Runner() as io:
        for i in range(n_row):
            for j in range(n_col):
                for n in range(n_cube):
                    # Definition of a cube as a convex shape
                    io.add_convex_shape('CubeCS' + str(n) + '_' + str(i) + '_' + str(j), [(-1.0, 1.0, -1.0),
                                                                                          (-1.0,
                                                                                           - 1.0, -1.0),
                                                                                          (-1.0,
                                                                                           - 1.0, 1.0),
                                                                                          (-1.0,
                                                                                           1.0, 1.0),
                                                                                          (1.0, 1.0, 1.0),
                                                                                          (1.0, 1.0, -1.0),
                                                                                          (1.0,
                                                                                           - 1.0, -1.0),
                                                                                          (1.0, -1.0, 1.0)])

        # Alternative to the previous convex shape definition.
        # io.add_primitive_shape('CubePrim', 'Box', (2, 2, 2))

        # Definition of the ground shape
        io.add_primitive_shape('Ground', 'Box', (200, 200, .5), insideMargin=0.04)

        # Definition of the left shape
        # io.add_primitive_shape('Left', 'Box', (100, 0.5, 50.))

        # Definition of the right shape
        # io.add_primitive_shape('Right', 'Box', (100, 0.5, 50.))

        # Definition of the rear shape
        # io.add_primitive_shape('Rear0', 'Box', (0.5, 100., 50.))

        # Definition of the front shape
        # io.add_primitive_shape('Front', 'Box', (100, 0.5, 50.))

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
                    io.add_object('cubeCS' + str(n) + '_' + str(i) + '_' + str(j), [Contactor('CubeCS' + str(n) + '_' + str(i) + '_' + str(j))],
                                  translation=[3.0 * i, 3.0 * j, 2.05 * (n + 1)],
                                  velocity=[10 * (1.0 + 2.0 * (random.random() - 1.0) / 2.0),
                                            10 * (1.0 + 2.0 * (random.random() - 1.0) / 2.0), 0, 1, 1, 1],
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

test = True
if test:
    nstep = 100
    step = 0.0005
else:
    nstep = 20000
    step = 0.0005

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = nstep * step
run_options["h"] = step

run_options["solver_options"] = options
run_options["multipoints_iterations"] = True
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = False

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

run_options["output_frequency"] = 100


with MechanicsHdf5Runner(mode='r+', io_filename=filename) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
