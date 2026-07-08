#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

import siconos.numerics as sn
import siconos.modeling as sm


edge_length = 0.1
plane_length = 2.0

velocity_init = -1.0
angular_velocity_init = 0.0

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # # Definition of a cube as a convex shape
    io.add_convex_shape('CubeCS1', [
        (-edge_length, edge_length, -edge_length),
        (-edge_length, -edge_length, -edge_length),
        (-edge_length, -edge_length, edge_length),
        (-edge_length, edge_length, edge_length),
        (edge_length, edge_length, edge_length),
        (edge_length, edge_length, -edge_length),
        (edge_length, -edge_length, -edge_length),
        (edge_length, -edge_length, edge_length)])

    io.add_convex_shape('CubeCS2', [
        (-edge_length, edge_length, -edge_length),
        (-edge_length, -edge_length, -edge_length),
        (-edge_length, -edge_length, edge_length),
        (-edge_length, edge_length, edge_length),
        (edge_length, edge_length, edge_length),
        (edge_length, edge_length, -edge_length),
        (edge_length, -edge_length, -edge_length),
        (edge_length, -edge_length, edge_length)])

    # Alternative to the previous convex shape definition.
    io.add_primitive_shape('CubePrim1', 'Box',
                           (2 * edge_length, 2 * edge_length, 2 * edge_length))

    io.add_primitive_shape('CubePrim2', 'Box',
                           (2 * edge_length, 2 * edge_length, 2 * edge_length))

    # Alternative to the previous convex shape definition.
    io.add_primitive_shape('SpherePrim1', 'Sphere', (edge_length,))

    io.add_primitive_shape('SpherePrim2', 'Sphere', (edge_length,))

    # select contactor

    contactors = ['CubePrim1', 'CubePrim2']
    contactors = ['CubePrim1', 'CubeCS2']
    contactors = ['CubeCS1', 'CubeCS2']

    contactors = ['SpherePrim1', 'SpherePrim2']
    test = 'Sphere_Sphere'
    contactors = ['SpherePrim1', 'CubePrim2']

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.5)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('object1', [Contactor(contactors[0])], translation=[0, 0, 2],
                  velocity=[0, 0, - velocity_init, angular_velocity_init,
                            angular_velocity_init, angular_velocity_init],
                  mass=1)

    io.add_object('object2', [Contactor(contactors[1])],
                  translation=[0, 0, 2 + 3 * edge_length],
                  velocity=[0, 0, velocity_init, angular_velocity_init,
                            angular_velocity_init, angular_velocity_init],
                  mass=1)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


test = True
if test :
    step = 2000
    hstep = 0.001
else:
    step = 20000
    hstep = 0.001


# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = step * hstep
run_options["h"] = hstep
run_options["theta"] = 1.0

run_options["solver_options"] = options
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = True
run_options["with_timer"] = False

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

run_options["output_frequency"] = 1
run_options["gravity_scale"] = 0.1
run_options["multipoints_iterations"] = True


with MechanicsHdf5Runner(mode='r+', collision_margin=0.05) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    io.run(run_options)
