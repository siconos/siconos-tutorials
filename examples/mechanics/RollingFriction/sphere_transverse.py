#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#
import siconos.numerics as sn
import siconos.simulation
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a sphere
    io.add_primitive_shape(
        "Sphere", "Sphere", (2,), insideMargin=0.2, outsideMargin=0.0
    )

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (100, 100, 0.5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    # io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl(
        "contact_rolling", e=0.9, mu=0.3, mu_r=0.5
    )
    # io.add_Newton_impact_friction_nsl('contact', e= 0.9, mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object(
        "sphere",
        [Contactor("Sphere")],
        translation=[0, 0, 3],
        velocity=[10, 10, 0, 1, 10, 0],
        mass=1,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, 49, 0])


# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 10
run_options["h"] = 0.005
# run_options["theta"] = 1.0

run_options["solver_options"] = options

# run_options['Newton_options']=simu.LINEAR
run_options["Newton_options"] = siconos.simulation.NONLINEAR
run_options["Newton_max_iter"] = 20


run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = False


run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = None


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
