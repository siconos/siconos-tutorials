#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor
import siconos.mechanics.collision.bullet
import siconos.numerics as sn
import siconos.simulation

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape("Disk", "Disk", (2,), insideMargin=0.0, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape(
        "Ground", "Box2d", (10, 1), insideMargin=0.0, outsideMargin=0.0
    )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl("contact", mu=0.1, mu_r=0.1, e=0.5)
    # io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.5)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object(
        "disk",
        [Contactor("Disk")],
        translation=[0, 5.0],
        velocity=[0, 0, 2.5],
        mass=1.0,
        inertia=2.0,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, -0.5])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

bullet_options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = siconos.mechanics.collision.bullet.TwoD
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3


options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_2D_NSGS)

# options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 6
run_options["h"] = 1e-3
# run_options["theta"] = 1.0
run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options

# run_options['Newton_options']=siconos.simulation.LINEAR
run_options["Newton_options"] = siconos.simulation.NONLINEAR
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = False

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

run_options["output_frequency"] = None


with MechanicsHdf5Runner(mode="r+") as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
