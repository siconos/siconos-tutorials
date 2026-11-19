#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#
import math
import os

import siconos.modeling as sm
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# ctest sets SICONOS_EXAMPLES_SHORT_RUN: shorter simulation for the tests,
# full simulation otherwise.
short_run = "SICONOS_EXAMPLES_SHORT_RUN" in os.environ

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a sphere
    io.add_primitive_shape(
        "Sphere", "Sphere", (2,), insideMargin=0.2, outsideMargin=0.0
    )

    # Definition of the ground shape
    io.add_primitive_shape(
        "Ground", "Box", (10, 10, 0.1), insideMargin=0.05, outsideMargin=0.0
    )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.0)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object(
        "sphere",
        [Contactor("Sphere")],
        translation=[0, 0, 2.5],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    io.add_object(
        "sphere2",
        [Contactor("Sphere")],
        translation=[-2.5, 0.0, 5.5],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    io.add_object(
        "sphere3",
        [Contactor("Sphere")],
        translation=[1.5, 0.0, 7.5],
        velocity=[0, 0, 0, 0, 0, 0],
        time_of_birth=0.02,
        time_of_death=1.5,
        mass=1,
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[2, 0, 0.3])

    angle = math.pi / 2.0
    io.add_object(
        "ground-tob",
        [Contactor("Ground")],
        translation=[2.1, 0, -0.3],
        orientation=[math.cos(angle / 2.0), 0, math.sin(angle / 2.0), 0.0],
        time_of_birth=0.01,
        time_of_death=0.7,
    )

    angle = -math.pi / 3.0
    io.add_object(
        "ground-tob-2",
        [Contactor("Ground")],
        translation=[3.0, 0, -0.3],
        orientation=[math.cos(angle / 2.0), 0, math.sin(angle / 2.0), 0.0],
        time_of_birth=0.55,
        time_of_death=2.0,
    )

angle = math.pi / 4.0


def apply_gravity(body):
    g = 9.81
    weight = [
        body.scalarMass * g * math.sin(angle),
        0.0,
        -body.scalarMass * g * math.cos(angle),
    ]
    body.setConstantFext(weight, sm.copy_t)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.perturbationIterations = 0
bullet_options.minimumPointsPerturbationThreshold = 0

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 2.1 if short_run else 4.0
run_options["h"] = 1e-3


run_options["solver_options"] = options
run_options["bullet_options"] = bullet_options
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 20
run_options["output_frequency"] = None

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

with MechanicsHdf5Runner(mode="r+", set_external_forces=apply_gravity) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
