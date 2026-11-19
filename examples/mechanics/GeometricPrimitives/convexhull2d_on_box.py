#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

import os

import numpy as np
import siconos.mechanics.collision.bullet
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor

# ctest sets SICONOS_EXAMPLES_SHORT_RUN: shorter simulation for the tests,
# full simulation otherwise.
short_run = "SICONOS_EXAMPLES_SHORT_RUN" in os.environ


density = 1000.0


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    #  # Definition of a convex hull
    vertices = np.array([[0.0, 0.0], [2.0, 0.0], [0, 2.0]])

    # computation of the centroid
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # move the vertices to center the center of mass at 0.0
    vertices = np.array(vertices)[:] - cm[:]
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # computation of inertia and volume
    inertia, area = ch2d.inertia(cm)
    mass = area * density
    inertia = inertia * density
    print("inertia,mass", inertia, mass)

    io.add_convex_shape("ConvexHull", vertices)
    io.add_object(
        "convexhull",
        [Contactor("ConvexHull")],
        translation=[0, 1.0],
        velocity=[0, 0, 0.0],
        mass=mass,
        inertia=inertia,
    )

    # Definition of a second convex hull
    vertices = [[0.0, 0.0], [2.0, 0.0], [0, 2.0], [1, 3.0], [3, 2.0]]

    # computation of the centroid
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # move the vertices to center the center of mass at 0.0
    vertices = np.array(vertices)[:] - cm[:]
    ch2d = ConvexHull2d(vertices)
    cm = ch2d.centroid()

    # computation of inertia and volume
    inertia, area = ch2d.inertia(cm)
    mass = area * density
    inertia = inertia * density
    print("inertia,mass", inertia, mass)

    io.add_convex_shape("ConvexHull2", vertices)
    io.add_object(
        "convexhull2",
        [Contactor("ConvexHull2")],
        translation=[0, 9.0],
        velocity=[0, 0, 0.0],
        mass=mass,
        inertia=inertia,
    )

    # Definition of the ground shape
    io.add_primitive_shape(
        "Ground", "Box2d", (15, 1), insideMargin=0.0, outsideMargin=0.0
    )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object("ground", [Contactor("Ground")], translation=[0, -0.5])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.5)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

bullet_options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = siconos.mechanics.collision.bullet.TwoD
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3


options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 1.0 if short_run else 10.0
run_options["h"] = 0.001


run_options["solver_options"] = options
run_options["bullet_options"] = bullet_options
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 1
run_options["output_frequency"] = 1

run_options["verbose"] = True
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options["numerics_verbose"] = True
run_options["numerics_verbose_level"] = 0

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
