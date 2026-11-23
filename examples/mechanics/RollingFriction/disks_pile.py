#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#
import math
import random

import siconos.mechanics.collision.bullet
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

diameter = 0.02
depth = 0.1

density = 1300

volume = math.pi * (diameter / 2.0) ** 2 * depth

mass = volume * density

inertia = 1 / 4.0 * mass * (diameter / 2.0) ** 2
margin_ratio = 1e-05

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a sphere
    io.add_primitive_shape(
        "Disk",
        "Disk",
        (diameter / 2.0,),
        insideMargin=diameter * margin_ratio,
        outsideMargin=diameter * margin_ratio,
    )

    # Definition of the ground shape

    box_x_scale = 30 * diameter
    box_y_scale = 2 * diameter

    # We use a convex hull rather than a box primitive.
    # With the box primitive, the outside margin
    # are not taken into account. With multipoints_iterations=False,
    # i.e, only one contact point
    # the spheres are slowly penetrating the ground.
    # This is a well-known issue with bullet collision
    # engine between Sphere and Box

    # io.add_primitive_shape('Ground', 'Box', (50*diameter,50*diameter , diameter))

    vertices = [
        (0 * box_x_scale, 0 * box_y_scale),
        (0 * box_x_scale, 1 * box_y_scale),
        (1 * box_x_scale, 0 * box_y_scale),
        (1 * box_x_scale, 1 * box_y_scale),
    ]

    io.add_convex_shape("ConvexHull", vertices, outsideMargin=diameter * margin_ratio)

    test = True
    if test is True:
        n_disks = 10
    else:
        n_disks = 250

    delta_tob = math.sqrt(2.0 * diameter / 9.81)
    print("delta_tob", delta_tob)
    T = (n_disks * delta_tob) * 1.1
    print("T", T)
    for s in range(n_disks):
        tob = s * delta_tob
        trans = [0.5 * diameter * random.random(), 12 * diameter]
        io.add_object(
            "disk_" + str(s),
            [Contactor("Disk")],
            translation=trans,
            velocity=[0, 0, 0],
            mass=mass,
            inertia=inertia,
            time_of_birth=tob,
        )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object(
        "ground",
        [Contactor("ConvexHull")],
        translation=[-box_x_scale / 2.0, -box_y_scale / 2.0],
    )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    # io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_rolling_friction_nsl(
        "contact_rolling", e=0.0, mu=0.3, mu_r=1e-03
    )
    # io.add_Newton_impact_friction_nsl('contact', e= 0.0, mu=0.3)

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_3D_NSGS)
# options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4
# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

bullet_options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
bullet_options.worldScale = 1000.0
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = siconos.mechanics.collision.bullet.TwoD

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
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

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = 10

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
