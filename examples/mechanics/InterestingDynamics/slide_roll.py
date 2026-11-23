#!/usr/bin/env python

# Various object types sliding, rolling, and sitting still.
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.mechanics.collision.tools import Contactor

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0

# Control the number of perturbations applied to generate multipoint
# surface-surface contact manifolds.  Default is 5 and 5.
bullet_options.perturbationIterations = 5
bullet_options.minimumPointsPerturbationThreshold = 5

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    # Definition of a tetrahedron as a convex shape.
    # Bottom purposely not even.
    import numpy

    pts = numpy.array(
        [(-1.0, 1.0, -1.0), (1.0, -1.0, -0.5), (-1.0, -1.0, -0.7), (0.0, 0.0, 1.0)]
    )
    io.add_convex_shape("Tetra", pts - pts.mean(0))

    io.add_primitive_shape("Cyl", "Cylinder", [1, 1])

    io.add_primitive_shape("Cube", "Box", [1, 1, 1])

    io.add_primitive_shape("Ball", "Sphere", [1])

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (100, 100, 1))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl(
        "contact", mu=0.3, e=0.0, collision_group1=1, collision_group2=2
    )

    # computation of inertia and volume
    ch = ConvexHull(pts)
    inertia, volume = ch.inertia(ch.centroid())

    # Copies of each object type, still and thrown horizontally.
    x = -20
    y = -40
    vel = 20
    spacing = 4

    io.add_object(
        "tetra1",
        [Contactor("Tetra", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
        inertia=inertia,
    )

    x += spacing
    io.add_object(
        "tetra2",
        [Contactor("Tetra", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, vel, 0, 0, 0, 0],
        mass=1,
        inertia=inertia,
    )

    x += spacing
    io.add_object(
        "cyl1",
        [Contactor("Cyl", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "cyl2",
        [Contactor("Cyl", collision_group=1)],
        translation=[x, y, 2],
        orientation=[(1, 1, 1), numpy.pi / 4],
        velocity=[0, 0, 0, 0, 0, 1],
        mass=1,
    )

    x += spacing
    io.add_object(
        "cyl3",
        [Contactor("Cyl", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, vel, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "cyl4",
        [Contactor("Cyl", collision_group=1)],
        translation=[x, y, 2],
        orientation=[(0, 0, 1), numpy.pi / 2],
        velocity=[0, vel, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "ball1",
        [Contactor("Ball", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "ball2",
        [Contactor("Ball", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, vel, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "box1",
        [Contactor("Cube", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    io.add_object(
        "box2",
        [Contactor("Cube", collision_group=1)],
        translation=[x, y, 2],
        velocity=[0, vel, 0, 0, 0, 0],
        mass=1,
    )

    x += spacing
    stack_height = 3
    for i in range(stack_height):
        io.add_object(
            f"box{i + 3}",
            [Contactor("Cube", collision_group=1 + (i % 2))],
            translation=[x, y, i + 1],
            velocity=[0, 0, 0, 0, 0, 0],
            mass=1,
        )

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object(
        "ground", [Contactor("Ground", collision_group=2)], translation=[0, 0, 0]
    )

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
test = True
if test:
    T = 0.5
else:
    T = 20.0

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.005


run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = False

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = None

with MechanicsHdf5Runner(mode="r+") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    # print(pydoc.render_doc(io.run, "Help on %s"))

    io.run(run_options)
