#!/usr/bin/env python
# The Dzhanibekov effect
# http://mathoverflow.net/questions/81960/the-dzhanibekov-effect-an-exercise-in-mechanics-or-fiction-explain-mathemat
#
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.mechanics.collision.tools import Contactor

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 0.01

with MechanicsHdf5Runner() as io:
    # Definition of a tetrahedron as a convex shape
    io.add_primitive_shape("Body1", "Cylinder", (1, 6))
    io.add_primitive_shape("Body2", "Box", (2, 0.4, 11))

    io.add_object(
        "roo",
        [Contactor("Body1"), Contactor("Body2", relative_translation=[0, 3, 0])],
        translation=[0, 0, 4],
        # a small perturbation on z axis
        velocity=[0, 0, 0, 0, 2, 0.0001],
        mass=1,
        inertia=[1, 10, 11],
    )
test = True
if test:
    T = 10
else:
    T = 100

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 5e-3

# run_options["theta"] = 1.0
run_options["bullet_options"] = bullet_options

run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["violation_verbose"] = False
run_options["with_timer"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = None
run_options["set_external_forces"] = lambda x: None

with MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
