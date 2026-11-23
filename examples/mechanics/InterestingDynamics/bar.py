#!/usr/bin/env python

import math
import pickle

# import scipy.constants as constants
import numpy as np
import siconos.modeling as sm
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

# WARNING : in 3D by default z-axis is upward
# this is very important to direct PLANx objects

dim = 3

unscaled_bar_length = 1.5
aspect_ratio = 100.0
unscaled_bar_height = unscaled_bar_length / aspect_ratio
unscaled_bar_width = unscaled_bar_length / aspect_ratio

unscaled_volume = unscaled_bar_length * unscaled_bar_height * unscaled_bar_width
unscaled_density = 1000

unscaled_mass = unscaled_volume * unscaled_density

print("unscaled_mass", unscaled_mass)

scale = 1.0 / unscaled_bar_length * 1.0

density = unscaled_density / (scale**3)

bar_height = unscaled_bar_height * scale
bar_length = unscaled_bar_length * scale
bar_width = unscaled_bar_width * scale

body_collection = {}
body_collection["plan_id"] = {}
id_plan = 0
# scale =1
# mass :3.375000e-01
# Inertia :
# 3.600000e-04, 0.000000e+00, 0.000000e+00,
# 0.000000e+00, 1.195312e-01, 0.000000e+00,
# 0.000000e+00, 0.000000e+00, 1.195312e-01,

# create some bodies
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    volume = bar_height * bar_length * bar_width
    mass = volume * density
    print("mass", mass)
    print("scale", scale)
    # raw_input()
    # Definition of a cube as a convex shape
    io.add_convex_shape(
        "Bar",
        [
            (-bar_length, bar_width, -bar_height),
            (-bar_length, -bar_width, -bar_height),
            (-bar_length, -bar_width, bar_height),
            (-bar_length, bar_width, bar_height),
            (bar_length, bar_width, bar_height),
            (bar_length, bar_width, -bar_height),
            (bar_length, -bar_width, -bar_height),
            (bar_length, -bar_width, bar_height),
        ],
    )

    angle = math.pi / 4.0
    ori = [math.cos(angle / 2.0), 0.0, math.sin(angle / 2.0), 0]
    print("ori initial", ori)
    io.add_object(
        "bar",
        [Contactor("Bar")],
        translation=[0, 0, 3.0 * scale],
        orientation=ori,
        velocity=[0, 0, 0, 0, 0.0, 0],
        mass=mass,
    )

    # Definition of the ground shape
    io.add_primitive_shape("Ground", "Box", (5 * scale, 5 * scale, 0.1 * scale))
    angleground = 0.0
    math.pi / 4.0
    origround = [math.cos(angleground / 2.0), 0.0, math.sin(angleground / 2.0), 0]
    io.add_object(
        "ground", [Contactor("Ground")], translation=[0, 0, 0.0], orientation=origround
    )

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl("contact", mu=0.3)

    print(body_collection)

    with open("body_collection.dict", "wb") as f:
        pickle.dump(body_collection, f)


step = 1000
hstep = 0.0001

gravity_scale = 1.0 / scale


body_weight = []
body_mext = []


def apply_forces(body):
    # g = constants.g / gravity_scale
    # weight = np.asarray([0, 0, -body.scalarMass() * g])
    # body.setConstantFext(weight)
    mext = [0.0, 0.0, 100.0]
    constant = 1e6  # unstable for 1e8
    mext = np.asarray(
        [
            scale * constant * math.cos(angle) / (aspect_ratio**2) * 0.90,
            0,
            scale * constant * math.sin(angle) / (aspect_ratio**2),
        ]
    )
    body.setConstantMext(mext, sm.alias_t)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-14


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = step * hstep
run_options["h"] = hstep
run_options["theta"] = 1.0

run_options["solver_options"] = options
run_options["Newton_max_iter"] = 10

run_options["verbose"] = True
run_options["violation_verbose"] = True
run_options["with_timer"] = True

run_options["numerics_verbose"] = False
run_options["numerics_verbose_level"] = 0

run_options["output_frequency"] = 10
run_options["gravity_scale"] = 1.0 / scale
run_options["multipoints_iterations"] = True
run_options["set_external_forces"] = apply_forces

with MechanicsHdf5Runner(mode="r+", collision_margin=0.01) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
