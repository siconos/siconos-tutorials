import os
import random

import chute
import rocas
import siconos.integrators
import siconos.numerics as sn
from siconos.io.friction_contact_trace import FrictionContactTraceParams
from siconos.io.mechanics_run import MechanicsHdf5Runner

random.seed(0)

unscaled_cube_size = 0.1
unscaled_planthickness_ = unscaled_cube_size
unscaled_density = 2500

scale = 1.0 / unscaled_cube_size * 10.0

cube_size = unscaled_cube_size * scale
planethickness_ = unscaled_planthickness_ * scale
density = unscaled_density / (scale**3)

box_height = 3.683
box_length = 6.900
box_width = 3.430

planethickness_ = 0.2

test = True
if test:
    n_layer = 2
    n_row = 2
    n_col = 2
    step = 200
    hstep = 1e-3
    itermax = 1000
else:
    n_layer = 200
    n_row = 2
    n_col = 16
    step = 20000
    hstep = 1e-4
    itermax = 10000

# Create solver options
with MechanicsHdf5Runner(mode="w") as io:
    ch = chute.create_chute(
        io,
        box_height=box_height,
        box_length=box_length,
        box_width=box_width,
        planethickness_=planethickness_,
        scale=1,
        trans=[-0.6, -1.8, -1],
    )

    rcs = rocas.create_rocas(
        io,
        n_layer=n_layer,
        n_row=n_row,
        n_col=n_col,
        x_shift=2.0,
        roca_size=0.1,
        top=3,
        rate=0.02,
        density=density,
    )

    io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=0.01)


dump_probability = 0.02
theta = 1.0
tolerance = 1e-03
if not os.path.exists("Chute"):
    os.mkdir("Chute")

sid = sn.solver_ids.SICONOS_GLOBAL_FRICTION_3D_ADMM

options = sn.SolverOptions(sid)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = itermax
options.dparam[sn.params.SICONOS_DPARAM_TOL] = tolerance
fileName = "./Chute/Chute"
title = "Chute"
description = f"""
Chute with 6400 polyhedra with Bullet collision detection
Moreau TimeStepping: h={hstep}, theta = {theta}
One Step non smooth problem: {sid}, maxiter={itermax}, tol={tolerance}
"""
mathInfo = ""

friction_contact_trace_params = FrictionContactTraceParams(
    dump_itermax=9000,
    dump_probability=None,
    fileName=fileName,
    title=title,
    description=description,
    mathInfo=mathInfo,
)

with MechanicsHdf5Runner(mode="r+", collision_margin=0.01) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    if test:
        io.run(
            gravity_scale=1.0 / scale,
            t0=0,
            T=step * hstep,
            h=hstep,
            multipoints_iterations=True,
            theta=1.0,
            Newton_max_iter=1,
            output_frequency=10,
            osi=siconos.integrators.MoreauJeanGOSI,
            solver_options=options,
        )
    else:
        io.run(
            gravity_scale=1.0 / scale,
            t0=0,
            T=step * hstep,
            h=hstep,
            multipoints_iterations=True,
            theta=1.0,
            Newton_max_iter=1,
            output_frequency=10,
            osi=siconos.integrators.MoreauJeanGOSI,
            solver_options=options,
            friction_contact_trace_params=friction_contact_trace_params,
        )
