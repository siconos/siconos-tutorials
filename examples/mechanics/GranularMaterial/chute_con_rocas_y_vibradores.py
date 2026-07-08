# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
import siconos.numerics as sn

# import chute_con_vibrator_rear_up as chute_con_vibradores
# import chute_con_vibrador_bottom as chute_con_vibradores
# import chute_con_vibrador_middle as chute_con_vibradores
import chute_con_vibrador_top as chute_con_vibradores

import rocas
import random
import sys
import numpy

if len(sys.argv) < 2:
    dist = "uniform"
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]

if dist not in ["uniform", "double", "exp"]:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 1.5:
    print("mu = [0.1 .. 1.5]")
    sys.exit(1)

fn = "chute_con_rocas_y_vibradores-{0}-mu-{1}.hdf5".format(dist, mu)

random.seed(0)
numpy.random.seed(0)

cube_size = 0.1
planthickness_ = cube_size
density = 2500

box_height = 3.683
box_length = 6.900
box_width = 3.430

planethickness_ = 0.2

test = True
if test:
    n_layer = 10
    n_row = 2
    n_col = 2
    T = 0.5
    hstep = 1e-3
else:
    n_layer = 200
    n_row = 2
    n_col = 16
    T = 45
    hstep = 1e-3


with MechanicsHdf5Runner(mode="w", io_filename=fn) as io:

    io.add_Newton_impact_friction_nsl(
        "contact_rocas", mu=float(mu), e=0.01, collision_group1=0, collision_group2=0
    )
    io.add_Newton_impact_friction_nsl(
        "contact_chute_rocas",
        mu=float(mu),
        e=0.01,
        collision_group1=1,
        collision_group2=0,
    )
    io.add_Newton_impact_friction_nsl(
        "contact_vibratores_rocas",
        mu=float(mu),
        e=0.01,
        collision_group1=2,
        collision_group2=0,
    )

    ch = chute_con_vibradores.create_chute(
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
        x_shift=1.5,
        roca_size=0.13,
        top=3,
        rate=0.2,
        density=density,
        distribution=(dist, 0.12),
    )

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = hstep

run_options["solver_options"] = options
run_options["Newton_max_iter"] = 1

run_options["verbose"] = True
run_options["verbose_progress"] = True
run_options["with_timer"] = False

run_options['numerics_verbose']=False
run_options['numerics_verbose_level']=0

run_options["output_frequency"] = 10
run_options["multipoints_iterations"] = True

with MechanicsHdf5Runner(mode="r+", io_filename=fn) as io:
    io.run(run_options)
    
