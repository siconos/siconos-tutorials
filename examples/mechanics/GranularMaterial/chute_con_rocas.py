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
import siconos.simulation as simu
import siconos.nonsmooth_formulations as nsf
from siconos.mechanics.collision.bullet import SiconosBulletOptions
import numpy

import chute
import rocas
import random

import sys

if len(sys.argv) < 2:
    dist = "uniform"
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]

if dist not in ["uniform", "double", "exp"]:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 2.0:
    print("mu = [0.1 .. 2.0]")
    sys.exit(1)

random.seed(0)

box_height = 3.683
box_length = 6.900
box_width = 3.430


density = 2500
planethickness_ = 0.2
cube_size = 0.1

rock_shape = "roca"
# rock_shape = 'esfera'
# rock_shape = 'cubo'

# hdf5 file name
if rock_shape == "roca":
    fn = "chute_con_rocas-{0}-mu-{1}.hdf5".format(dist, mu)
else:
    fn = "chute_con_rocas-{0}-mu-{1}_{2}.hdf5".format(dist, mu, rock_shape)

test = True
if test:
    n_layer = 10
    n_row = 2
    n_col = 2
    T = 3.0
    hstep = 1e-3
else:
    n_layer = 200
    n_row = 6
    n_col = 16
    T = 20.0
    hstep = 1e-4

N = int(T / hstep)
# A hook to remove body that to far from the hopper at the end of the iteration


class death_hook:
    def __init__(self):
        pass

    def initialize(self, io):
        self._io = io
        pass

    def call(self, step):
        if step % 100 == 0:
            print("\n call death hook at step", step, "/", N)
        else:
            return
        # # First way (slow) : loop over the bodies
        # nsds = self._io._nsds
        # print(dir(nsds))
        # print(nsds.displayDynamicalSystems())
        # nds= self._io._nsds.getNumberOfDS()
        # print('nds =', nds)
        # allds=self._io._nsds.dynamicalSystemsVector()

        # for ds in allds:
        #     #ds.display()
        #     z = ds.q()[2]
        #     print('z', z)
        #     if (z <= -2):
        #         self._io._interman.removeBody(ds)
        #         self._io._nsds.removeDynamicalSystem(ds)
        #         print('remove ds number ', ds.number(), ' with height = ', z)

        # second way (faster) :  direct access to nsds positions
        positions = self._io._io.positions(self._io._nsds)
        # print(positions, positions.shape)

        height_filter = -2.0

        if positions.shape[0] > 0:
            # pass
            z = positions[3, :]

            if numpy.linalg.norm(z) > 1e10:
                print("WARNING Huge height")

            # print('z', z)
            # We search for the ds index that are below a given criteria
            ds_idx = numpy.nonzero(z < height_filter)[0]
            for i in ds_idx:
                n_ds = int(positions[0, i])
                ds = self._io._nsds.dynamicalSystem(n_ds)
                self._io._interman.removeBody(ds)
                self._io._nsds.removeDynamicalSystem(ds)
                print("remove ds number ", ds.number(), " with height = ", z[i])


dh = death_hook()

# Create solver options
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-3
sn.numerics_set_verbose(2)
with MechanicsHdf5Runner(mode="w", io_filename=fn) as io:
    ch = chute.create_chute(
        io,
        box_height=box_height,
        box_length=box_length,
        box_width=box_width,
        planethickness_=planethickness_,
        scale=1,
        trans=[-0.9, -1.8, -1],
    )

    rcs = rocas.create_rocas(
        io,
        n_layer=n_layer,
        n_row=n_row,
        n_col=n_col,
        x_shift=1.8,
        roca_size=0.1,
        top=3,
        rate=0.2,
        density=density,
        rock_shape=rock_shape,
    )

    io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=0.01)

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = hstep
run_options["theta"] = 1.0

bullet_options = SiconosBulletOptions()
bullet_options.perturbationIterations = 0
# bullet_options.minimumPointsPerturbationThreshold = 0.

run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options
run_options["constraint_activation_threshold"] = 1e-05

run_options["end_run_iteration_hook"] = dh


# run_options['Newton_options']=simu.LINEAR
run_options["Newton_options"] = simu.NONLINEAR
run_options["Newton_max_iter"] = 10
run_options["Newton_tolerance"] = 1e-8
run_options["display_Newton_convergence"] = True

run_options["skip_last_update_output"] = True
run_options["skip_reset_lambdas"] = True

run_options["osns_assembly_type"] = nsf.REDUCED_DIRECT

run_options["verbose"] = True
run_options["violation_verbose"] = True
run_options["with_timer"] = True
run_options["explode_computeOneStep_in_python"] = False
run_options["explode_computeOneStepNSProblem_in_python"] = False

# run_options['numerics_verbose']=True
# run_options['numerics_verbose_level']=0

run_options["output_frequency"] = 100

run_options["time_stepping"] = None


with MechanicsHdf5Runner(mode="r+", io_filename=fn) as io:
    io.run(run_options)
