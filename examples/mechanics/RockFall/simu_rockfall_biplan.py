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
#
import argparse
import random
from pathlib import Path

import siconos.integrators
import siconos.io.rocks_generator as rg
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)

# from convexhull_modif import ConvexHull
from siconos.mechanics.collision.bullet import SiconosBulletOptions

# siconos
from siconos.mechanics.collision.tools import Contactor

random.seed(42)
###############################################################################
# Definition
###############################################################################

# Simu parameters---------------------------------------------------------------
import os

# time
short_run = "SICONOS_EXAMPLES_SHORT_RUN" in os.environ
time = 1.0 if short_run else 20.0


# Soil param -------------------------------------------------------------------
# concat

zone_id = [1, 2]

nb_zones = len(zone_id)

# zone 1:
e1 = 0.0
mu1 = 0.8
mur1 = 0.4

# zone 2:
e2 = 0.0
mu2 = 0.7
mur2 = 0.1

e_c = [e1, e2]
mu_c = [mu1, mu2]
mu_r_c = [mu1, mur2]

mnt_file_prefix = "./data/mnt_nord_sud_id"
mnt_raster_file = "./data/mnt_nord_sud.asc"

# passage argument --output

parser = argparse.ArgumentParser(description="Simulation Rockfall")

parser.add_argument(
    "--nblocks",
    type=int,
    default=10,
    help="number of blocks for each run",
)

parser.add_argument(
    "--output",
    type=str,
    default=Path(__file__).with_suffix(".h5").name,
    help="Chemin complet du fichier HDF5 de sortie",
)
args = parser.parse_args()
if not short_run:
    time += args.nblocks

# time = 14.2
# output = f"time is {time}..."
# input(f"time is {time}...")


# -- Rock shape config --
rock_config = rg.RockShapeConfig(
    nb_pts=40,
    y_aspect_ratio=1.1,
    z_aspect_ratio=1.2,
    volume_min=1.0,
    volume_max=2.0,
    density=2650,
)

# -- Drop config --
raster_dep_file = "./data/zone_dep_mnt_nord_sud.stl"
drop_config = rg.RocksDropConfig(
    drop_zone=raster_dep_file,
    number_of_rocks=args.nblocks,
    height_fall_min=2.0,
    height_fall_max=3.0,
)


class deleterock:
    def __init__(self):
        pass

    def initialize(self, io):
        self._io = io

    def call(self, step):
        print("call death hook at step", step)
        # list contact points
        contact_points = self._io._io.contactPoints(self._io._nsds, 1)
        if contact_points is not None:
            for cp in contact_points:
                # print("cp", cp)
                inter_id = int(cp[22])
                print("inter_id", inter_id)
                ds1_id = int(cp[23])
                ds2_id = int(cp[24])
                print("ds id :", ds1_id, ds2_id)

                """
                if ds1_id == ds2_id:
                    print("interaction with a static object")
                    #inter = self._io._nsds.interaction(inter_id)
                    inter = siconos.modeling.interactions
                    contact_r = inter.relation()
                    print("contact_r.bodyShapeRecordA", contact_r.distance())
                    print("contact_r.bodyShapeRecordA", contact_r.bodyShapeRecordA)
                    print("contact_r.bodyShapeRecordB", contact_r.bodyShapeRecordB)

                    print("contact_r.bodyShapeRecordB.staticBody")
                    print(
                        "contact_r.bodyShapeRecordB.staticBody",
                        contact_r.bodyShapeRecordB.staticBody,
                    )
                    print(
                        "contact_r.bodyShapeRecordB.staticBody.number",
                        contact_r.bodyShapeRecordB.staticBody.number,
                    )
                    print(
                        "contact_r.bodyShapeRecordB.ds", contact_r.bodyShapeRecordB.ds
                    )

                    print("contact_r.bodyShapeRecordA.staticBody")
                    print(
                        "contact_r.bodyShapeRecordA.staticBody",
                        contact_r.bodyShapeRecordA.staticBody,
                    )
                    print(
                        "contact_r.bodyShapeRecordA.staticBody.ds.number()",
                        contact_r.bodyShapeRecordA.ds.number(),
                    )

                    if contact_r.bodyShapeRecordB.staticBody.number == static_body_id:
                        self._io._interman.removeStaticBody(
                            contact_r.bodyShapeRecordB.staticBody
                        )
                        # remove the body from that list of static object completely
                        for s in self._io._static:
                            # print(self._io._static[s])
                            # print(self._io._static[s]['number'])
                            if self._io._static[s]["number"] == static_body_id:
                                s_remove = s
                        self._io._static.pop(s_remove)
        else:
            # print('no contact points')
            pass

            # second way (faster) :  direct access to nsds positions
        # positions  = self._io._io.positions(self._io._nsds)
        # if positions is not None:
        #     z = positions[:,3]
        #     # We search for the ds index that are below a given criteria
        #     ds_idx = numpy.nonzero(z < -2)[0]
        #     for i in ds_idx :
        #         n_ds = int(positions[i,0])
        #         ds = self._io._nsds.dynamicalSystem(n_ds)
        #         self._io._interman.removeBody(ds)
        #         self._io._nsds.removeDynamicalSystem(ds)
        #         print('remove ds number ', ds.number(), ' with height = ', z[i])
        """


###############################################################################
# Simulation
###############################################################################

with MechanicsHdf5Runner(io_filename=args.output) as io:
    # soil ---------------------------------------------------------------------
    print("Resample Soil mask... \n")
    for i in range(nb_zones):
        dem_file = mnt_file_prefix + str(zone_id[i]) + ".stl"
        # if i==1:
        # dem_file='mnt_2m_PR_simple_'+str(i)+'_mod.stl'
        # else:
        # dem_file='mnt_2m_PR_simple_'+str(i)+'.stl'
        mesh1 = io.add_mesh_from_file(
            f"dem_TIN_shape{i}",
            dem_file,
            scale=1,
            insideMargin=0.0,
            outsideMargin=0.0,
        )
        io.add_object(
            f"dem_TIN{i}",
            [Contactor(f"dem_TIN_shape{i}", collision_group=i)],
            translation=[0.0, 0.0, 0.0],
        )

    # blocks
    rg.generate_random_blocks(io, drop_config, rock_config)

    # contact laws
    for i in range(nb_zones):
        io.add_Newton_impact_rolling_friction_nsl(
            f"contact_soil_{int(i)}",
            e=e_c[i - 1],
            mu=mu_c[i - 1],
            mu_r=mu_r_c[i - 1],
            collision_group1=100,
            collision_group2=i,
        )

# TEST OPTIONS DE BASE EXAMPLE TUTORIALS

np_timestep = 0.001


bullet_options = SiconosBulletOptions()
# bullet_options.worldScale = .1
# bullet_options.contactBreakingThreshold = 0.4


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-5


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = time
run_options["h"] = np_timestep

run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options


# run_options['skip_last_update_output']=True
# run_options['skip_reset_lambdas']=True
# run_options['osns_assembly_type']= siconos.nonsmooth_formulations.REDUCED_DIRECT

run_options["Newton_options"] = siconos.simulation.NONLINEAR

run_options["Newton_max_iter"] = 50
run_options["Newton_tolerance"] = 1e-8

run_options["verbose"] = True
run_options["with_timer"] = False
# run_options["explode_computeOneStep_in_python"] = True
run_options["explode_computeOneStepNSProblem_in_python"] = True

# run_options['violation_verbose'] = True
run_options["output_frequency"] = 1

run_options["time_stepping"] = None

dr = deleterock()

run_options["end_run_iteration_hook"] = dr

with MechanicsHdf5Runner(mode="r+", io_filename=args.output) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.

    # We run the simulation with relatively low precision because we
    # are interested mostly just in the location of contacts with the
    # height field, but we are not evaluating the performance of
    # individual contacts here.
    io.run(run_options)
