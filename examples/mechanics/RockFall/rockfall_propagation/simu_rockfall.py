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
#

# passage argument --output
import argparse
import shutil
from pathlib import Path

import siconos.integrators
import siconos.io.rocks_generator as rg
import siconos.numerics as sn
import siconos.simulation as simu
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)

# from convexhull_modif import ConvexHull
from siconos.mechanics.collision.bullet import SiconosBulletOptions

# siconos
from siconos.mechanics.collision.tools import Contactor

###############################################################################
# Definition
###############################################################################

# Simu parameters---------------------------------------------------------------

# time
time = 35.0  # 35.0  # sec  ok avec 1 dt

# Warning arg parsing will work only when the script is called with python
# not with 'siconos'.
parser = argparse.ArgumentParser(description="Simulation Rockfall")

parser.add_argument("--nblocks", type=int, default=100, help="nombre blocs par run")

parser.add_argument(
    "--output",
    type=str,
    default=Path(__file__).with_suffix(".h5").name,
    help="Chemin complet du fichier HDF5 de sortie",
)
args = parser.parse_args()


do_build_config = False

if do_build_config:
    # Soil param -------------------------------------------------------------------
    # concat
    zone_id = [1, 2, 12, 13, 15, 16, 17, 18, 19]
    zone_id = [1, 2, 12, 13, 15, 16, 17, 18]

    # zone 1 == none !!
    nb_zones = 9
    nb_zones = 8

    # zone 13:
    e13 = 0.0
    mu13 = 0.6
    mur13 = 0.35
    # mur13 = 0.1
    # zone 1:
    e1 = 0.0
    mu1 = 0.7
    mur1 = 0.39
    # mur1 = 0.15
    # zone 2:
    e2 = e1
    mu2 = mu1
    mur2 = mur1
    # zone 12:
    e12 = 0.0
    mu12 = 0.7
    mur12 = 0.42
    # mur12 = 0.2
    # zone 15:
    e15 = e12
    mu15 = mu12
    mur15 = mur12
    # zone 16:
    e16 = 0.0
    mu16 = 0.8
    mur16 = 0.55
    # mur16 = 0.25
    # zone 17:
    e17 = e1
    mu17 = mu1
    mur17 = mur1

    # zone 18:
    e18 = e17
    mu18 = mu17
    mur18 = mur17
    # zone 19:
    e19 = e18
    mu19 = mu18
    mur19 = mur18

    e_c = [e1, e2, e12, e13, e15, e16, e17, e18, e19]
    mu_c = [mu1, mu2, mu12, mu13, mu15, mu16, mu17, mu18, mu19]
    mu_r_c = [mu1, mur2, mur12, mur13, mur15, mur16, mur17, mur18, mur19]

    mnt_file_prefix = "./data/dem_red_site_reduit_simple_id"

    mnt_raster_file = "./data/dem_red_site_reduit.asc"

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
    # raster_dep_file = "./data/zones_dep_red_site_reduit_simplifie.stl"
    raster_dep_file = "./data/dem_red_site_reduit_simple_zone_dep_calib.stl"
    drop_config = rg.RocksDropConfig(
        drop_zone=raster_dep_file,
        number_of_rocks=args.nblocks,
        height_fall_min=2.0,
        height_fall_max=3.0,
    )
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

    # Save initial configuration (repr. purpose)
    config_backup = Path(args.output).stem + "_start.h5"
    shutil.copy(args.output, config_backup)

    run_config = args.output
else:
    # Copy backup file to restart simu from an existing setup
    config_backup = Path(args.output).stem + "_start.h5"
    run_config = args.output
    shutil.copy(config_backup, run_config)


###############################################################################
# Simulation
###############################################################################

# LIST OF NUMERICAL PARAMETERS

# Siconos options
np_timestep = 1e-3  # timetstep
np_theta = 0.5  # for the timestep scheme (theta=1 : implicit ; theta=0 : explicit)
np_gamma = 0.0  # velocity prediction : $q_k=h\gamma+v_k$ #on s'en bat les couilles
np_Newton_max_iter = 20  # iterations of the newton algorithm : 20
# (to be set at 2 or 3 -> reach tolerance 10-5 in 2 iterations)
np_Newton_tolerance = 1e-8  # tolerance ofthe newton algorithm (default : 1e-10)
# but we decrease it for these applications
np_constraint_activation_threshold = (
    0.0  # epsilon with $q_{k+1}-q_kh\gamma-v_k<\epsilon$
)
np_freq_output = 1  # frequency of output data (number of timestep)

# bullet_options
np_bull_contactBreakingThreshold = (
    0.2  # distance at whcih contacts are deactivated (and activated ?)
)
# solver_options (solver dependant)
np_solv_iter_max = 10000
np_solv_tol = 1e-6
# np_solv_freezing=20

# Siconos options
"""        d['h']=0.0005
        d['multipoints_iterations']=None # deprecated
        d['theta']=0.50001
        d['gamma']=0.0
        d['Newton_options']=sk.SICONOS_TS_NONLINEAR
        d['Newton_max_iter']=20
        d['Newton_tolerance']=1e-10
        d['Newton_warning_on_nonconvergence']=True
        d['Warning_nonsmooth_solver']=True
        d['set_external_forces']=None
        d['solver_options']=None
        d['solver_options_pos']=None
        d['osnspb_max_size']=0
        d['exit_tolerance']=None
        d['projection_itermax']=20
        d['projection_tolerance']=1e-8
        d['projection_tolerance_unilateral']=1e-8
        d['numerics_verbose']=False
        d['numerics_verbose_level']=0
        d['violation_verbose']=False
        d['verbose']=True
        d['verbose_progress']=True
        d['output_frequency']=None
        d['output_backup']=False
        d['output_backup_frequency']=None
        d['friction_contact_trace_params']=None
        d['output_contact_index_set']=1
        d['osi']=sk.MoreauJeanOSI
        d['constraint_activation_threshold']=0.0
        d['explode_Newton_solve']=False
        d['explode_computeOneStep']=False
        d['display_Newton_convergence']=False
        d['start_run_iteration_hook']=None
        d['before_next_step_iteration_hook']=None
        d['end_run_iteration_hook']=None
        d['skip_last_update_output']=False
        d['skip_last_update_input']=False
        d['skip_reset_lambdas']=False
        d['osns_assembly_type']= None
        d['output_contact_forces']=True,
        d['output_contact_info']=True,
        d['output_contact_work']=True,
"""

# bullet_options
"""
siconosBulletOptions::SiconosBulletOptions()
  : dimension(SICONOS_BULLET_3D)
  , contactBreakingThreshold(0.02) "distance à laquelle le contact est supprimé
  ,(et activé je pense ... jes suis sur)"
  , contactProcessingThreshold(0.03) " not used"
  , worldScale(1.0) ( a mette grand quand les bodies sont petits
  , pour se ramerner à échelle métrique)
  , useAxisSweep3(false) : je sais pas
  , clearOverlappingPairCache(false) : si true refait tt
  , la detection de contact à chaque pas de temsp
  , perturbationIterations(3) : à mettre à 0 pour avoir un pt de contact
  , minimumPointsPerturbationThreshold(3) : je sais pas - A VOIR
  , enableSatConvex(false) - pas utilisé (autre type algo de détection de bullet)
  , enablePolyhedralContactClipping(false) - pas utilisé
  , (autre type algo de détection de bullet)
  , Depth2D(0.04) - profondeur car solveur 3D de bullet utilsié pour le 2D
"""

# solver_options (solver dependant)
"""
#solver_options (solver dependant)
# NSGS:
#  https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/doxygen/Friction__cst_8h.html
    - iter_max : (nbre d iteration du NSGS) : 1000 par défaut
     (généralement atteint la tol bien avant) #tous les solveurs le partagent
    - tol : 10-4 (suffisante sauf pour les forces) #tous les solveurs le partagent
    - freezing : frequence de pas de Gauss Seidel à laquelle les contacts
     dont les valeurs de p sont stabilisés
      (en plus ca accélère la convergence du Gauss Seidel)
# Create solver options
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10
"""
# END LIST OF NUMERICAL PARAMETERS

# setting of the numerical parameters#bullet

bullet_options = SiconosBulletOptions()
bullet_options.contactBreakingThreshold = np_bull_contactBreakingThreshold

# solver
options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = np_solv_iter_max
options.dparam[sn.params.SICONOS_DPARAM_TOL] = np_solv_tol
# options.iparam[sn.params.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = np_solv_freezing;


class deleterock:
    def __init__(self):
        pass

    def initialize(self, io):
        self._io = io

    def call(self, step):
        # print("call death hook at step", step)
        positions = self._io.get_io_array(self._io._io.positions(self._io._nsds))
        # velocities = self._io._io.velocities(self._io._nsds)
        nsds = self._io._nsds
        min_kinetic = 1e-1
        if positions is not None:
            ds_idx = positions[:, 0]
            # print('ds_idx',ds_idx)
            for i in ds_idx:
                n_ds = int(i)
                neds = nsds.dynamicalSystem(n_ds)
                kinetic = neds.computeKineticEnergy()
                # print('kinetic',kinetic,i)
                if kinetic < min_kinetic:
                    self._io._interman.removeBody(neds)
                    self._io._nsds.removeDynamicalSystem(neds)
                    print(
                        "remove ds number ",
                        neds.number(),
                        " with kinetic enrgy = ",
                        kinetic,
                    )
                    # input()


# run
run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = time
run_options["h"] = np_timestep
run_options["theta"] = np_theta
run_options["gamma"] = np_gamma
run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options
run_options["Newton_options"] = simu.NONLINEAR  # WARNING : BEFORE IT WAS "LINEAR"
run_options["Newton_max_iter"] = np_Newton_max_iter
run_options["Newton_tolerance"] = np_Newton_tolerance
run_options["Newton_warning_on_nonconvergence"] = True
run_options["Warning_nonsmooth_solver"] = True
# run_options['osnspb_max_size']=0 #### WARNING : ASK VINCENT
# run_options['exit_tolerance']=None #### WARNING : ASK VINCENT
run_options["skip_last_update_output"] = True  # A REVOIR#######
run_options["skip_reset_lambdas"] = True  # A REVOIR#######
# run_options['osns_assembly_type']= nsf.REDUCED_DIRECT #####A REVOIR#######

run_options["numerics_verbose"] = True
run_options["numerics_verbose_level"] = 0
run_options["violation_verbose"] = True
run_options["verbose"] = True
run_options["verbose_progress"] = False
run_options["explode_computeOneStep_in_python"] = False
run_options["explode_computeOneStepNSProblem_in_python"] = False
run_options["display_Newton_convergence"] = True

run_options["output_frequency"] = np_freq_output

run_options["output_contact_index_set"] = 1  # ???
run_options["osi"] = siconos.integrators.MoreauJeanOSI

run_options["constraint_activation_threshold"] = np_constraint_activation_threshold

run_options["output_contact_forces"] = True
run_options["output_contact_info"] = True
run_options["output_contact_work"] = False

dr = deleterock()
run_options["end_run_iteration_hook"] = dr

# unused options
"""
run_options['start_run_iteration_hook']=None
run_options['before_next_step_iteration_hook']=None
run_options['end_run_iteration_hook']=None
run_options['skip_last_update_output']=False
run_options['skip_last_update_input']=False
run_options['skip_reset_lambdas']=False
run_options['osns_assembly_type']= None
"""

# input()

# solve

with MechanicsHdf5Runner(mode="r+", io_filename=run_config) as io:
    io.run(run_options)
