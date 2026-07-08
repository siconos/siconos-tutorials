# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
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

import siconos.numerics as sn
import siconos.nonsmooth_formulations
import siconos.mechanics.collision.tools
import siconos.io.mechanics_run
import siconos.mechanics.collision.bullet

# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.

# Creation of the hdf5 file for input/output
with siconos.io.mechanics_run.MechanicsHdf5Runner() as io:

    width, depth, height = 1, 1, 1
    io.add_primitive_shape("Box", "Box", [width, depth, height])

    k = 0
    sep = 0.01

    def make_stack(X, Y, N, M, W):
        global k
        z = height / 2.0
        while W > 0:
            for i in range(N):
                for j in range(M):
                    x = (i - N / 2.0) * (width + sep) + X
                    y = (j - M / 2.0) * (depth + sep) + Y
                    io.add_object(
                        "box%03d" % k,
                        [siconos.mechanics.collision.tools.Contactor("Box")],
                        translation=[x, y, z],
                        mass=1.0,
                    )
                    k += 1
            N = N - 1 if N > 1 else N
            M = M - 1 if M > 1 else M
            W = W - 1
            z += height + sep

    # A column
    make_stack(0, -10, 1, 1, 5)

    # A pyramid
    make_stack(0, 0, 5, 5, 5)

    # A wall
    make_stack(0, 10, 1, 5, 5)

    # Definition of the ground
    io.add_primitive_shape("Ground", "Box", (50, 50, 0.1))
    io.add_object(
        "ground", [siconos.mechanics.collision.tools.Contactor("Ground")], [0, 0, -0.05]
    )

    # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[30,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Fremond_impact_friction_nsl('contact', mu=0.3)


bullet_options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.4
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

run_options = siconos.io.mechanics_run.MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 1.0
run_options["h"] = 0.01
# run_options['theta'] = 1.

run_options["Newton_max_iter"] = 1
run_options["Newton_tolerance"] = 1e-6

run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options


# run_options['skip_last_update_output'] = True
# run_options['skip_reset_lambdas'] = True
run_options["osns_assembly_type"] = siconos.nonsmooth_formulations.REDUCED_DIRECT

run_options["verbose"] = True
run_options["with_timer"] = True
run_options['explode_computeOneStep_in_python'] = True
# run_options['explode_computeOneStepNSProblem_in_python']=True

# run_options['violation_verbose'] = True
run_options["output_frequency"] = 1

test = True
if test:
    run_options["T"] = 0.2
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-06
else:
    run_options["T"] = 10.0

# Load and run the simulation
with siconos.io.mechanics_run.MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
