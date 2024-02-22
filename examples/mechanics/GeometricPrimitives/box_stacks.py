#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import siconos.numerics as sn
import siconos.kernel as sk
# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    width, depth, height = 1, 1, 1
    io.add_primitive_shape('Box', 'Box', [width, depth, height])

    k = 0
    sep = 0.01

    def make_stack(X, Y, N, M, W):
        global k
        z = height/2.0
        while W > 0:
            for i in range(N):
                for j in range(M):
                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-M/2.0)*(depth+sep) + Y
                    io.add_object('box%03d' % k, [Contactor('Box')],
                                  translation=[x, y, z],
                                  mass=1.0)
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
    io.add_primitive_shape('Ground', 'Box', (50, 50, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0, 0, -0.05])

    # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[30,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)


bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.
bullet_options.contactBreakingThreshold = 0.4
bullet_options.perturbationIterations = 3.
bullet_options.minimumPointsPerturbationThreshold = 3.

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10

run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=1.0
run_options['h']=0.01
#run_options['theta']=1.

run_options['Newton_max_iter'] =1
run_options['Newton_tolerance'] =1e-6

run_options['bullet_options']=bullet_options
run_options['solver_options']=options


#run_options['skip_last_update_output']=True
#run_options['skip_reset_lambdas']=True
run_options['osns_assembly_type']= sk.REDUCED_DIRECT

run_options['verbose']=True
run_options['with_timer']=True
#run_options['explode_Newton_solve']=True
#run_options['explode_computeOneStep']=True

#run_options['violation_verbose'] = True
run_options['output_frequency']=1

test=True
if test:
    run_options['T']=1.0
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-06
else:
    run_options['T']=10.0

# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(run_options)
