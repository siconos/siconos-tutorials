#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk
# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    width, depth, height = 1, 1, 1
    io.add_primitive_shape('Box', 'Box', [width, depth, height])
    io.add_primitive_shape('Sphere', 'Sphere', (width/2.,))
    io.add_primitive_shape('Capsule', 'Capsule', (width/4., depth/4.))
    io.add_primitive_shape('Cone', 'Cone', (width/2., depth/2.))
    io.add_primitive_shape('Cylinder', 'Cylinder', (width/2., depth/2.))

    shape_list = ['Box', 'Sphere', 'Capsule', 'Cone', 'Cylinder']

    k = 0
    sep = 0.01

    def make_stack(X, Y, N, M, W):
        global k
        z = height/2.0
        for l in range(W):
        #while W >0 :
            for i in range(N):
                for j in range(M):
                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-M/2.0)*(depth+sep) + Y
                    
                    
                    # select an object in the shape list

                    obj_idx = (i+2*k+17*l) % 5
                    #print('obj-idx', obj_idx)
                    
                    io.add_object('object%03d' % k, [Contactor(shape_list[obj_idx])],
                                  translation=[x, y, z],
                                  mass=1.0)
                    k += 1
            # N = N - 1 if N > 1 else N
            # M = M - 1 if M > 1 else M
            # W = W - 1
            z += height + sep

    make_stack(0, 0, 4, 4, 4)


    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (50, 50, 0.5))
    
    io.add_object('ground',  [Contactor('Ground')], translation= [0, 0, -0.5])

    # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[30,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 1.0
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3

test=True
if test:
    T=1.0
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
else:
    T=10.0
# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=T,
           h=0.001,
           theta=0.5,
           Newton_max_iter=1,
           solver_options=options,
           bullet_options=bullet_options,
           output_frequency=1)

