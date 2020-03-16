#!/usr/bin/env python

import os

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner



import siconos.numerics as sn
import siconos.kernel as sk

from siconos.io.FrictionContactTrace import FrictionContactTraceParams

# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.


test=True
if test:
    step = 125
    hstep = 1e-2
    itermax=100
    tolerance = 1e-03
    size_stack=2
else:
    step = 125
    hstep = 1e-2
    itermax=1000
    tolerance = 1e-12
    size_stack=5




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
                                  translation=[x,y,z],
                                  mass=1.0)
                    k += 1
            N = N - 1 if N > 1 else N
            M = M - 1 if M > 1 else M
            W = W - 1
            z += height + sep

    # A column
    make_stack(0, -10, 1, 1, size_stack)

    # A pyramid
    make_stack(0, 0, size_stack, size_stack, size_stack)

    # A wall
    make_stack(0, 10, 1, size_stack, size_stack)

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


    
solver = sn.SICONOS_GLOBAL_FRICTION_3D_ADMM

dump_probability = .02
theta = 0.50

if not os.path.exists('box_stacks'):
    os.mkdir('box_stacks')

fileName = "./Box_stacks/Box_Stacks"
title = "Box_stacks"
description = """
Box stacks with Bullet collision detection
Moreau TimeStepping: h={0}, theta = {1}
One Step non smooth problem: {2}, maxiter={3}, tol={4}
""".format(hstep, theta, sn.solver_options_id_to_name(solver),
           itermax,
           tolerance)

mathInfo = ""

friction_contact_trace_params = FrictionContactTraceParams(
    dump_itermax=20, dump_probability=None,
    fileName=fileName, title=title,
    description=description, mathInfo=mathInfo)

options = sk.solver_options_create(solver)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = itermax
options.dparam[sn.SICONOS_DPARAM_TOL] = tolerance
    
# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=step*hstep,
           h=hstep,
           theta=theta,
           Newton_max_iter=1,
           solver_options=options,
           output_frequency=1,
           osi=sk.MoreauJeanGOSI,
           friction_contact_trace_params=friction_contact_trace_params)
