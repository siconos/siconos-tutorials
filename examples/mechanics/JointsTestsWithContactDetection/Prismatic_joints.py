#!/usr/bin/env python
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

with MechanicsHdf5Runner() as io:
    io.add_primitive_shape('Cube', 'Box', (1, 1, 1))
    io.add_primitive_shape('Ground', 'Box', (10, 10, .5))
    io.add_Newton_impact_friction_nsl('contact', e=0.0, mu=0.2)
    io.add_object('ground', [Contactor('Ground')], translation=[0, 0, -0.25])

    io.add_object('cube1', [Contactor('Cube')], translation=[-1, 0, 0.5], mass=0.1)
    io.add_object('cube2', [Contactor('Cube')], translation=[ 1, 0, 0.5], mass=0.1,
                  velocity=[5,0,5,0,0,0])

    io.add_joint('joint1', 'cube1', 'cube2', None, [[0,0,1]], 'PrismaticJointR')


    
test =True
if test:
    T=0.1
else:
    T=3.
    
options = sk.solver_options_create(sn.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4

with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=T,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           solver_options=options,
           )
