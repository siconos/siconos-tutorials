#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk


with MechanicsHdf5Runner() as io:
    io.add_primitive_shape('BigBox', 'Box', (1, 1, 1))
    io.add_primitive_shape('SmallBox', 'Box', (0.3, 0.3, 0.3))
    io.add_Newton_impact_friction_nsl('contact', e=0.0, mu=0.0)

    io.add_object('heavy1', [Contactor('SmallBox')], translation=[0, 0, 4],
                  mass=1.0,
                  velocity=[0, 0, 15 ,0 ,0 ,0 ])
    io.add_object('heavy2', [Contactor('SmallBox')], translation=[0, 0, -4],
                  mass=1.0,
                  velocity=[0, 0, 15 ,0 ,0 ,0 ])

    io.add_object('big', [Contactor('BigBox')], translation=[0, 0, 0],
                  mass=0.1,
                  velocity=[0,0,0,0,0,0])

    # We create a prismatic joint with friction on its axis.  In fact
    # as there is no "normal force" in a 1-D friction, it is expressed
    # as a lower- and upper-bound of permissible velocity by means of
    # a RelayNSL, here between [-0.1, 0.1].

    io.add_relay_nsl('friction', lb=-0.1, ub=0.1)

    io.add_joint('joint1', 'big', None, None, [[0,0,1]], 'PrismaticJointR',
                 friction='friction')

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4
 
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=8,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=2,
           solver_options=options,
           )
