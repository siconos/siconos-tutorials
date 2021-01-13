
from __future__ import print_function
import os,sys
import numpy as np
import math
import siconos.kernel as Kernel

np.set_printoptions(precision=3)

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.joints import cast_PrismaticJointR
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.kernel import SiconosVector, BlockVector

import siconos.numerics as sn
import siconos.kernel as sk

# An example of applying force to the axis of a joint, and applying
# spring and virtual damping by measuring position and velocity along
# the same axis.

# Note: This example is to demonstrate external measurement of joint
# positions and application of forces to dynamical systems attached to
# joints.  In practice it is better to use internal forces (fInt,
# mInt) to model joint spring-dampers, see folder
# JointsTestsWithInternalForces, and extra Relations with associated
# Non-Smooth Laws to model non-linearities such as joint stops and
# friction, see JointsTestsWithContactDetection.

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of two bars connected by a prismatic joint
    io.add_primitive_shape('Bar', 'Box', (1, 0.1, 0.1))
    io.add_object('bar1', [Contactor('Bar')], [0.05,0,2],
                 orientation=[(0,0,1),np.pi/2], mass=1.0, velocity=[0,0,0,0,0,1])
    io.add_object('bar2', [Contactor('Bar')], [-0.05,0,2],
                 orientation=[(0,0,1),np.pi/2], mass=1.0)
    io.add_joint('joint1', 'bar1', 'bar2', None, [[0,1,0]], 'PrismaticJointR', True)

    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (5, 5, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0,0,-0.05])

class Ctrl(object):
    def initialize(self, io):
        self.count = 0
        self.topo = io._nsds.topology()
        self.ds1 = Kernel.cast_NewtonEulerDS(self.topo.getDynamicalSystem('bar1'))
        self.ds2 = Kernel.cast_NewtonEulerDS(self.topo.getDynamicalSystem('bar2'))
        
        self.joint1 = cast_PrismaticJointR(
            self.topo.getInteraction('joint1').relation())

        # Apply initial forces
        self.step()

    def step(self):
        self.count += 1

        # Make a temporary BlockVector containing both qs
        bv = BlockVector(self.ds1.q(), self.ds2.q())

        force1 = np.zeros(3)
        force2 = np.zeros(3)

        # Apply an impulse to watch the response
        if (self.count == 1000):
            tmp = np.array(self.joint1.normalDoF(bv, 0)) * -1000
            print('applying impulse', tmp)
            force1 +=  np.array(tmp)/2
            force2 += -np.array(tmp)/2

        # Get the position and use it to project a force vector
        # onto the DoF (spring force)
        pos = SiconosVector(1)
        self.joint1.computehDoF(0, bv, pos, 0)

        setpoint = 1.0
        pos_diff = setpoint - pos.getValue(0)
        spring_force = np.array(self.joint1.normalDoF(bv, 0)) * pos_diff * 100.0

        # Get the velocity of each body projected onto the DoF and
        # calculate their difference (damping force)
        vel1 = self.joint1.projectVectorDoF(self.ds1.linearVelocity(True), bv, 0)
        vel2 = self.joint1.projectVectorDoF(self.ds2.linearVelocity(True), bv, 0)
        vel_diff = vel1 - vel2
        damping_force = vel_diff * 10.0

        # Calculate total forces for each body
        force1 += - (spring_force + damping_force)/2
        force2 += + (spring_force + damping_force)/2

        print('applying spring-damper forces', force1, force2)

        self.ds1.setFExtPtr(force1)
        self.ds2.setFExtPtr(force2)

options = sk.solver_options_create(sn.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-12
        
# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=20,
           h=0.01,
           theta=0.5,
           Newton_max_iter=1,
           controller=Ctrl(),
           solver_options=options,
           output_frequency=1)
