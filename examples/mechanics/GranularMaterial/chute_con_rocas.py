#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

import chute
import rocas
import random

import sys

if (len(sys.argv) < 2):
    dist = 'uniform'
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]

if dist not in ['uniform', 'double', 'exp']:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 2.0:
    print("mu = [0.1 .. 2.0]")
    sys.exit(1)


# hdf5 file name
fn = 'chute_con_rocas-{0}-mu-{1}.hdf5'.format(dist, mu)

random.seed(0)

box_height = 3.683
box_length = 6.900
box_width = 3.430


density = 2500
plane_thickness = 0.2
cube_size = 0.1

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
    T = 20.
    hstep = 1e-4


# A hook to remove body that to far from the hopper at the end of the iteration
import numpy
class death_hook():
    def __init__(self):
        pass

    def initialize(self, io):
        self._io= io
        pass

    def call(self, step):
        #print('call death hook at step', step)

        # # First way (slow) : loop over the bodies
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
        positions  = self._io._io.positions(self._io._nsds)
        if positions is not None:
            z = positions[:,3]
            # We search for the ds index that are below a given criteria
            ds_idx = numpy.nonzero(z < -2)[0]
            for i in ds_idx :
                n_ds = int(positions[i,0])
                ds = self._io._nsds.dynamicalSystem(n_ds)
                self._io._interman.removeBody(ds)
                self._io._nsds.removeDynamicalSystem(ds)
                print('remove ds number ', ds.number(), ' with height = ', z[i])



dh = death_hook()

# Create solver options
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3

with MechanicsHdf5Runner(mode='w', io_filename=fn) as io:
    ch = chute.create_chute(io, box_height=box_height,
                            box_length=box_length,
                            box_width=box_width,
                            plane_thickness=plane_thickness,
                            scale=1, trans=[-0.9, -1.8, -1])

    rcs = rocas.create_rocas(io, n_layer=n_layer, n_row=n_row, n_col=n_col,
                             x_shift=1.8, roca_size=0.1, top=3,
                             rate=0.2, density=density)

    io.add_Newton_impact_friction_nsl('contact', mu=1.0, e=0.01)

with MechanicsHdf5Runner(mode='r+', io_filename=fn) as io:
    io.run(with_timer=False,
           explode_Newton_solve=False,
           t0=0,
           T=T,
           h=hstep,
           multipoints_iterations=True,
           theta=1.0,
           Newton_max_iter=1,
           solver_options=options,
           output_frequency=10,
           end_run_iteration_hook=dh)
