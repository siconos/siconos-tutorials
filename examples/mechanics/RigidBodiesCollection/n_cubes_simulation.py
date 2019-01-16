#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


import os


filename = 'n_cubes_scene.hdf5'
filename_py = 'n_cubes_scene.py'

if not os.path.isfile(filename):
    print('creation of the scene that is not existing')
    exec(compile(open(filename_py, "rb").read(), filename_py, 'exec'))
    filename = 'n_cubes_simulation.hdf5'

test=True
if test:
    nstep=100
    step=0.0005
else:
    nstep=20000
    step=0.0005

with MechanicsHdf5Runner(mode='r+',io_filename=filename) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=nstep*step,
           h=step,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=100)
