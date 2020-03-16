#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner

import os.path
import os
from pathlib import Path

if not os.path.isfile('cube_scene.hdf5'):
    cube_scene_path = Path(__file__).absolute().parent
    cube_scene = Path(cube_scene_path, 'cube_scene.py')
    print('The scene cube_scene.hdf5 is not created')
    #exec(open("./cube_scene.py").read())
    exec(open(cube_scene).read())

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+', io_filename='cube_scene.hdf5') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(output_frequency=100)
