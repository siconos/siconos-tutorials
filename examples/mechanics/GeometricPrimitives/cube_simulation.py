#!/usr/bin/env python
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

import os
import os.path
import runpy
from pathlib import Path

from siconos.io.mechanics_run import MechanicsHdf5Runner

if not os.path.isfile("cube_scene.hdf5"):
    cube_scene_path = Path(__file__).absolute().parent
    cube_scene = Path(cube_scene_path, "cube_scene.py")
    print("The scene cube_scene.hdf5 is not created")
    # exec(open("./cube_scene.py").read())
    runpy.run_path(str(cube_scene))  # creates cube_scene.hdf5

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode="r+", io_filename="cube_scene.hdf5") as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(output_frequency=100)
