# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
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
#
#
# Importing the necessary packages from python
import math
import os
import time

import numpy
import numpy as np
import siconos
import siconos.numerics as Numerics
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.mechanics.collision.tools import Contactor

# from numpy import seed
from stl import mesh

siconos.integrators.enable_solver_check()

# List of input paramters #######
# The length of the structure
L_structure = 14.0855

# ETAG projetile paramteres
m_projectile = 2600.0
impact_energy = 10.2e5
angle_impact = -0.0 * 3.14 / 180
L_b = 1.1

# The on-site experimental loacation
X_block = -1.05 * (0.5 * L_b)
Y_block = 0.5 * L_structure
Z_block = 1.7

# How much disk offset do you want?
disk_offset = 0.068
# The vertical play in total
VP_total = 0.071
# Definition of a non smooth law. friction
mu_t = 0.316
# Defining a different friction between concrete and soil
mu_t_cs = 0.307
# Definition of non-smooth law restitution coefficient for all contacts
e_t = 0.222
# Defining the contacts for steel-steel and steel-others as follows
mu_steel = 0.2
e_steel = 0
# This is just for the display: use 1 for the true model
display_f = 1.0

# Importing the block geomerty #######
# Defining the wall configuration
Nb_block_vert = 4.0
Nb_pattern = 19.0
#
# Pattern definition
Bottom_block = False
rot_block = math.pi / 4  # 0.785 #0.804
Orientation_list = [
    0,
    -rot_block,
    -rot_block,
    0,
    0,
    0,
    rot_block,
    rot_block,
    0,
    0,
    0,
    -rot_block,
    -rot_block,
    0,
    0,
    0,
    rot_block,
    rot_block,
    0,
]
#
# concrete block size
block_width = 1.56
block_depth = 0.76
block_height = 0.8

# block hole
block_hole_y1 = 0.38 - 0.5 * block_width
block_hole_y2 = (0.8 + 0.38) - 0.5 * block_width
disk_position_z1 = disk_offset - 0.5 * block_height
disk_position_z2 = 0.5 * block_height - disk_offset
disk_position_z3 = -disk_offset
disk_position_z4 = disk_offset
block_hole_diameter = 0.154
#
# block mass
mass_block = 1800.0
mass_half_block = (
    mass_block / 2
)  # Ritesh: This is for the four blocks present in the exterior!
#
# steel bar
bar_diam_ext = 0.1397
bar_diam_int = 0.1237
#
#
# Spherical vertical play
factor_bar_diam_ext = 0.5
the_vertical_play = VP_total / 6  # The total vertical play is divided into six segments
#
# How much vertical play per block do you want?
VP_box_depth = factor_bar_diam_ext * bar_diam_ext
#
# for box-box
VP_offset = 0.5 * the_vertical_play + VP_box_depth / 2
# for sphere-box
VP_offset_sphere = 0.5 * the_vertical_play + VP_box_depth
cylinder_radius = bar_diam_ext / 2.0
bar_length = (
    Nb_block_vert * block_height
)  # Ritesh: This represents the steel bar length equal to the structure height
mass_bar = 7800 * math.pi * (bar_diam_ext**2 - bar_diam_int**2) / 4 * bar_length
#
# Importing the function to create the pattern ##############


def make_pattern(
    io,
    X,
    Y,
    orientation,
    Nb_block_vert,
    mass_block,
    Bottom_block,
    num_block,
    inertia_block,
):
    if Nb_block_vert % 2 == 0:
        Nb_block_vert_pattern = int(Nb_block_vert / 2)
    else:
        if Bottom_block is True:
            Nb_block_vert_pattern = int(Nb_block_vert / 2) + 1
        else:
            Nb_block_vert_pattern = int(Nb_block_vert / 2)
    if Bottom_block is True:
        z_block = 0.5 * block_height
    else:
        z_block = 1.5 * block_height

    row_number = 1
    while row_number <= Nb_block_vert_pattern:
        x_block = X + (0.5 * block_depth + 0.02) * math.sin(orientation)
        y_block = Y + (0.5 * block_depth + 0.02) * math.cos(orientation)

        if z_block == 0.5 * block_height:
            create_block_bottom(
                io,
                x_block,
                y_block,
                z_block,
                num_block,
                orientation,
                mass_block,
                inertia_block,
            )
        elif z_block == 3.5 * block_height:
            create_block_top(
                io,
                x_block,
                y_block,
                z_block,
                num_block,
                orientation,
                mass_block,
                inertia_block,
            )
        else:
            create_block_middle(
                io,
                x_block,
                y_block,
                z_block,
                num_block,
                orientation,
                mass_block,
                inertia_block,
            )

        z_block = z_block + 2 * block_height
        row_number = row_number + 1
        num_block = num_block + 1
    if Bottom_block is True:
        Bottom_block = False
    else:
        Bottom_block = True
    return [num_block, Bottom_block]


def create_block_bottom(io, X, Y, Z, k, orientation, mass_block, inertia_block):
    #
    # create object
    io.add_object(
        f"block{k:03d}",
        [
            Contactor("Box_B", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, 0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "Block_curvature",
                relative_translation=[0, -0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[
                    0.0,
                    block_hole_y1,
                    -(block_height / 2 - VP_offset),
                ],
                collision_group=6,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[
                    0.0,
                    block_hole_y2,
                    -(block_height / 2 - VP_offset),
                ],
                collision_group=6,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[math.cos(orientation / 2), 0, 0, -math.sin(orientation / 2)],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_block,
        inertia=inertia_block,
    )


def create_block_middle(io, X, Y, Z, k, orientation, mass_block, inertia_block):
    #
    # create object
    io.add_object(
        f"block{k:03d}",
        [
            Contactor("Box_B", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, 0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "Block_curvature",
                relative_translation=[0, -0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z4],
                collision_group=5,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[math.cos(orientation / 2), 0, 0, -math.sin(orientation / 2)],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_block,
        inertia=inertia_block,
    )


def create_block_top(io, X, Y, Z, k, orientation, mass_block, inertia_block):
    #
    # create object
    io.add_object(
        f"block{k:03d}",
        [
            Contactor("Box_B", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, 0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "Block_curvature",
                relative_translation=[0, -0.4, 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y1, disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, block_hole_y2, disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    block_hole_y1,
                    (block_height / 2 - VP_offset_sphere),
                ],
                collision_group=7,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    block_hole_y2,
                    (block_height / 2 - VP_offset_sphere),
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[math.cos(orientation / 2), 0, 0, -math.sin(orientation / 2)],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_block,
        inertia=inertia_block,
    )


#
# Importing the function to create the bars #######
len_cyl1 = block_height
mass_cyl1 = len_cyl1 * mass_bar / bar_length
inertia_cyl1 = numpy.eye(3)
inertia_cyl1[0, 0] = (
    1
    / 12.0
    * mass_cyl1
    * (3 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2) + len_cyl1**2)
)
inertia_cyl1[1, 1] = (
    1 / 2.0 * mass_cyl1 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2)
)
inertia_cyl1[2, 2] = (
    1
    / 12.0
    * mass_cyl1
    * (3 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2) + len_cyl1**2)
)
len_cyl2 = 0.5 * block_height
mass_cyl2 = len_cyl2 * mass_bar / bar_length
inertia_cyl2 = numpy.eye(3)
inertia_cyl2[0, 0] = (
    1
    / 12.0
    * mass_cyl2
    * (3 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2) + len_cyl2**2)
)
inertia_cyl2[1, 1] = (
    1 / 2.0 * mass_cyl2 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2)
)
inertia_cyl2[2, 2] = (
    1
    / 12.0
    * mass_cyl2
    * (3 * ((0.5 * bar_diam_ext) ** 2 + (0.5 * bar_diam_int) ** 2) + len_cyl2**2)
)


def make_bars(io, X, Y, num_pattern):
    # create_shape
    k = num_pattern + 1
    io.add_object(
        f"bar1{k:03d}",
        [
            Contactor("Cyl2", collision_group=2),
            Contactor(
                "Box_Stop",
                relative_translation=[
                    0.0,
                    0.5
                    * (
                        bar_length / 16.0
                        - (VP_offset + 2 * (factor_bar_diam_ext * bar_diam_ext) + 0.002)
                    ),
                    0.0,
                ],
                relative_orientation=[1, 1, 0, 0],
                collision_group=4,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, 0.5 * (bar_length / 16.0 + VP_offset), 0.0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=6,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    -0.5 * (bar_length / 16.0 + VP_offset_sphere),
                    0.0,
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, bar_length / 16.0],
        orientation=[1, 1, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_cyl2,
        inertia=inertia_cyl2,
    )

    io.add_object(
        f"bar2{k:03d}",
        [
            Contactor("Cyl1", collision_group=3),
            Contactor(
                "Box_Stop",
                relative_translation=[
                    0.0,
                    0.5
                    * (
                        bar_length / 8.0
                        - (VP_offset + 2 * (factor_bar_diam_ext * bar_diam_ext) + 0.002)
                    ),
                    0.0,
                ],
                relative_orientation=[1, 1, 0, 0],
                collision_group=4,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, 0.5 * (bar_length / 8.0 + VP_offset), 0.0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=6,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    -0.5 * (bar_length / 8.0 + VP_offset_sphere),
                    0.0,
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, bar_length / 4.0],
        orientation=[1, 1, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_cyl1,
        inertia=inertia_cyl1,
    )

    io.add_object(
        f"bar3{k:03d}",
        [
            Contactor("Cyl1", collision_group=2),
            Contactor(
                "Box_Stop",
                relative_translation=[
                    0.0,
                    0.5
                    * (
                        bar_length / 8.0
                        - (VP_offset + 2 * (factor_bar_diam_ext * bar_diam_ext) + 0.002)
                    ),
                    0.0,
                ],
                relative_orientation=[1, 1, 0, 0],
                collision_group=4,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, 0.5 * (bar_length / 8.0 + VP_offset), 0.0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=6,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    -0.5 * (bar_length / 8.0 + VP_offset_sphere),
                    0.0,
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, bar_length / 2.0],
        orientation=[1, 1, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_cyl1,
        inertia=inertia_cyl1,
    )

    io.add_object(
        f"bar4{k:03d}",
        [
            Contactor("Cyl1", collision_group=3),
            Contactor(
                "Box_Stop",
                relative_translation=[
                    0.0,
                    0.5
                    * (
                        bar_length / 8.0
                        - (VP_offset + 2 * (factor_bar_diam_ext * bar_diam_ext) + 0.002)
                    ),
                    0.0,
                ],
                relative_orientation=[1, 1, 0, 0],
                collision_group=4,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, 0.5 * (bar_length / 8.0 + VP_offset), 0.0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=6,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    -0.5 * (bar_length / 8.0 + VP_offset_sphere),
                    0.0,
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, 3.0 * bar_length / 4.0],
        orientation=[1, 1, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_cyl1,
        inertia=inertia_cyl1,
    )

    io.add_object(
        f"bar5{k:03d}",
        [
            Contactor("Cyl2", collision_group=2),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, 0.5 * (bar_length / 16.0 + VP_offset), 0.0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=6,
            ),
            Contactor(
                "Sp_VP",
                relative_translation=[
                    0.0,
                    -0.5 * (bar_length / 16.0 + VP_offset_sphere),
                    0.0,
                ],
                collision_group=7,
            ),
        ],
        translation=[X, Y, bar_length - 0.2],
        orientation=[1, 1, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=0.1 * mass_cyl2,
        inertia=inertia_cyl2,
    )


#
# Importing the function to create a half block
def create_half_block_left(io, X, Y, Z, k, mass_half_block, inertia_half_block):
    # create object
    io.add_object(
        f"half_block{k:03d}",
        [
            Contactor("Box_half", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, -(0.2), 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z4],
                collision_group=5,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[0, 0, 0, 1],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_half_block,
        inertia=inertia_half_block,
    )


def create_half_block_left_bottom(io, X, Y, Z, k, mass_half_block, inertia_half_block):
    # create object
    io.add_object(
        f"half_block{k:03d}",
        [
            Contactor("Box_half", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, -(0.2), 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, -(0.2), disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, -(0.2), -(block_height / 2 - VP_offset)],
                collision_group=6,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[0, 0, 0, 1],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_half_block,
        inertia=inertia_half_block,
    )


def create_half_block_right(io, X, Y, Z, k, mass_half_block, inertia_half_block):
    # create object
    io.add_object(
        f"half_block{k:03d}",
        [
            Contactor("Box_half", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, (0.2), 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z4],
                collision_group=5,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[0, 0, 0, 1],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_half_block,
        inertia=inertia_half_block,
    )


def create_half_block_right_bottom(io, X, Y, Z, k, mass_half_block, inertia_half_block):
    # create object
    io.add_object(
        f"half_block{k:03d}",
        [
            Contactor("Box_half", collision_group=0),
            Contactor(
                "Block_curvature",
                relative_translation=[0, (0.2), 0],
                relative_orientation=[1, 1, 0, 0],
                collision_group=1,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z1],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z2],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z3],
                collision_group=5,
            ),
            Contactor(
                "disk",
                relative_translation=[0.0, (0.2), disk_position_z4],
                collision_group=5,
            ),
            Contactor(
                "Box_VP",
                relative_translation=[0.0, (0.2), -(block_height / 2 - VP_offset)],
                collision_group=6,
            ),
        ],
        translation=[X, Y, Z],
        orientation=[0, 0, 0, 1],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=mass_half_block,
        inertia=inertia_half_block,
    )


#
# Create the projectile
# create a rock
vertices_proj = [
    [-L_b / 2, -L_b / 4, -L_b / 4],
    [-L_b / 2, -L_b / 4, L_b / 4],
    [-L_b / 2, L_b / 4, -L_b / 4],
    [-L_b / 2, L_b / 4, L_b / 4],
    [L_b / 2, -L_b / 4, -L_b / 4],
    [L_b / 2, -L_b / 4, L_b / 4],
    [L_b / 2, L_b / 4, -L_b / 4],
    [L_b / 2, L_b / 4, L_b / 4],
    [-L_b / 4, L_b / 2, -L_b / 4],
    [-L_b / 4, L_b / 2, L_b / 4],
    [L_b / 4, L_b / 2, -L_b / 4],
    [L_b / 4, L_b / 2, L_b / 4],
    [-L_b / 4, -L_b / 2, -L_b / 4],
    [-L_b / 4, -L_b / 2, L_b / 4],
    [L_b / 4, -L_b / 2, -L_b / 4],
    [L_b / 4, -L_b / 2, L_b / 4],
    [-L_b / 4, -L_b / 4, L_b / 2],
    [-L_b / 4, L_b / 4, L_b / 2],
    [L_b / 4, -L_b / 4, L_b / 2],
    [L_b / 4, L_b / 4, L_b / 2],
    [-L_b / 4, -L_b / 4, -L_b / 2],
    [-L_b / 4, L_b / 4, -L_b / 2],
    [L_b / 4, -L_b / 4, -L_b / 2],
    [L_b / 4, L_b / 4, -L_b / 2],
]
#
# computation of the centroid
ch = ConvexHull(vertices_proj)
cm = ch.centroid()
#
# move the vertices to center the center of mass at 0.0
vertices_proj = numpy.array(vertices_proj)[:] - cm[:]
ch = ConvexHull(vertices_proj)
cm = ch.centroid()
#
# computation of inertia and volume
inertia_projectile_1, volume_projectile = ch.inertia(cm)
vx_block = numpy.sqrt(impact_energy * 2.0 / m_projectile)
density_projectile = m_projectile / volume_projectile
inertia_projectile = inertia_projectile_1 * density_projectile
#
#
# Let's start the simulation
t = time.time()
#
# This is where magic happens :)
with MechanicsHdf5Runner() as io:
    # Let's define all the model ingredients #####
    # The concrete block constituents
    io.add_primitive_shape(
        "Box_B",
        "Box",
        [block_depth, display_f * (block_width - block_depth), block_height],
    )
    io.add_primitive_shape(
        "Block_curvature", "Cylinder", (display_f * 0.5 * block_depth, block_height)
    )
    io.add_primitive_shape(
        "Box_half",
        "Box",
        [block_depth, display_f * 0.5 * (block_width - block_depth), block_height],
    )

    # The connector-sling connection replication hollow cylinders
    io.add_primitive_shape("Cyl1", "Cylinder", (cylinder_radius, bar_length / 4.0))
    io.add_primitive_shape("Cyl2", "Cylinder", (cylinder_radius, bar_length / 8.0))
    io.add_primitive_shape(
        "Box_Stop",
        "Box",
        [1.5 * block_hole_diameter, 1.5 * block_hole_diameter, VP_box_depth],
    )

    # The horizontal play thingi
    io.add_mesh_from_file(
        "disk", "data/circle_blocarme.stl", scale=1, insideMargin=0.0, outsideMargin=0.0
    )

    # The vertical play thingi
    io.add_primitive_shape(
        "Box_VP",
        "Box",
        [2.5 * block_hole_diameter, 2.5 * block_hole_diameter, VP_box_depth],
    )
    io.add_primitive_shape("Sp_VP", "Sphere", [factor_bar_diam_ext * bar_diam_ext])

    # The projectile vertices as a convex hull
    io.add_convex_shape("ConvexHull", vertices_proj)

    # The ground... Be humble... :)
    io.add_primitive_shape("Ground", "Box", (10 * block_depth, 25 * block_width, 0.1))

    # Sphere to mark boundary
    io.add_primitive_shape("Sphere", "Sphere", [0.05])

    # Let's create all the model assembly ##### ingrediants to cuisine :)
    # Creating block intertia matrix using the actual shape mesh convex hull
    my_mesh = mesh.Mesh.from_file("data/bloc_arme.stl")
    points = numpy.around(
        numpy.unique(
            my_mesh.vectors.reshape([int(my_mesh.vectors.size / 3), 3]), axis=0
        ),
        2,
    )
    vertices = points.tolist()
    ch = ConvexHull(vertices)
    cm = ch.centroid()
    vertices = numpy.array(vertices)[:] - cm[:]
    ch = ConvexHull(vertices)
    cm = ch.centroid()
    inertia_block1, volume_block = ch.inertia(cm)
    inertia_block = mass_block / volume_block * inertia_block1

    # Creating the half block intertia matrix using the actual half shape mesh convex hull
    # - Left side
    my_mesh_half_left = mesh.Mesh.from_file("data/bloc_arme_half_left.stl")
    points_half_left = numpy.around(
        numpy.unique(
            my_mesh_half_left.vectors.reshape(
                [int(my_mesh_half_left.vectors.size / 3), 3]
            ),
            axis=0,
        ),
        2,
    )
    vertices_half_left = points_half_left.tolist()
    ch_half_left = ConvexHull(vertices_half_left)
    cm_half_left = ch_half_left.centroid()
    vertices_half_left = numpy.array(vertices_half_left)[:] - cm_half_left[:]
    ch_half_left = ConvexHull(vertices_half_left)
    cm_half_left = ch_half_left.centroid()
    inertia_half_block_left_1, volume_block_half_left = ch_half_left.inertia(cm)
    inertia_half_block_left = (
        mass_half_block / volume_block_half_left * inertia_half_block_left_1
    )

    # Creating the half block intertia matrix using the actual half shape mesh convex hull
    #  - Right side
    my_mesh_half_right = mesh.Mesh.from_file("data/bloc_arme_half_right.stl")
    points_half_right = numpy.around(
        numpy.unique(
            my_mesh_half_right.vectors.reshape(
                [int(my_mesh_half_right.vectors.size / 3), 3]
            ),
            axis=0,
        ),
        2,
    )
    vertices_half_right = points_half_right.tolist()
    ch_half_right = ConvexHull(vertices_half_right)
    cm_half_right = ch_half_right.centroid()
    vertices_half_right = numpy.array(vertices_half_right)[:] - cm_half_right[:]
    ch_half_right = ConvexHull(vertices_half_right)
    cm_half_right = ch_half_right.centroid()
    inertia_half_block_right_1, volume_block_half_right = ch_half_right.inertia(cm)
    inertia_half_block_right = (
        mass_half_block / volume_block_half_right * inertia_half_block_right_1
    )

    # Creating blocks pattern as on site ####
    k = 1
    num_pattern = 0
    num_block = 1
    X_pattern = 0.5 * block_depth
    Y_pattern = 0.38
    while num_pattern < Nb_pattern:
        orientation = Orientation_list[num_pattern]
        [num_block, Bottom_block] = make_pattern(
            io,
            X_pattern,
            Y_pattern,
            orientation,
            Nb_block_vert,
            mass_block,
            Bottom_block,
            num_block,
            inertia_block,
        )
        make_bars(io, X_pattern, Y_pattern, num_pattern)
        num_pattern = num_pattern + 1
        X_pattern = X_pattern + (block_width - block_depth) * math.sin(
            orientation
        )  # I know its counterintuitive... :) #Howeever fair
        Y_pattern = Y_pattern + (block_width - block_depth) * math.cos(orientation)
    # The last bar
    make_bars(io, X_pattern, Y_pattern, num_pattern)

    # The four half blocks in bottom and third layers at the left and right extreme edges
    create_half_block_left_bottom(
        io,
        block_depth / 2,
        (0.2 - 0.02),
        0.4,
        1,
        mass_half_block,
        inertia_half_block_left,
    )
    create_half_block_left(
        io,
        block_depth / 2,
        (0.2 - 0.02),
        2.0,
        2,
        mass_half_block,
        inertia_half_block_left,
    )
    create_half_block_right_bottom(
        io,
        block_depth / 2,
        L_structure - (0.2 - 0.02),
        0.4,
        3,
        mass_half_block,
        inertia_half_block_right,
    )
    create_half_block_right(
        io,
        block_depth / 2,
        L_structure - (0.2 - 0.02),
        2.0,
        4,
        mass_half_block,
        inertia_half_block_right,
    )

    # Creating the projectile
    obj_projectile = io.add_object(
        "convexhull2",
        [Contactor("ConvexHull", collision_group=8)],
        translation=[X_block, Y_block, Z_block],
        velocity=[
            vx_block * numpy.cos(angle_impact),
            0,
            vx_block * numpy.sin(angle_impact),
            0,
            0,
            0,
        ],
        mass=m_projectile,
        inertia=inertia_projectile,
    )

    obj_projectile_id = obj_projectile.attrs["id"]
    # print(obj_projectile.attrs['id'])
    # print(obj_projectile.keys())
    # print(obj_projectile.attrs.items())

    # for k in obj_projectile.keys():
    #     print(obj_projectile[k])

    # def print_attrs(name, obj):
    #     print(name)
    #     for key, val in obj.attrs.items():
    #         print("    %s: %s" % (key, val))

    # obj_projectile.visititems(print_attrs)

    # for k, data in obj_projectile.items():
    #     print(k, obj_projectile[k], data, data.attrs)

    # input()
    # Creating ground
    io.add_object(
        "ground",
        [Contactor("Ground", collision_group=9)],
        translation=[0, L_structure / 2, -0.05],
    )

    # Creating spheres to mark the boundaries
    io.add_object(
        "sphere1",
        [Contactor("Sphere", collision_group=10)],
        translation=[0, 0, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )
    io.add_object(
        "sphere2",
        [Contactor("Sphere", collision_group=10)],
        translation=[0, L_structure, 0],
        velocity=[0, 0, 0, 0, 0, 0],
        mass=1,
    )

    # The collision groups are assigned as follows:
    # 0 = Block_cube_portion
    # 1 = Block_curvature_portion
    # 2 = Cylinderical_bars_odd
    # 3 = Cylinderical_bars_even
    # 4 = Bar_bar_penetration_stopper
    # 5 = Disks_of_horizontal_play
    # 6 = Box_of_vertical_play
    # 7 = Sphere_of_vertical_play
    # 8 = Projectile_ETAG
    # 9 = Ground
    # 10 = Spheres_marking_strucuture_edges

    # Let's the component socialise and communicate the mechanics between each other and
    # around the structure's domain ######
    # Establishing the foundation - Putting blocks and cylinderical on ground - here,
    # the cube portion of the block is considered invisible to the ground
    io.add_Newton_impact_friction_nsl(
        "contact_block(curvature)_ground",
        mu=mu_t_cs,
        e=e_t,
        collision_group1=1,
        collision_group2=9,
    )
    io.add_Newton_impact_friction_nsl(
        "contact_cylinder_ground",
        mu=mu_steel,
        e=e_steel,
        collision_group1=2,
        collision_group2=9,
    )

    # newly added contact between projectile and ground
    io.add_Newton_impact_friction_nsl(
        "contact_projectile_ground",
        mu=mu_t_cs,
        e=e_t,
        collision_group1=8,
        collision_group2=9,
    )

    # Establishing connection between blocks - here, the cube portion of the block is
    # considered invisible to each other
    io.add_Newton_impact_friction_nsl(
        "contact_block_block", mu=mu_t, e=e_t, collision_group1=1, collision_group2=1
    )

    # Establishing the contact between cylinderical bars and blocks - through disks
    # incorporating the horizontal play
    io.add_Newton_impact_friction_nsl(
        "contact_bar(odd)_disks",
        mu=mu_steel,
        e=e_steel,
        collision_group1=2,
        collision_group2=5,
    )
    io.add_Newton_impact_friction_nsl(
        "contact_bar(even)_disks",
        mu=mu_steel,
        e=e_steel,
        collision_group1=3,
        collision_group2=5,
    )

    # Making sure that cylinderical bars stay intact in their position without penetrating
    # onto another - thanks to the Bar_bar_penetration_stopper
    io.add_Newton_impact_friction_nsl(
        "contact_bar_stopper(box)",
        mu=0.0,
        e=e_steel,
        collision_group1=7,
        collision_group2=4,
    )

    # Establishing the vertical play by introducing box VP to spherical VP
    io.add_Newton_impact_friction_nsl(
        "contact_VP(sphere)_VP(box)",
        mu=0.0,
        e=e_steel,
        collision_group1=7,
        collision_group2=6,
    )

    # Establishing the contact between the projectile and blocks
    io.add_Newton_impact_friction_nsl(
        "contact_block(cube)_proj",
        mu=mu_t,
        e=e_t,
        collision_group1=0,
        collision_group2=8,
    )
    io.add_Newton_impact_friction_nsl(
        "contact_block(curvature)_proj",
        mu=mu_t,
        e=e_t,
        collision_group1=1,
        collision_group2=8,
    )

# Let's run the model ##########
# Here are the model run options
bullet_options = SiconosBulletOptions()
options = Numerics.SolverOptions(Numerics.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[Numerics.params.SICONOS_IPARAM_MAX_ITER] = 800
options.dparam[Numerics.params.SICONOS_DPARAM_TOL] = 1e-4
options.iparam[Numerics.params.SICONOS_NSGS_FREEZING_CONTACT] = 20

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 1.0
run_options["h"] = 2.5e-4
# run_options['T']= 0.00025 * 500
run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options
run_options["Newton_options"] = siconos.simulation.LINEAR  # TO BE UPDATED

run_options["Newton_max_iter"] = 10
run_options["Newton_tolerance"] = 1e-10
run_options["Newton_warning_on_nonconvergence"] = True
run_options["Warning_nonsmooth_solver"] = True


# run_options['skip_last_update_output']=True
run_options["skip_reset_lambdas"] = True
run_options["osns_assembly_type"] = (
    siconos.nonsmooth_formulations.REDUCED_DIRECT
)  # TO BE UPDATED
# run_options['Newton_max_iter']=1
run_options["theta"] = 0.50001
run_options["verbose"] = True
run_options["with_timer"] = False
run_options["explode_computeOneStep_in_python"] = False
run_options["explode_computeOneStepNSProblem_in_python"] = False
run_options["output_frequency"] = 10
run_options["time_stepping"] = None
# run_options['output_contact_index_set'] = 0 # This is to detect all contacts in the assembly

run_options["constraint_activation_threshold"] = 1e-5


class kinetic_iteration_hook:
    def __init__(self, obj_projectile_id):
        self._io = None
        self._kinetic_sum = []
        self._kinetic_sum_wall = []
        self._normal_work_sum = []
        self._potential_sum = []
        self._tangent_work_sum = []
        self._friction_work_sum = []
        self._normal_work_negative_sum = []
        self._tangent_work_negative_sum = []

        self._obj_projectile_id = obj_projectile_id
        # self._frequency = 1/self._period
        # print('init hook')
        # input()

    def initialize(self, io):
        self._io = io

    def call(self, step):
        # print(' end hook step', step)
        # nsds = self._io._nsds
        positions = self._io._io.positions(self._io._nsds)
        # velocities = self._io._io.velocities(self._io._nsds)
        nsds = self._io._nsds
        ds_idx = positions[0]
        # print(ds_idx)
        kinetic_sum = 0.0
        kinetic_sum_wall = 0.0
        potential_sum = 0.0
        for i in ds_idx:
            n_ds = int(i)

            neds = nsds.dynamicalSystem(n_ds)
            # neds= siconos.modeling.cast_NewtonEulerDS(ds) # TO BE UPDATED
            kinetic = neds.computeKineticEnergy()  # TO BE UPDATED
            time = self._io._simulation.startingTime()
            # print('kinetic', kinetic, 'time ', self._io._simulation.startingTime())
            if n_ds == self._obj_projectile_id:
                print("kinetic obj ", kinetic)
            else:
                kinetic_sum_wall = kinetic_sum_wall + kinetic

            kinetic_sum = kinetic + kinetic_sum
            m = neds.scalarMass  # TO BE UPDATED

            pos_z = neds.q()[2]  # TO BE UPDATED
            mass = m
            potential = -mass * 9.81 * pos_z
            # print('potential', potential, 'time ', self._io._simulation.startingTime())
            potential_sum = potential + potential_sum

        self._kinetic_sum.append([time, kinetic_sum])
        self._kinetic_sum_wall.append([time, kinetic_sum_wall])
        self._potential_sum.append([time, potential_sum])
        print("kinetic_sum", kinetic_sum)
        print("potential_sum", potential_sum)
        #        input()

        cf_work = self._io._io.contactContactWork(self._io._nsds, 1)
        if cf_work is not None:
            # print('cf_work', cf_work)

            normal_work = cf_work[:, 1]
            normal_work_negative = np.where(normal_work < 0, normal_work, 0)
            # print('normal_work_negative', normal_work_negative)

            tangent_work = cf_work[:, 2]
            tangent_work_negative = np.where(tangent_work < 0, tangent_work, 0)
            # print('tangent_work_negative', tangent_work_negative)

            normal_work_negative_sum = np.sum(normal_work_negative)
            tangent_work_negative_sum = np.sum(tangent_work_negative)
            friction_work_sum = np.sum(cf_work[:, 3])

            normal_work_sum = np.sum(normal_work)
            tangent_work_sum = np.sum(tangent_work)
            friction_work_sum = np.sum(cf_work[:, 3])

            print("normal_work_sum", normal_work_sum)
            print("tangent_work_sum", tangent_work_sum)
            print("friction_work_sum", friction_work_sum)

            self._normal_work_sum.append(normal_work_sum)
            self._tangent_work_sum.append(tangent_work_sum)

            self._normal_work_negative_sum.append(normal_work_negative_sum)
            self._tangent_work_negative_sum.append(tangent_work_negative_sum)

            self._friction_work_sum.append(friction_work_sum)
        else:
            self._normal_work_sum.append(0.0)
            self._tangent_work_sum.append(0.0)

            self._normal_work_negative_sum.append(0.0)
            self._tangent_work_negative_sum.append(0.0)

            self._friction_work_sum.append(0.0)
        # input()


before_next_step_hook = kinetic_iteration_hook(obj_projectile_id)


run_options["before_next_step_iteration_hook"] = before_next_step_hook

# This gets things going
with MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)

# print('kinetic energy', before_next_step_hook._kinetic_sum)
# print('kinetic energy wall ', before_next_step_hook._kinetic_sum_wall)
directory = "nonlinear_10"
directory = "nonlinear_spurious"
directory = "nonlinear_activation_1e-05"

if not os.path.exists(directory):
    os.makedirs(directory)

np.save(
    os.path.join(directory, "kinetic_energy.npy"), before_next_step_hook._kinetic_sum
)
np.save(
    os.path.join(directory, "kinetic_energy_wall.npy"),
    before_next_step_hook._kinetic_sum_wall,
)
np.save(
    os.path.join(directory, "potential_energy.npy"),
    before_next_step_hook._potential_sum,
)


np.save(
    os.path.join(directory, "normal_work.npy"), before_next_step_hook._normal_work_sum
)
np.save(
    os.path.join(directory, "tangent_work.npy"), before_next_step_hook._tangent_work_sum
)
np.save(
    os.path.join(directory, "friction_work.npy"),
    before_next_step_hook._friction_work_sum,
)


np.save(
    os.path.join(directory, "normal_work_negative.npy"),
    before_next_step_hook._normal_work_negative_sum,
)
np.save(
    os.path.join(directory, "tangent_work_negative.npy"),
    before_next_step_hook._tangent_work_negative_sum,
)


elapsed = round((time.time() - t) / 60, 2)
print("The simulation time is", elapsed, " minutes")
