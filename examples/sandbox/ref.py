#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2024 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
# make
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
import numpy as np
import siconos.pb11_template as tp

ndof = 3
initial_position = np.array([1, 0, 0], dtype=np.float64)

ball = tp.ClassA(initial_position)


fext = np.zeros_like(initial_position)
fext[:] = 122
ball.setConstantVector2(fext)

print(ball.vector2())

fext[1] = 12
print(ball.vector2())

## ball.vector2()[2] = 1.2  Read-only !

# print(ball.vector2())
# print(fext)

fext = 12
print(ball.vector2())
print(fext)


# ndof =

ball2 = tp.ClassA(initial_position)


def external_forces(time, fext):
    for i in range(ndof):
        fext[i] = time + i


ball2.setComputeVector2Function(external_forces)

ball2.computeVector2(1.0)


print(ball2.vector2())

ball2.computeVector2(14.0)


print(ball2.vector2())


# ndof = 100000000

q0 = np.zeros(ndof, dtype=np.float64)

ball3 = tp.ClassA(q0)


fext2 = np.zeros_like(q0)
fext2[:] = 122
ball3.setConstantVector2(fext2)
ball3.setComputeVector2Function(external_forces)


ball3.computeVector2(1.0)

pos = np.zeros(ndof)
ball3.computeMatrix1(pos, 1)  # Work but does nothing

mass = np.zeros((ndof, ndof), dtype=np.float64, order="F")
ball3.setConstantMatrix1(mass)

print(ball3.matrix1)

#  ball3.matrix1[1,1] = 12 # Read-only Must fail !!

print(ball3.matrix1)


ball3.computeMatrix1(pos, 1)
print(ball3.matrix1)


def mass(time, pos, mass_storage):
    mass_storage[...] = 4


ball3.setComputeMatrix1Function(mass)

ball3.computeMatrix1(pos, 1.0)

print(ball3.matrix1)
