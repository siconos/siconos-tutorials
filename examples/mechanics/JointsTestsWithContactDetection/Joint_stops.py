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


from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn

with MechanicsHdf5Runner() as io:
    io.add_primitive_shape("Box", "Box", (1, 1, 1))
    io.add_primitive_shape("Ground", "Box", (10, 30, 0.5))
    io.add_Newton_impact_friction_nsl("contact", e=0.5, mu=0.5)
    io.add_Newton_impact_nsl("stop", e=0.8)
    # not a contact NSL by default (group1,2==-1)
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, -0.25])

    # We enable self-collide on the objects and their common joint to
    # show how contact behaviour and stops can interact.

    io.add_object(
        "cube1",
        [Contactor("Box")],
        translation=[-1, 5, 2.5],
        mass=0.1,
        velocity=[0, 0, 20, 0, 0, 0],
        allow_self_collide=True,
    )

    io.add_object(
        "cube2",
        [Contactor("Box")],
        translation=[-1, 5, 4.0],
        mass=0.1,
        velocity=[0, 0, 0, 0, 0, 0],
        allow_self_collide=True,
    )

    # We add two stops to each joint, one positive and the other
    # negative The initial velocity causes the cube to fly upwards,
    # hit the stop, and then bounce back downwards and hit the other
    # stop.  The stops are positioned at 2.0 and -2.0 and are on axis
    # 0 (the only DoF for prismatic joint.)  The prismatic vectors are
    # diagonal for joint1, and vertical for joint2, to demonstrate the
    # joint dynamics.

    io.add_joint(
        "joint1",
        "cube1",
        None,
        None,
        [[0, 0.707, 0.707]],
        "PrismaticJointR",
        nslaws="stop",
        stops=[[0, 2.0, -1], [0, -2.0, 1]],
    )

    io.add_joint(
        "joint2",
        "cube1",
        "cube2",
        None,
        [[0, 0, 1]],
        "PrismaticJointR",
        nslaws="stop",
        stops=[[0, 2.0, -1], [0, -2.0, 1]],
        allow_self_collide=True,
    )

options = sn.SolverOptions(sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-4
# projection_itermax=3,
# projection_tolerance=1e-5,
# projection_tolerance_unilateral=1e-5,

with MechanicsHdf5Runner(mode="r+") as io:
    io.run(
        t0=0,
        T=5,
        h=0.001,
        theta=0.50001,
        Newton_max_iter=2,
        solver_options=options,
        # time_stepping=siconos.simulation.TimeSteppingDirectProjection,
        # osi=siconos.integrators.MoreauJeanDirectProjectionOSI,
    )
