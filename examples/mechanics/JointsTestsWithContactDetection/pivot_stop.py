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
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
import siconos.simulation
import siconos.integrators


import siconos.numerics as sn
import numpy as np

# Configuration of two stops, at 3/4 pi on either side of the initial position.
# The values are [axis, position, direction], where direction must be -1 or 1.
stops = [[0, -np.pi * 3 / 4, 1], [0, np.pi * 3 / 4, -1]]

# Initial rotational velocity of the body attached to the joint.
twist = 10

with MechanicsHdf5Runner() as io:
    # self-collide property: if set to true, collisions between bodies
    # connected by joints are allowed.  In this example, the joint
    # keeps the bodies in an overlapping state, which leads to
    # ill-posed problems for the solver.
    self_collide = False

    # A "bar" connected to a "post".  A "knob" is attached to the post
    # for visual reference -- otherwise no rotation can be seen, since
    # it is cylindrical.
    io.add_primitive_shape("Bar", "Box", (1, 0.1, 0.1))
    io.add_primitive_shape("Post", "Cylinder", (0.05, 1))
    io.add_primitive_shape("Knob", "Box", (0.2, 0.05, 0.05))
    io.add_primitive_shape("Ground", "Box", (4, 4, 0.5))

    # Ground is defined for visual reference
    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, 0])

    # We define a contact law even though this simulation should not
    # feature contact.
    io.add_Newton_impact_friction_nsl("contact", e=0.7, mu=0.02)

    # This law is used to specify "bouncy" stops on joint1.
    io.add_Newton_impact_nsl("stop", e=0.8)

    # Very high friction on the second joint causes "almost-fixed"
    # behaviour, resulting in a bounce when the first joint hits the
    # stop, with energy partially absorbed by a small movement of the
    # second joint.
    io.add_relay_nsl("friction", lb=-3.0, ub=3.0)

    # The objects, with self-collision disabled as noted above.
    io.add_object(
        "bar",
        [Contactor("Bar")],
        translation=[0.45, 0.45, 3],
        mass=10,
        allow_self_collide=self_collide,
        velocity=[0, 0, 0, 0, twist, 0],
    )
    io.add_object(
        "post",
        [Contactor("Post"), Contactor("Knob", relative_translation=[0.1, 0, 0])],
        translation=[0, 0, 3],
        mass=1,
        allow_self_collide=self_collide,
    )

    # Connect the two bodies by a pivot joint
    io.add_joint(
        "joint1",
        "bar",
        "post",
        [[-0.45, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        allow_self_collide=self_collide,
        nslaws="stop",
        stops=stops,
        absolute=False,
    )

    # Joint from "bar" to the world reference frame, to keep things from falling.
    io.add_joint(
        "joint2",
        "post",
        None,
        [[0, 0, 0]],
        [[0, 1, 0]],
        "PivotJointR",
        absolute=False,
        friction="friction",
    )

    # For fully fixed behaviour replace with a FixedJointR.
    # io.add_joint('joint2', 'post', None, None, None, 'FixedJointR')


# We define a "controller" here to show how to measure the angle of
# the joint using computehDoF.
class Ctrl(object):
    def initialize(self, io):
        self.nsds = io._nsds
        self.topo = self.nsds.topology()
        self.joint1_inter = self.topo.getInteraction("joint1")
        self.joint1 = self.joint1_inter.relation()
        self.bar = self.topo.getDynamicalSystem("bar")
        self.post = self.topo.getDynamicalSystem("post")
        self.y = np.zeros(5)
        self.yDoF = np.zeros(1)
        self.jachq = np.zeros((1, 14))

    def step(self):
        self.joint1.computeh(self.bar.q(), self.post.q(), self.y)
        self.joint1.computehDoF(self.bar.q(), self.post.q(), self.yDoF)
        self.joint1.computeJachqDoF(
            self.joint1_inter, self.bar.q(), self.post.q(), self.jachq, 0
        )
        print("joint angle", self.yDoF)


options = sn.SolverOptions(sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10
options.update_internal(1, sn.solver_ids.SICONOS_FRICTION_3D_ONECONTACT_NSN)


options = None


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 5.0
run_options["h"] = 0.001
run_options["theta"] = 0.50001
run_options["Newton_max_iter"] = 1
run_options["solver_options"] = options
run_options["controller"] = Ctrl()
run_options["projection_itermax"] = 3
run_options["projection_tolerance"] = 1e-5
run_options["projection_tolerance_unilateral"] = 1e-5
run_options["time_stepping"] = siconos.simulation.TimeSteppingDirectProjection
run_options["osi"] = siconos.integrators.MoreauJeanDirectProjectionOSI
run_options["numerics_verbose"] = False

# Run the simulation
with siconos.io.mechanics_run.MechanicsHdf5Runner(mode="r+") as io:
    io.run(run_options)
