"""Checks of the results of an example without a reference trajectory.

Python twin of toolbox/ResultChecks.hpp (see its documentation): when the results depend
on an external library (contact detection by Bullet, ...) or on round-off, compare
features of the trajectory with explicit tolerances and check physical properties,
instead of comparing the whole trajectory with a reference file.

Usage, from an example directory two levels below examples/::

    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "toolbox"))
    from result_checks import Checks

    chk = Checks("BulletBouncingBox")
    chk.near("first impact time", t_impact, (2 * z0 / g) ** 0.5, 2 * h)
    sys.exit(chk.report())  # prints the table, 1 if a check failed
"""

import numpy as np


class Checks:
    def __init__(self, name):
        self.name = name
        self.rows = []

    def near(self, what, value, expected, tol):
        """|value - expected| <= tol"""
        self.rows.append(
            (
                what,
                value,
                abs(value - expected) <= tol,
                f"expected {expected:.6g} +/- {tol:.6g}",
            )
        )

    def at_least(self, what, value, bound):
        self.rows.append((what, value, value >= bound, f">= {bound:.6g}"))

    def at_most(self, what, value, bound):
        self.rows.append((what, value, value <= bound, f"<= {bound:.6g}"))

    def ok(self):
        return all(r[2] for r in self.rows)

    def report(self):
        """Print the table of the checks; return 0 if all passed, 1 otherwise."""
        print(f"\n====> Checks of the results of {self.name}")
        for what, value, ok, expected in self.rows:
            print(
                f"  [{' ok ' if ok else 'FAIL'}] {what:45s} {value:<14.6g} {expected}"
            )
        print("====> " + ("all checks passed" if self.ok() else "SOME CHECKS FAILED"))
        return 0 if self.ok() else 1


class MaxIncrease:
    """Running maximum of the increase of a quantity between two successive pushes."""

    def __init__(self):
        self.last = None
        self.value = 0.0

    def push(self, x):
        if self.last is not None:
            self.value = max(self.value, x - self.last)
        self.last = x

    def restart(self, x):
        """next reference value, without counting the change (e.g. a body was added)"""
        self.last = x


# Rigid bodies, Siconos (NewtonEuler) conventions: q = (x, y, z, q0, q1, q2, q3),
# twist = (v, omega), v in the inertial frame, omega in the body frame.


def _third_row(q):
    w, x, y, z = q[3:7]
    return np.array([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)])


def box_lowest_z(q, a, b, c):
    """z of the lowest corner of a box of half sizes (a, b, c), centred at q"""
    return q[2] - np.abs(_third_row(q)) @ np.array([a, b, c])


def vertical_axis_cosine(q):
    """max_j |R_2j|: 1 when one of the body axes is vertical (box lying on a face)"""
    return np.abs(_third_row(q)).max()


def energy(q, twist, mass_matrix, g):
    """1/2 twist^T M twist + m g z, M = body.totalInertiaMatrix (6x6). Use the matrix of
    the body: the Bullet collision manager replaces the inertia given to a RigidBodyDS by
    the one of its contactors (useContactorInertia, true by default)."""
    twist = np.asarray(twist)
    return 0.5 * twist @ (mass_matrix @ twist) + mass_matrix[0, 0] * g * q[2]


def centre_distance(q1, q2):
    return float(np.linalg.norm(np.asarray(q1[:3]) - np.asarray(q2[:3])))
