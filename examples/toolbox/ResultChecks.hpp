/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef SICONOS_EXAMPLES_RESULT_CHECKS_HPP
#define SICONOS_EXAMPLES_RESULT_CHECKS_HPP

/** Checks of the results of an example without a reference trajectory.

    A row-by-row comparison with a reference file (siconos::algebra::io::compareRefFile)
    is the right test when the results are determined by the kernel and numerics only.
    It is not when they also depend on an external library (contact detection by Bullet
    or OpenCascade: margins, contact points, version) or on round-off (a contact at a
    gap of 1e-17): a contact detected one step later shifts the whole trajectory.

    For these examples, check instead:
    - features of the trajectory that are robust to such shifts, with explicit
      tolerances: time of the first impact (+/- a few steps), restitution of the first
      impact, final position at rest;
    - physical properties that need no reference: bounded penetration, energy that does
      not increase, no overlap of bodies, bodies at rest at the end.

    Usage (include it with a path relative to the example, e.g.
    #include "../../toolbox/ResultChecks.hpp"):

    \code
    siconos::examples::Checks chk("BulletBouncingBox");
    chk.near("first impact time", t_impact, std::sqrt(2 * z0 / g), 2 * h);
    chk.atLeast("lowest point of the box", min_z, -0.5 - 0.1);
    chk.atMost("energy increase in one step", max_increase, 1e-6 * E0);
    return chk.report();  // prints the table, returns 1 if a check failed
    \endcode

    The Python twin is toolbox/result_checks.py.
*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace siconos {
namespace examples {

class Checks {
 public:
  explicit Checks(std::string name) : _name(std::move(name)) {}

  /** |value - expected| <= tol */
  void near(const std::string& what, double value, double expected, double tol) {
    add(what, value, std::abs(value - expected) <= tol,
        "expected " + fmt(expected) + " +/- " + fmt(tol));
  }
  /** value >= bound */
  void atLeast(const std::string& what, double value, double bound) {
    add(what, value, value >= bound, ">= " + fmt(bound));
  }
  /** value <= bound */
  void atMost(const std::string& what, double value, double bound) {
    add(what, value, value <= bound, "<= " + fmt(bound));
  }
  /** a condition that has no single value (e.g. an event that must occur) */
  void isTrue(const std::string& what, bool ok, const std::string& detail = "") {
    _rows.push_back({what, detail, "", ok});
  }

  bool ok() const {
    return std::all_of(_rows.begin(), _rows.end(), [](const Row& r) { return r.ok; });
  }

  /** Print the table of the checks; return 0 if all passed, 1 otherwise. */
  int report() const {
    std::printf("\n====> Checks of the results of %s\n", _name.c_str());
    for (const auto& r : _rows)
      std::printf("  [%s] %-45s %-14s %s\n", r.ok ? " ok " : "FAIL", r.what.c_str(),
                  r.value.c_str(), r.expected.c_str());
    std::printf("====> %s\n", ok() ? "all checks passed" : "SOME CHECKS FAILED");
    return ok() ? 0 : 1;
  }

 private:
  struct Row {
    std::string what, value, expected;
    bool ok;
  };
  static std::string fmt(double x) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", x);
    return buf;
  }
  void add(const std::string& what, double value, bool ok, const std::string& expected) {
    _rows.push_back({what, fmt(value), expected, ok});
  }
  std::string _name;
  std::vector<Row> _rows;
};

/** Running maximum of the increase of a quantity between two successive calls
    (e.g. an energy, which must not increase). */
class MaxIncrease {
 public:
  void push(double x) {
    if (_started) _max = std::max(_max, x - _last);
    _last = x;
    _started = true;
  }
  double value() const { return _max; }
  /** next reference value, without counting the change (e.g. a body was added) */
  void restart(double x) {
    _last = x;
    _started = true;
  }

 private:
  double _last{0.}, _max{0.};
  bool _started{false};
};

/** Rigid body with the Siconos (NewtonEuler) conventions:
    q = (x, y, z, q0, q1, q2, q3) (position, unit quaternion), twist = (v, omega) with the
    velocity v in the inertial frame and the angular velocity omega in the body frame. */
namespace rigid_body {

/** z of the lowest corner of a box of half sizes (a, b, c), centred at q */
template <class Q>
double boxLowestZ(const Q& q, double a, double b, double c) {
  const double w = q(3), x = q(4), y = q(5), z = q(6);
  // third row of the rotation matrix
  const double r20 = 2 * (x * z - w * y), r21 = 2 * (y * z + w * x),
               r22 = 1 - 2 * (x * x + y * y);
  return q(2) - (std::abs(r20) * a + std::abs(r21) * b + std::abs(r22) * c);
}

/** max_j |R_{2j}| (R: rotation of the body): 1 when one of the body axes is vertical,
    e.g. a box lying on a face */
template <class Q>
double verticalAxisCosine(const Q& q) {
  const double w = q(3), x = q(4), y = q(5), z = q(6);
  return std::max({std::abs(2 * (x * z - w * y)), std::abs(2 * (y * z + w * x)),
                   std::abs(1 - 2 * (x * x + y * y))});
}

/** kinetic + gravity potential energy: 1/2 twist^T M twist + m g z, with M the 6x6 mass
    matrix of the body (NewtonEulerDS::totalInertiaMatrix()). Use the matrix of the body,
    not the inertia given to its constructor: e.g. the Bullet collision manager replaces
    the inertia of a RigidBodyDS by the one of its contactors (useContactorInertia(),
    true by default; margins included). */
template <class Q, class V, class M>
double energy(const Q& q, const V& twist, const M& massMatrix, double g) {
  return 0.5 * twist.dot(massMatrix * twist) + massMatrix(0, 0) * g * q(2);
}

/** distance between the centres of two bodies */
template <class Q1, class Q2>
double centreDistance(const Q1& q1, const Q2& q2) {
  return std::hypot(q1(0) - q2(0), q1(1) - q2(1), q1(2) - q2(2));
}

}  // namespace rigid_body
}  // namespace examples
}  // namespace siconos

#endif
