/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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
#ifndef ROCKING_R_HPP
#define ROCKING_R_HPP

#include <BlockVector.hpp>
#include <LagrangianScleronomousR.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>

namespace user_defined {
class RockingBlockR : public siconos::modeling::LagrangianScleronomousR {
 public:
  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  RockingBlockR() : LagrangianScleronomousR{} { hasHdot_ = true; };

  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) {
    double q1 = q(1);
    double q2 = q(2);
    y(0) = q1 - 0.5 * LengthBlock * sin(q2) - 0.5 * HeightBlock * cos(q2);
  }

  void computeH(const siconos::algebra::BlockVector& q) {
    double q2 = q(2);
    H_view_->setValue(0, 0, 0.0);
    H_view_->setValue(0, 1, 1.0);
    H_view_->setValue(
        0, 2, -0.5 * LengthBlock * cos(q2) + 0.5 * HeightBlock * sin(q2));
  }

  void computejacobianhOver_q_dot(const siconos::algebra::BlockVector& q,
                                  const siconos::algebra::BlockVector& qdot) {
    double q2 = q(2);
    double qdot2 = qdot(2);
    Hdot_->setValue(0, 0, 0.0);
    Hdot_->setValue(0, 1, 0.0);
    Hdot_->setValue(
        0, 2, (0.5 * LengthBlock * sin(q2) + 0.5 * HeightBlock * cos(q2)) * qdot2);
  }
};

class RockingBlockR2 : public siconos::modeling::LagrangianScleronomousR {
 public:
  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  RockingBlockR2() : LagrangianScleronomousR{} { hasHdot_ = true; };

  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) {
    double q1 = q(1);
    double q2 = q(2);
    y(0) = q1 + 0.5 * LengthBlock * sin(q2) - 0.5 * HeightBlock * cos(q2);
  }

  void computeH(const siconos::algebra::BlockVector& q) {
    double q2 = q(2);
    H_view_->setValue(0, 0, 0.0);
    H_view_->setValue(0, 1, 1.0);
    H_view_->setValue(0, 2,
                                    0.5 * LengthBlock * cos(q2) + 0.5 * HeightBlock * sin(q2));
  }

  void computeDotJachqcomputejacobianhOver_q_dot(const siconos::algebra::BlockVector& q,
                                                 const siconos::algebra::BlockVector& qdot) {
    double q2 = q(2);
    double qdot2 = qdot(2);
    Hdot_->setValue(0, 0, 0.0);
    Hdot_->setValue(0, 1, 0.0);
    Hdot_->setValue(
        0, 2, (-0.5 * LengthBlock * sin(q2) + 0.5 * HeightBlock * cos(q2)) * qdot2);
  }
};
}  // namespace user_defined
#endif
