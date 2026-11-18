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

/*!\file BouncingBallNETS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, O. Bonnefon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>

#ifdef WITH_FC3D
#define R_CLASS NewtonEulerVelocityR
#else
#define R_CLASS NewtonEulerNormalVelocityR
#endif

#ifndef USER_REL
#define USER_REL

namespace user_defined {
class my_NewtonEulerR : public siconos::modeling::R_CLASS {
  double _sBallRadius;

 public:
  my_NewtonEulerR(double radius) : R_CLASS(), _sBallRadius(radius){};

  virtual void computeOutput(double time, siconos::modeling::Interaction& inter,
                             siconos::algebra::blocks::size_type derivativeNumber) override {
    const auto& ds_vars = inter.read_dynamical_systems_variables();
    if (derivativeNumber == 0) {
      auto q0 = ds_vars[siconos::tools::enum_to_index(ds_var::q0)]->vector(0);
      auto q1 = ds_vars[siconos::tools::enum_to_index(ds_var::q0)]->vector(1);
      if (q1)
        computeh(*q0, *q1, *inter.y(0));
      else
        computeh(*q0, std::nullopt, *inter.y(0));
    } else {
      R_CLASS::computeOutput(time, inter, derivativeNumber);
    }
  }

  void computeh(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
                const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override {
    assert(q2);
    double height = q1(0) - _sBallRadius - (*q2)(0);
    y(0) = height;
    nc_ << 1., 0., 0.;
    contactPoint1_(0) = q1(0) - _sBallRadius;
    contactPoint1_(1) = q1(1);
    contactPoint1_(2) = q1(2);

    contactPoint2_ = q2->head(3);

    // printf("my_NewtonEulerR N, Pc\n");
    // siconos::algebra::print(nc_);
    // siconos::algebra::print(contactPoint1_);
    // siconos::algebra::print(contactPoint2_);
    // std::cout << "my_NewtonEulerR:: computeh ends" << std::endl;
  }
};
}  // namespace user_defined
#endif