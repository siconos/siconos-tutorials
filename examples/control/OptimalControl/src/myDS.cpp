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

#include "myDS.h"

// #define SICONOS_DEBUG

user_defined::MyDS::MyDS(Eigen::Ref<siconos::algebra::SiconosVector> x0)
    : FirstOrderNonLinearDS(x0, siconos::algebra::copy_t) {
  setComputefVectorFunction(
      [this](const Eigen::Ref<const siconos::algebra::SiconosVector>& state, double time,
             Eigen::Ref<siconos::algebra::MapVectorType> result) {
        siconos::algebra::SiconosVector X{2};
        X(0) = (state(0) - 2.0);
        X(1) = (state(1) + 1.0);
        auto QX = *Q * X;

        siconos::algebra::SiconosVector P{2};
        P(0) = state(2);
        P(1) = state(3);
        auto K1P = *K1 * P;
        auto alphatmp = alpha(time, state);

        result(0) = alphatmp(0);
        result(1) = alphatmp(1);
        result(2) = -QX(0) + K1P(0);
        result(3) = -QX(1) + K1P(1);
      });

  setComputeJacobianfOver_xFunction(
      [this](const Eigen::Ref<const siconos::algebra::SiconosVector>& state, double time,
             Eigen::Ref<siconos::algebra::MapType> result) {
        auto jacXalpha = JacobianXalpha(time, state);
        result.setZero();
        result(0, 0) = jacXalpha(0, 0);
        result(0, 1) = jacXalpha(0, 1);
        result(1, 0) = jacXalpha(1, 0);
        result(1, 1) = jacXalpha(1, 1);
        result(2, 0) = -(*Q)(0, 0);
        result(2, 1) = -(*Q)(0, 1);
        result(2, 2) = (*K1)(0, 0);
        result(2, 3) = (*K1)(0, 1);
        result(3, 0) = -(*Q)(1, 0);
        result(3, 1) = -(*Q)(1, 1);
        result(3, 2) = (*K1)(1, 0);
        result(3, 3) = (*K1)(1, 1);
      });

  Q = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  Q->setIdentity();
  K1 = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, 1.0 / 2.0);
  K1->setValue(1, 0, -1.0 / 2.0);
  K1->setValue(1, 1, +1.0);
}

siconos::algebra::SiconosVector user_defined::alpha(
    double t, const Eigen::Ref<const siconos::algebra::SiconosVector>& state) {
  siconos::algebra::SiconosVector res{2};
  res(0) = 1.0 / 2.0 * state(1) + 1.0 / 2.0;
  res(1) = -1.0 / 2.0 * state(0) - state(1);
  return res;  // RVO
}

siconos::algebra::SiconosMatrix user_defined::JacobianXalpha(
    double t, const Eigen::Ref<const siconos::algebra::SiconosVector>& state) {
  siconos::algebra::SiconosMatrix res{2, 2};
  res(0, 0) = 0.0;
  res(0, 1) = 1.0 / 2.0;
  res(1, 0) = -1.0 / 2.0;
  res(1, 1) = -1.0;

#ifdef SICONOS_DEBUG
  std::cout << "JacXalpha\n" << std::endl;
  ;
  siconos::algebra::print(res);
#endif
  return res;
}
