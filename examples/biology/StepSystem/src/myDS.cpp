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
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <siconos_debug.h>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

user_defined::MyDS::MyDS(Eigen::Ref<siconos::algebra::SiconosVector> x0)
    : FirstOrderNonLinearDS(x0, siconos::algebra::alias_t) {
  setComputefVectorFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& x,
                               double time,
                               Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result(0) = -4.5 * x(0);
    result(1) = -1.5 * x(1);
  });

  setComputeJacobianfOver_xFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& x, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
        result(0, 0) = -4.5;
        result(1, 1) = -1.5;
      });
}
