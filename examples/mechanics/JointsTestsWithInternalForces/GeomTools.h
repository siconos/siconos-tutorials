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

#ifndef GEOMTOOLS_H
#define GEOMTOOLS_H
#include <boost/math/quaternion.hpp>

namespace geomtools {
/* Given a position of a point in the Inertial Frame and the configuration vector q of a solid
 * returns a position in the spatial frame.
 */
void fromInertialToSpatialFrame(double *positionInInertialFrame,
                                double *positionInSpatialFrame, auto q) {
  double q0 = q->getValue(3);
  double q1 = q->getValue(4);
  double q2 = q->getValue(5);
  double q3 = q->getValue(6);

  boost::math::quaternion<double> quatQ(q0, q1, q2, q3);
  boost::math::quaternion<double> quatcQ(q0, -q1, -q2, -q3);
  boost::math::quaternion<double> quatpos(
      0, positionInInertialFrame[0], positionInInertialFrame[1], positionInInertialFrame[2]);
  boost::math::quaternion<double> quatBuff;

  // perform the rotation
  quatBuff = quatQ * quatpos * quatcQ;

  positionInSpatialFrame[0] = quatBuff.R_component_2() + q->getValue(0);
  positionInSpatialFrame[1] = quatBuff.R_component_3() + q->getValue(1);
  positionInSpatialFrame[2] = quatBuff.R_component_4() + q->getValue(2);
}
void tipTrajectories(auto q, double *traj, double length) {
  double positionInInertialFrame[3];
  double positionInSpatialFrame[3];
  // Output the position of the tip of beam1
  positionInInertialFrame[0] = length / 2;
  positionInInertialFrame[1] = 0.0;
  positionInInertialFrame[2] = 0.0;

  fromInertialToSpatialFrame(positionInInertialFrame, positionInSpatialFrame, q);
  traj[0] = positionInSpatialFrame[0];
  traj[1] = positionInSpatialFrame[1];
  traj[2] = positionInSpatialFrame[2];

  // std::cout <<  "positionInSpatialFrame[0]" <<  positionInSpatialFrame[0]<<std::endl;
  // std::cout <<  "positionInSpatialFrame[1]" <<  positionInSpatialFrame[1]<<std::endl;
  // std::cout <<  "positionInSpatialFrame[2]" <<  positionInSpatialFrame[2]<<std::endl;

  positionInInertialFrame[0] = -length / 2;
  fromInertialToSpatialFrame(positionInInertialFrame, positionInSpatialFrame, q);
  traj[3] = positionInSpatialFrame[0];
  traj[4] = positionInSpatialFrame[1];
  traj[5] = positionInSpatialFrame[2];
}

}  // namespace geomtools

#endif
