/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <iostream>
#include <math.h>

#include "SCConst.h" // Simulation parameters
#include "SiconosException.hpp"

#undef restrict
#define restrict __restrict

// plugins for smooth equations of motion
SICONOS_EXPORT void mass(unsigned int sizeOfq, double *restrict q,
                         double *restrict mass, unsigned int sizeZ,
                         double *restrict z) {
  // columnwise definition of mass matrix
  mass[0] = parameters::J1 +
            (0.25 * parameters::m1 + parameters::m2 + parameters::m3) *
                parameters::l1 * parameters::l1;
  mass[1] = (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
            parameters::l2 * cos(q[1] - q[0]);
  mass[2] = 0.;

  mass[3] = (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
            parameters::l2 * cos(q[1] - q[0]);
  mass[4] = parameters::J2 + (0.25 * parameters::m2 + parameters::m3) *
                                 parameters::l2 * parameters::l2;
  mass[5] = 0.;

  mass[6] = 0.;
  mass[7] = 0.;
  mass[8] = parameters::J3;
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, double *restrict q,
                         double *restrict velocity, double *restrict FGyr,
                         unsigned int sizeZ, double *restrict z) {
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  FGyr[0] = (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
            parameters::l2 * sin(q[0] - q[1]) * velocity[1] * velocity[1];
  FGyr[1] = -(0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
            parameters::l2 * sin(q[0] - q[1]) * velocity[0] * velocity[0];
  FGyr[2] = 0.;
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, double *restrict q,
                                  double *restrict velocity,
                                  double *restrict jacob, unsigned int sizeZ,
                                  double *restrict z) {
  // Jacobian of nonlinear inertia terms with respect to q (columnwise)
  jacob[0] = (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * cos(q[0] - q[1]) * velocity[1] * velocity[1];
  jacob[1] = -(0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * cos(q[0] - q[1]) * velocity[0] * velocity[0];
  jacob[2] = 0.;

  jacob[3] = -(0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * cos(q[0] - q[1]) * velocity[1] * velocity[1];
  jacob[4] = (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * cos(q[0] - q[1]) * velocity[0] * velocity[0];
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void jacobianFGyrqDot(unsigned int sizeOfq, double *restrict q,
                                     double *restrict velocity,
                                     double *restrict jacob, unsigned int sizeZ,
                                     double *restrict z) {
  // Jacobian of nonlinear inertia terms with respect to velocity (columnwise)
  jacob[0] = 0.;
  jacob[1] = -2. * (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * sin(q[0] - q[1]) * velocity[0];
  jacob[2] = 0.;

  jacob[3] = 2. * (0.5 * parameters::m2 + parameters::m3) * parameters::l1 *
             parameters::l2 * sin(q[0] - q[1]) * velocity[1];
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, double *restrict q,
                         double *restrict velocity, double *restrict fInt,
                         unsigned int sizeZ, double *restrict z) {
  // internal forces (negative in h according to LagrangianDS)
  fInt[0] = (0.5 * parameters::m1 + parameters::m2 + parameters::m3) *
            parameters::gravity * parameters::l1 * cos(q[0]);
  fInt[1] = (0.5 * parameters::m2 + parameters::m3) * parameters::gravity *
            parameters::l2 * cos(q[1]);
  fInt[2] = 0.;
}

SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq,
                                  double *restrict q, double *restrict velocity,
                                  double *restrict jacob, unsigned int sizeZ,
                                  double *restrict z) {
  // Jacobian of internal forces with respect to q (columnwise)
  jacob[0] = -(0.5 * parameters::m1 + parameters::m2 + parameters::m3) *
             parameters::gravity * parameters::l1 * sin(q[0]);
  jacob[1] = 0.;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = -(0.5 * parameters::m2 + parameters::m3) * parameters::gravity *
             parameters::l2 * sin(q[1]);
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void jacobianFIntqDot(double time, unsigned int sizeOfq,
                                     double *restrict q,
                                     double *restrict velocity,
                                     double *restrict jacob, unsigned int sizeZ,
                                     double *restrict z) {
  // Jacobian of internal forces with respect to velocity (columnwise)
  jacob[0] = 0.;
  jacob[1] = 0.;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void g1(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict g,
                       unsigned int sizeZ, double *restrict z) {
  g[0] = 0.5 * parameters::d -
         (parameters::l1 * sin(q[0]) + parameters::l2 * sin(q[1]) -
          parameters::a * sin(q[2]) + parameters::b * cos(q[2])); // normal
  // printf("a = %e\t", a);
  // printf("b = %e\t",b);
  // printf("c = %e\n", c);
  // printf("parameters::l1 = %e\t",parameters::l1);
  // printf("parameters::l2 = %e\n", parameters::l2);
  // printf("q[0] = %e\t", q[0]);
  // printf("q[1] = %e\t", q[1]);
  // printf("q[2] = %e\n", q[2]);
  // printf("g1[0] = %e\n", g[0]);

  if (sizeOfY > 1)
    g[1] = parameters::l1 * cos(q[0]) + parameters::l2 * cos(q[1]) -
           parameters::a * cos(q[2]) - parameters::b * sin(q[2]); // tangential
}

SICONOS_EXPORT void W1(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict W,
                       unsigned int sizeZ, double *restrict z) {
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1) {
    W[0] = -parameters::l1 * cos(q[0]);
    W[1] = -parameters::l2 * cos(q[1]);
    W[2] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  } else if (sizeOfY == 2) {
    W[0] = -parameters::l1 * cos(q[0]);
    W[1] = -parameters::l1 * sin(q[0]);

    W[2] = -parameters::l2 * cos(q[1]);
    W[3] = -parameters::l2 * sin(q[1]);

    W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
    W[5] = parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  } else
    THROW_EXCEPTION("W1 - not implemented!");
}

SICONOS_EXPORT void W1dot(unsigned int sizeOfq, double *restrict q,
                          unsigned int sizeOfqdot, double *restrict qdot,
                          double *restrict S2, unsigned int sizeOfZ,
                          double *restrict z) {
  int sizeOfY = 1;
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1) {
    S2[0] = parameters::l1 * sin(q[0]) * qdot[0];
    S2[1] = parameters::l2 * sin(q[1]) * qdot[1];
    S2[2] = (-parameters::a * sin(q[2]) + parameters::b * cos(q[2])) * qdot[2];
  }
  //   // else if (sizeOfY == 2)
  //   // {
  //   //   W[0] = -parameters::l1 * cos(q[0]);
  //   //   W[1] = -parameters::l1 * sin(q[0]);

  //   //   W[2] = -parameters::l2 * cos(q[1]);
  //   //   W[3] = -parameters::l2 * sin(q[1]);

  //   //   W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  //   //   W[5] = parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  //   }
  //   else
  //     THROW_EXCEPTION("W1 - not implemented!");
}

SICONOS_EXPORT void g2(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict g,
                       unsigned int sizeZ, double *restrict z) {
  g[0] = 0.5 * parameters::d -
         (parameters::l1 * sin(q[0]) + parameters::l2 * sin(q[1]) +
          parameters::a * sin(q[2]) + parameters::b * cos(q[2])); // normal
  if (sizeOfY > 1)
    g[1] = parameters::l1 * cos(q[0]) + parameters::l2 * cos(q[1]) +
           parameters::a * cos(q[2]) - parameters::b * sin(q[2]); // tangential
}

SICONOS_EXPORT void W2(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict W,
                       unsigned int sizeZ, double *restrict z) {
  // Jacobian of g2 (columnwise)
  if (sizeOfY == 1) {
    W[0] = -parameters::l1 * cos(q[0]);
    W[1] = -parameters::l2 * cos(q[1]);
    W[2] = -parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  } else if (sizeOfY == 2) {
    W[0] = -parameters::l1 * cos(q[0]);
    W[1] = -parameters::l1 * sin(q[0]);

    W[2] = -parameters::l2 * cos(q[1]);
    W[3] = -parameters::l2 * sin(q[1]);

    W[4] = -parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
    W[5] = -parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  } else
    THROW_EXCEPTION("W2 - not implemented!");
}
SICONOS_EXPORT void W2dot(unsigned int sizeOfq, double *restrict q,
                          unsigned int sizeOfqdot, double *restrict qdot,
                          double *restrict S2, unsigned int sizeOfZ,
                          double *restrict z) {
  int sizeOfY = 1;
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1) {
    S2[0] = parameters::l1 * sin(q[0]) * qdot[0];
    S2[1] = parameters::l2 * sin(q[1]) * qdot[1];
    S2[2] = (-parameters::a * sin(q[2]) + parameters::b * cos(q[2])) * qdot[2];
  }
  //   // else if (sizeOfY == 2)
  //   // {
  //   //   W[0] = -parameters::l1 * cos(q[0]);
  //   //   W[1] = -parameters::l1 * sin(q[0]);

  //   //   W[2] = -parameters::l2 * cos(q[1]);
  //   //   W[3] = -parameters::l2 * sin(q[1]);

  //   //   W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  //   //   W[5] = parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  //   }
  //   else
  //     THROW_EXCEPTION("W1 - not implemented!");
}

SICONOS_EXPORT void g3(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict g,
                       unsigned int sizeZ, double *restrict z) {
  g[0] = 0.5 * parameters::d + parameters::l1 * sin(q[0]) +
         parameters::l2 * sin(q[1]) - parameters::a * sin(q[2]) -
         parameters::b * cos(q[2]); // normal
  if (sizeOfY > 1)
    g[1] = parameters::l1 * cos(q[0]) + parameters::l2 * cos(q[1]) -
           parameters::a * cos(q[2]) + parameters::b * sin(q[2]); // tangential
}

SICONOS_EXPORT void W3(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict W,
                       unsigned int sizeZ, double *restrict z) {
  // Jacobian of g3 (columnwise)
  if (sizeOfY == 1) {
    W[0] = parameters::l1 * cos(q[0]);
    W[1] = parameters::l2 * cos(q[1]);
    W[2] = -parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  } else if (sizeOfY == 2) {
    W[0] = parameters::l1 * cos(q[0]);
    W[1] = -parameters::l1 * sin(q[0]);

    W[2] = parameters::l2 * cos(q[1]);
    W[3] = -parameters::l2 * sin(q[1]);

    W[4] = -parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
    W[5] = parameters::a * sin(q[2]) + parameters::b * cos(q[2]);
  } else
    THROW_EXCEPTION("W3 - not implemented!");
}
SICONOS_EXPORT void W3dot(unsigned int sizeOfq, double *restrict q,
                          unsigned int sizeOfqdot, double *restrict qdot,
                          double *restrict S2, unsigned int sizeOfZ,
                          double *restrict z) {
  int sizeOfY = 1;
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1) {
    S2[0] = -parameters::l1 * sin(q[0]) * qdot[0];
    S2[1] = -parameters::l2 * sin(q[1]) * qdot[1];
    S2[2] = (parameters::a * sin(q[2]) + parameters::b * cos(q[2])) * qdot[2];
  }
  //   // else if (sizeOfY == 2)
  //   // {
  //   //   W[0] = -parameters::l1 * cos(q[0]);
  //   //   W[1] = -parameters::l1 * sin(q[0]);

  //   //   W[2] = -parameters::l2 * cos(q[1]);
  //   //   W[3] = -parameters::l2 * sin(q[1]);

  //   //   W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  //   //   W[5] = parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  //   }
  //   else
  //     THROW_EXCEPTION("W1 - not implemented!");
}
SICONOS_EXPORT void g4(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict g,
                       unsigned int sizeZ, double *restrict z) {
  g[0] = 0.5 * parameters::d + parameters::l1 * sin(q[0]) +
         parameters::l2 * sin(q[1]) + parameters::a * sin(q[2]) -
         parameters::b * cos(q[2]); // normal
  if (sizeOfY > 1)
    g[1] = parameters::l1 * cos(q[0]) + parameters::l2 * cos(q[1]) +
           parameters::a * cos(q[2]) + parameters::b * sin(q[2]); // tangential
}

SICONOS_EXPORT void W4(unsigned int sizeOfq, double *restrict q,
                       unsigned int sizeOfY, double *restrict W,
                       unsigned int sizeZ, double *restrict z) {
  // Jacobian of g4 (columnwise)
  if (sizeOfY == 1) {
    W[0] = parameters::l1 * cos(q[0]);
    W[1] = parameters::l2 * cos(q[1]);
    W[2] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  } else if (sizeOfY == 2) {
    W[0] = parameters::l1 * cos(q[0]);
    W[1] = -parameters::l1 * sin(q[0]);

    W[2] = parameters::l2 * cos(q[1]);
    W[3] = -parameters::l2 * sin(q[1]);

    W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
    W[5] = -parameters::a * sin(q[2]) + parameters::b * cos(q[2]);
  } else
    THROW_EXCEPTION("W4 - not implemented!");
}
SICONOS_EXPORT void W4dot(unsigned int sizeOfq, double *restrict q,
                          unsigned int sizeOfqdot, double *restrict qdot,
                          double *restrict S2, unsigned int sizeOfZ,
                          double *restrict z) {
  int sizeOfY = 1;
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1) {
    S2[0] = -parameters::l1 * sin(q[0]) * qdot[0];
    S2[1] = -parameters::l2 * sin(q[1]) * qdot[1];
    S2[2] = (parameters::a * sin(q[2]) + parameters::b * cos(q[2])) * qdot[2];
  }
  //   // else if (sizeOfY == 2)
  //   // {
  //   //   W[0] = -parameters::l1 * cos(q[0]);
  //   //   W[1] = -parameters::l1 * sin(q[0]);

  //   //   W[2] = -parameters::l2 * cos(q[1]);
  //   //   W[3] = -parameters::l2 * sin(q[1]);

  //   //   W[4] = parameters::a * cos(q[2]) + parameters::b * sin(q[2]);
  //   //   W[5] = parameters::a * sin(q[2]) - parameters::b * cos(q[2]);
  //   }
  //   else
  //     THROW_EXCEPTION("W1 - not implemented!");
}
