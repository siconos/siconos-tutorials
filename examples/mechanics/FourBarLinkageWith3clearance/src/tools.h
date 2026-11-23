/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <SiconosVector.hpp>
#include <cmath>
#include <iostream>
#include <numbers>  // for pi

namespace user {

// parameters according to Table 1
// geometrical characteristics
constexpr double l1 = 1.0;
constexpr double l2 = 4.0;
constexpr double l3 = 2.5;
constexpr double l0 = 3.0;
constexpr double r2 = 0.06;
constexpr double gravity = 9.81;
constexpr double r4 = 0.06;
constexpr double r6 = 0.06;
// inertial properties
constexpr double m1 = 1.0;
constexpr double m2 = 1.0;
constexpr double m3 = 1.0;
constexpr double J1 = m1 * l1 * l1 / 12.0;
constexpr double I1 = m1 * l1 * l1 / 3.0;
constexpr double J2 = m2 * l2 * l2 / 12.0;
constexpr double I3 = m3 * l3 * l3 / 3.0;
constexpr double J3 = m3 * l3 * l3 / 12.0;

constexpr double Ix1 = m1 * l1 * l1 / 12.0;
constexpr double Ix2 = m2 * l2 * l2 / 12.0;
constexpr double Ix3 = m3 * l3 * l3 / 12.0;
constexpr double Jx1 = 0.5 * (0.25 * m1 * l1 * l1 + Ix1 + m2 * l1 * l1);
constexpr double Jx2 = 0.5 * (0.25 * m2 * l2 * l2 + Ix2);
constexpr double Jx3 = 0.5 * (0.25 * m3 * l3 * l3 + Ix3);
constexpr double P1 = 0.5 * m2 * l1 * l2;

double fcnExpression1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double result1;
  result1 = sqrt(pow((q(3) - 0.5 * l2 * cos(q(1)) - l1 * cos(q(0))), 2) +
                 pow((q(4) - 0.5 * l2 * sin(q(1)) - l1 * sin(q(0))), 2));
  return result1;
}

double fcnExpression2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double result2;
  result2 = sqrt(pow((l0 + q(5) + 0.5 * l3 * cos(q[2]) - 0.5 * l2 * cos(q(1)) - q(3)), 2) +
                 pow((-q(4) + q(6) + 0.5 * l3 * sin(q[2]) - 0.5 * l2 * sin(q(1))), 2));
  return result2;
}

double fcnExpression3(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double result3;
  result3 =
      sqrt((pow((q(5) - 0.5 * l3 * cos(q[2])), 2) + pow((q(6) - 0.5 * l3 * sin(q[2])), 2)));
  return result3;
}

double Fphi(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double phi, k1, k2, k3;
  k1 = -2.0 * l1 * l3 * sin(q(0));
  k2 = 2.0 * l3 * (l0 - l1 * cos(q(0)));
  k3 = l0 * l0 + l1 * l1 - l2 * l2 + l3 * l3 - 2 * l0 * l1 * cos(q(0));
  phi = 2.0 * (std::numbers::pi + atan2((-k1 - sqrt(k1 * k1 + k2 * k2 - k3 * k3)), (k3 - k2)));
  // cout << "phi" << phi << endl;
  return phi;
}

double Falfa(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double alfa;
  double v = Fphi(q);
  alfa = atan2(-l1 * sin(q(0)) + l3 * sin(v), l0 - l1 * cos(q(0)) + l3 * cos(v));
  // cout << "alfa " << alfa << endl;
  return alfa;
}

double FS1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double s1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  return (l1 * sin(v - q(0))) / (l2 * sin(v1 - v));
}

double FS2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double s2;
  double v = Fphi(q);
  double v1 = Falfa(q);
  s2 = (l1 * sin(v1 - q(0))) / (l3 * sin(v1 - v));
  return s2;
}

double Fc1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double c1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  c1 = cos(q(0) - v1);
  return c1;
}

double Fdtc1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double dtc1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  dtc1 = sin(v1 - q(0));
  return dtc1;
}

double Fdac1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double dac1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  dac1 = -sin(v1 - q(0));
  return dac1;
}

double Fdtg(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double dtg;
  double v = Fphi(q);
  double v1 = Falfa(q);
  dtg = -(0.5 * m1 * l1 + m2 * l1) * gravity * cos(q(0));
  return dtg;
}

double Fdag(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double dag;
  double v = Fphi(q);
  double v1 = Falfa(q);
  dag = -0.5 * m2 * l2 * gravity * cos(v1);
  return dag;
}

double Fdpg(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double dpg;
  double v = Fphi(q);
  double v1 = Falfa(q);
  dpg = -0.5 * m3 * l3 * gravity * cos(v);
  return dpg;
}

double Ft1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t1 = (-l1 * cos(v - q(0))) / (l2 * sin(v1 - v));
  // cout << "t1 " << t1 << endl;
  return t1;
}

double Fta1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta1 = (l1 * cos(v - q(0)) * cos(v1 - v)) / (l2 * sin(v1 - v) * sin(v1 - v));
  return ta1;
}

double Ftp1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp1;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp1 = (2.0 * l1 * cos(v - q(0))) / (l2 * cos(2.0 * v1 - 2.0 * v));
  return tp1;
}

double Ft2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t2;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t2 = (-l1 * sin(v - q(0)) * cos(v1 - v)) / (l2 * sin(v1 - v) * sin(v1 - v));
  // cout << "t2 " << t2 << endl;
  return t2;
}

double Ftt2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tt2;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tt2 = (l1 * cos(v - q(0)) * cos(v1 - v)) / (l2 * sin(v1 - v) * sin(v1 - v));
  return tt2;
}

double Fta2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta2;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta2 = (l1 * (-sin(v1 - v) * sin(v1 - v) + 2.0) * sin(v - q(0))) /
        (l2 * sin(v1 - v) * sin(v1 - v) * sin(v1 - v));
  return ta2;
}

double Ftp2(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp2;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp2 = (l1 * (sin(v1 - v) * sin(v1 - v) * sin(v - q(0)) -
               sin(v1 - v) * cos(v1 - v) * cos(v - q(0)) - 2.0 * sin(v - q(0)))) /
        (l2 * sin(v1 - v) * sin(v1 - v) * sin(v1 - v));
  return tp2;
}

double Ft3(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t3;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t3 = (-2.0 * l1 * sin(v1 - q(0))) / (l2 * cos(2.0 * v1 - 2.0 * v) - l2);
  // cout << "t3 " << t3 << endl;
  return t3;
}

double Ftt3(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tt3;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tt3 = (2.0 * l1 * cos(-v1 + q(0))) / (l2 * cos(-2.0 * v1 + 2.0 * v) - l2);
  return tt3;
}

double Fta3(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta3;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta3 = (l1 * (-4.0 * sin(-2.0 * v1 + 2.0 * v) * sin(-v1 + q(0)) -
               2.0 * cos(-2.0 * v1 + 2.0 * v) * cos(-v1 + q(0)) + 2.0 * cos(-v1 + q(0)))) /
        (l2 * (-sin(-2.0 * v1 + 2.0 * v) - 2.0 * cos(-2.0 * v1 + 2.0 * v) + 2.0));
  return ta3;
}

double Ftp3(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp3;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp3 = (4.0 * l1 * l2 * sin(-2.0 * v1 + 2.0 * v) * sin(-v1 + q(0))) /
        ((l2 * cos(-2.0 * v1 + 2.0 * v) - l2) * (l2 * cos(-2.0 * v1 + 2.0 * v) - l2));
  return tp3;
}

double Ft4(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t4;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t4 = (-l1 * cos(v1 - q(0))) / (l3 * sin(v1 - v));
  return t4;
}

double Ftt4(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tt4;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tt4 = (-l1 * sin(v1 - q(0))) / (l3 * sin(v1 - v));
  return tt4;
}

double Fta4(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta4;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta4 = (l1 * (-sin(-v1 + q(0)) * sin(v1 - v) + cos(-v1 + q(0)) * cos(v1 - v))) /
        (l3 * sin(v1 - v) * sin(v1 - v));
  return ta4;
}

double Ftp4(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp4;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp4 = (-l1 * cos(-v1 + q(0)) * cos(v1 - v)) / (l3 * sin(v1 - v) * sin(v1 - v));
  return tp4;
}

double Ft5(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t5;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t5 = (2.0 * l1 * sin(v - q(0))) / (l3 * cos(2.0 * v1 - 2.0 * v) - l3);
  return t5;
}

double Ftt5(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tt5;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tt5 = (-2.0 * l1 * cos(v - q(0))) / (l3 * cos(2.0 * v1 - 2.0 * v) - l3);
  return tt5;
}

double Fta5(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta5;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta5 = (-4.0 * l1 * sin(-2.0 * v1 + 2.0 * v) * sin(v - q(0))) /
        (l3 * (-sin(-2.0 * v1 + 2.0 * v) * sin(-2.0 * v1 + 2.0 * v) -
               2.0 * cos(-2.0 * v1 + 2.0 * v) + 2.0));
  return ta5;
}

double Ftp5(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp5;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp5 = (2.0 * l1 *
         (2.0 * sin(-2.0 * v1 + 2.0 * v) * sin(v - q(0)) +
          cos(-2.0 * v1 + 2.0 * v) * cos(v - q(0)) - cos(v - q(0)))) /
        (l3 * (-sin(-2.0 * v1 + 2.0 * v) * sin(-2.0 * v1 + 2.0 * v) -
               2.0 * cos(-2.0 * v1 + 2.0 * v) + 2.0));
  return tp5;
}

double Ft6(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double t6;
  double v = Fphi(q);
  double v1 = Falfa(q);
  t6 = (l1 * sin(v1 - q(0)) * cos(v1 - v)) / (l3 * sin(v1 - v) * sin(v1 - v));
  return t6;
}

double Ftt6(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tt6;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tt6 = (-l1 * cos(-v1 + q(0)) * cos(v1 - v)) / (l3 * sin(v1 - v) * sin(v1 - v));
  return tt6;
}

double Fta6(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ta6;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ta6 = (l1 * (-sin(-v1 + q(0)) * sin(v1 - v) * sin(v1 - v) + 2.0 * sin(-v1 + q(0)) +
               sin(v1 - v) * cos(-v1 + q(0)) * cos(v1 - v))) /
        (l3 * sin(v1 - v) * sin(v1 - v) * sin(v1 - v));
  return ta6;
}

double Ftp6(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double tp6;
  double v = Fphi(q);
  double v1 = Falfa(q);
  tp6 = (l1 * (sin(v1 - v) * sin(v1 - v) - 2.0) * sin(-v1 + q(0))) /
        (l3 * sin(v1 - v) * sin(v1 - v) * sin(v1 - v));
  return tp6;
}

double Fft9(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft9;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Ta1 = Fta1(q);
  double Tp1 = Ftp1(q);
  ft9 = -s11 + (Ta1 * s11) + (Tp1 * s21);
  return ft9;
}

double Fft10(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft10;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Tt2 = Ftt2(q);
  double Ta2 = Fta2(q);
  double Tp2 = Ftp2(q);
  ft10 = Tt2 + (Ta2 * s11) + (Tp2 * s21);
  return ft10;
}

double Fft11(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft11;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double T1 = Ft1(q);
  double T2 = Ft2(q);
  double T3 = Ft3(q);
  ft11 = T1 + (T2 * s11) + (T3 * s21);
  return ft11;
}

double Fft12(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft12;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Tt3 = Ftt3(q);
  double Ta3 = Fta3(q);
  double Tp3 = Ftp3(q);
  ft12 = Tt3 + (Ta3 * s11) + (Tp3 * s21);
  return ft12;
}

double Fft13(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft13;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double T4 = Ft4(q);
  double T5 = Ft5(q);
  double T6 = Ft6(q);
  ft13 = T4 + (T5 * s11) + (T6 * s21);
  return ft13;
}

double Fft14(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft14;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Tt4 = Ftt4(q);
  double Ta4 = Fta4(q);
  double Tp4 = Ftp4(q);
  ft14 = Tt4 + (Ta4 * s11) + (Tp4 * s21);
  return ft14;
}

double Fft15(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft15;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Tt5 = Ftt5(q);
  double Ta5 = Fta5(q);
  double Tp5 = Ftp5(q);
  ft15 = Tt5 + (Ta5 * s11) + (Tp5 * s21);
  return ft15;
}

double Fft16(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft16;
  double v = Fphi(q);
  double v1 = Falfa(q);
  double s11 = FS1(q);
  double s21 = FS2(q);
  double Tt6 = Ftt6(q);
  double Ta6 = Fta6(q);
  double Tp6 = Ftp6(q);
  ft16 = Tt6 + (Ta6 * s11) + (Tp6 * s21);
  return ft16;
}

double Fft17(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft17;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ft17 = -cos(v1 - q(0));
  return ft17;
}
double Fft18(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft18;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ft18 = cos(v1 - q(0));
  return ft18;
}

double Fft19(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft19;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ft19 = cos(v1 - q(0));
  return ft19;
}
double Fft20(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  double ft20;
  double v = Fphi(q);
  double v1 = Falfa(q);
  ft20 = -cos(v1 - q(0));
  return ft20;
}

// plugins for smooth equations of motion
double MASS1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q) {
  // columnwise definition of mass matrix
  double mass1;
  double s11 = FS1(q);
  double s21 = FS2(q);
  double c11 = Fc1(q);
  // double st1 = Fdts1(q);
  // double sa1 = Fdas1(q);
  // double sp1 = Fdps1(q);
  // double st2 = Fdts2(q);
  // double sa2 = Fdas2(q);
  // double sp2 = Fdps2(q);
  // double ct1 = Fdtc1(q);
  // double ca1 = Fdac1(q);
  // double gt1 = Fdtg(q);
  // double ga1 = Fdag(q);
  // double gp1 = Fdpg(q);

  mass1 = 2 * (J1 + (J2 * s11 * s11) + (J3 * s21 * s21) + (P1 * c11 * s11));
  return mass1;
  // cout << "mass" << mass[0] << endl;
}

double NonNL1(const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
              const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity) {
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  double nonnl1;
  double s11 = FS1(q);
  double s21 = FS2(q);
  double c11 = Fc1(q);
  double T7 = Fdtc1(q);
  double T8 = Fdac1(q);
  double T11 = Fft11(q);
  double T13 = Fft13(q);
  double ct1 = Fdtc1(q);
  double ca1 = Fdac1(q);
  double gt1 = Fdtg(q);
  double ga1 = Fdag(q);
  double gp1 = Fdpg(q);

  nonnl1 =
      (2 * J2 * s11 * T11 + 2 * J3 * s21 * T13 + P1 * (c11 * T11 + s11 * (T7 + s11 * T8))) *
      velocity(0) * velocity(0);
  return nonnl1;
  // NNL[0] = (2*J2*s11*(st1+s11*sa1+s21*sp1) + 2*J3*s21*(st2+s11*sa2+s21*sp2) +
  // P1*(c11*(st1+s11*sa1+s21*sp1) + s11*(T7 + s11*T8)))*velocity[0]*velocity[0]; cout << "NNL"
  // << NNL[0] << endl;
}
}  // namespace user
