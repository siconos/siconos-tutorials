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

/* !\file SMCElectroPneumatic.cpp
  \brief Simulation of an Electropneumatic setup controlled with a Twisting algorithm
  O. Huber
  */

#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

// main program
int main(int argc, char* argv[]) {
  // User-defined parameters
  unsigned ndof = 4;         // Number of degrees of freedom of your system
  double t0 = 0.0;           // Starting time
  double T = 100.0;          // Total simulation time
  double h = 1.0e-3;         // Time step for simulation
  double hControl = 1.0e-1;  // Time step for control
  double y0 = 0.0;
  double v0 = 0.0;
  double pP0 = 501570.871800124;
  double pN0 = 501570.871800124;

  Vector param{3};
  param << 1e-2, 2. / 3., 10.;
  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:

  // Matrix declaration
  // For the DynamycalSystem
  Vector x0{ndof};
  x0 << pP0, pN0, v0, y0;

  // Dynamical Systems
  auto plant =
      std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(x0, siconos::algebra::copy_t);

  plant->setComputefVectorFunction(
      [&param](const Eigen::Ref<const siconos::algebra::SiconosVector>& x, double time,
               Eigen::Ref<siconos::algebra::MapVectorType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = 0.0054 * v;
        double x1 = 0.0045 * y;
        double x2 = 1.0 / (x1 + 0.000374193137473672);
        double x3 = 1.0 / (-x1 + 0.000374193137473672);
        result(0) = -pp * x0 * x2 +
                    133.641661725 * x2 *
                        (6.95536270702858e-32 * pp * pp * pp * pp * pp -
                         1.54181893025638e-25 * pp * pp * pp * pp +
                         1.21454472133664e-19 * pp * pp * pp - 4.47823017035973e-14 * pp * pp +
                         6.78962532904015e-9 * pp + 8.53342339766871e-5);
        result(1) = pn * x0 * x3 +
                    133.641661725 * x3 *
                        (6.95536270702858e-32 * pn * pn * pn * pn * pn -
                         1.54181893025638e-25 * pn * pn * pn * pn +
                         1.21454472133664e-19 * pn * pn * pn - 4.47823017035973e-14 * pn * pn +
                         6.78962532904015e-9 * pn + 8.53342339766871e-5);
        result(2) =
            -0.00132352941176471 * pn + 0.00132352941176471 * pp - 14.7058823529412 * v;
        result(3) = v;
      });

  plant->setComputeJacobianfOver_xFunction(
      [&param](const Eigen::Ref<const siconos::algebra::SiconosVector>& x, double time,
               Eigen::Ref<siconos::algebra::MapType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = 0.0045 * y;
        double x1 = x0 + 0.000374193137473672;
        double x2 = 1.0 / (x1);
        double x3 = 0.0054 * x2;
        double x4 = pp * pp;
        double x5 = pp * pp * pp;
        double x6 = pp * pp * pp * pp;
        double x7 = 2.43e-5 * v;
        double x8 = 1. / (x1 * x1);
        double x9 = -x0 + 0.000374193137473672;
        double x10 = 1.0 / (x9);
        double x11 = 0.0054 * x10;
        double x12 = pn * pn;
        double x13 = pn * pn * pn;
        double x14 = pn * pn * pn * pn;
        double x15 = 1. / (x9 * x9);
        result.setZero();
        result(0, 0) = -v * x3 + 133.641661725 * x2 *
                                     (-8.95646034071946e-14 * pp + 3.64363416400992e-19 * x4 -
                                      6.1672757210255e-25 * x5 + 3.47768135351429e-31 * x6 +
                                      6.78962532904015e-9);
        result(0, 2) = -pp * x3;
        result(0, 3) =
            pp * x7 * x8 -
            0.6013874777625 * x8 *
                (6.95536270702858e-32 * pp * pp * pp * pp * pp + 6.78962532904015e-9 * pp -
                 4.47823017035973e-14 * x4 + 1.21454472133664e-19 * x5 -
                 1.54181893025638e-25 * x6 + 8.53342339766871e-5);
        result(1, 1) = v * x11 + 133.641661725 * x10 *
                                     (-8.95646034071946e-14 * pn + 3.64363416400992e-19 * x12 -
                                      6.1672757210255e-25 * x13 + 3.47768135351429e-31 * x14 +
                                      6.78962532904015e-9);
        result(1, 2) = pn * x11;
        result(1, 3) =
            pn * x15 * x7 +
            0.6013874777625 * x15 *
                (6.95536270702858e-32 * pn * pn * pn * pn * pn + 6.78962532904015e-9 * pn -
                 4.47823017035973e-14 * x12 + 1.21454472133664e-19 * x13 -
                 1.54181893025638e-25 * x14 + 8.53342339766871e-5);
        result(2, 0) = 0.00132352941176471;
        result(2, 1) = -0.00132352941176471;
        result(2, 2) = -14.7058823529412;
        result(3, 2) = 1;
      });

  // -------------
  // --- Model process ---
  // -------------
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(t0, T, h);

  simLsodar->addDynamicalSystem(plant);
  // ------------------
  // --- Simulation ---
  // ------------------
  // TimeDiscretisation
  // Control stuff
  // For the Sensor
  auto sensorC = std::make_shared<Matrix>(4, 4);
  sensorC->setIdentity();
  auto sensor = std::make_shared<siconos::control::LinearSensor>(plant, sensorC);
  // add the sliding mode controller
  auto twisting = std::make_shared<siconos::control::Twisting>(sensor, hControl);

  twisting->setComputehFunction(
      [&param](const siconos::algebra::BlockVector& x, double t,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
               Eigen::Ref<siconos::algebra::MapVectorType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double l1 = lamb(0);
        double l2 = lamb(1);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = 0.1 * t;
        double x1 = v - 0.004 * cos(x0);
        double x2 = sin(x0);
        result(0) = alpha * (-0.04 * x2 + y) + x1;
        result(1) = alpha * x1 - 0.00132352941176471 * pn + 0.00132352941176471 * pp -
                    14.7058823529412 * v + 0.0004 * x2;
      });

  twisting->setComputegFunction(
      [&param](const siconos::algebra::BlockVector& x, double t,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
               siconos::algebra::BlockVector& result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double l1 = lamb(0);
        double l2 = lamb(1);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = beta * l2 + l1;
        double x1 = 0.0045 * y;
        double x2 = 133.641661725 * G * x0 / (x1 + 0.000374193137473672);
        double x3 = pp * pp;
        double x4 = pp * pp * pp;
        double x5 = pp * pp * pp * pp;
        double x6 = pp * pp * pp * pp * pp;
        double x7 = G * x0;
        int x8 = x7 >= 0.0;
        int x9 = x7 < 0.0;
        double x10 = 133.641661725 * G * x0 / (-x1 + 0.000374193137473672);
        double x11 = pn * pn;
        double x12 = pn * pn * pn;
        double x13 = pn * pn * pn * pn;
        double x14 = pn * pn * pn * pn * pn;
        if (x8) {
          result(0) = x2 * (-6.27567493976828e-5 * pp + 4.32959552123787e-10 * x3 -
                            1.36184599900831e-15 * x4 + 1.96424682992079e-21 * x5 -
                            1.10088063828155e-27 * x6 + 14.462054466034);
          result(1) = -x10 * (0.000105066798974752 * pn - 4.59203872247267e-10 * x11 +
                              1.14855459001783e-15 * x12 - 1.38377833290208e-21 * x13 +
                              6.42465433109007e-28 * x14 - 6.59951125132568);
        };
        if (x9) {
          result(0) = x2 * (0.000105066798974752 * pp - 4.59203872247267e-10 * x3 +
                            1.14855459001783e-15 * x4 - 1.38377833290208e-21 * x5 +
                            6.42465433109007e-28 * x6 - 6.59951125132568);

          result(1) = -x10 * (-6.27567493976828e-5 * pn + 4.32959552123787e-10 * x11 -
                              1.36184599900831e-15 * x12 + 1.96424682992079e-21 * x13 -
                              1.10088063828155e-27 * x14 + 14.462054466034);
        };
        result(2) = 0.;
        result(3) = 0.;
      });

  twisting->setComputeJacobianhOver_stateFunction(
      [&param](const siconos::algebra::BlockVector& x, double time,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
               Eigen::Ref<siconos::algebra::MapType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double l1 = lamb(0);
        double l2 = lamb(1);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        result.setZero();
        result(0, 2) = 1.;
        result(0, 3) = alpha;
        result(1, 0) = 0.00132352941176471;
        result(1, 1) = -0.00132352941176471;
        result(1, 2) = alpha - 14.7058823529412;
      });

  twisting->setComputeJacobiangOver_stateFunction(
      [&param](const siconos::algebra::BlockVector& x, double t,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
               Eigen::Ref<siconos::algebra::MapType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double l1 = lamb(0);
        double l2 = lamb(1);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = beta * l2 + l1;
        double x1 = 0.0045 * y;
        double x2 = x1 + 0.000374193137473672;
        double x3 = 133.641661725 * G * x0 / x2;
        double x4 = pp * pp;
        double x5 = pp * pp * pp;
        double x6 = pp * pp * pp * pp;
        double x7 = G * x0;
        int x8 = x7 >= 0.0;
        int x9 = x7 < 0.0;
        double x10 = 0.6013874777625 * G * x0 / x2 * x2;
        double x11 = pp * pp * pp * pp * pp;
        double x12 = -x1 + 0.000374193137473672;
        double x13 = 133.641661725 * G * x0 / x12;
        double x14 = pn * pn;
        double x15 = pn * pn * pn;
        double x16 = pn * pn * pn * pn;
        double x17 = 0.6013874777625 * G * x0 / x12 * x12;
        double x18 = pn * pn * pn * pn * pn;
        result.setZero();
        if (x8) {
          result(0, 0) = x3 * (8.65919104247574e-10 * pp - 4.08553799702493e-15 * x4 +
                               7.85698731968318e-21 * x5 - 5.50440319140777e-27 * x6 -
                               6.27567493976828e-5);

          result(0, 3) = -x10 * (-6.27567493976828e-5 * pp - 1.10088063828155e-27 * x11 +
                                 4.32959552123787e-10 * x4 - 1.36184599900831e-15 * x5 +
                                 1.96424682992079e-21 * x6 + 14.462054466034);

          result(1, 1) = -x13 * (-9.18407744494533e-10 * pn + 3.44566377005348e-15 * x14 -
                                 5.5351133316083e-21 * x15 + 3.21232716554503e-27 * x16 +
                                 0.000105066798974752);

          result(1, 3) = -x17 * (0.000105066798974752 * pn - 4.59203872247267e-10 * x14 +
                                 1.14855459001783e-15 * x15 - 1.38377833290208e-21 * x16 +
                                 6.42465433109007e-28 * x18 - 6.59951125132568);
        };
        if (x9) {
          result(0, 0) = x3 * (-9.18407744494533e-10 * pp + 3.44566377005348e-15 * x4 -
                               5.5351133316083e-21 * x5 + 3.21232716554503e-27 * x6 +
                               0.000105066798974752);
          result(0, 3) = -x10 * (0.000105066798974752 * pp + 6.42465433109007e-28 * x11 -
                                 4.59203872247267e-10 * x4 + 1.14855459001783e-15 * x5 -
                                 1.38377833290208e-21 * x6 - 6.59951125132568);
          result(1, 1) = -x13 * (8.65919104247574e-10 * pn - 4.08553799702493e-15 * x14 +
                                 7.85698731968318e-21 * x15 - 5.50440319140777e-27 * x16 -
                                 6.27567493976828e-5);
          result(1, 3) = -x17 * (-6.27567493976828e-5 * pn + 4.32959552123787e-10 * x14 -
                                 1.36184599900831e-15 * x15 + 1.96424682992079e-21 * x16 -
                                 1.10088063828155e-27 * x18 + 14.462054466034);
        };
      });

  twisting->setComputeJacobiangOver_lambdaFunction(
      [&param](const siconos::algebra::BlockVector& x, double time,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
               Eigen::Ref<siconos::algebra::MapType> result) {
        double pp = x(0);
        double pn = x(1);
        double v = x(2);
        double y = x(3);
        double l1 = lamb(0);
        double l2 = lamb(1);
        double G = param(0);
        double beta = param(1);
        double alpha = param(2);
        double x0 = 0.0045 * y;
        double x1 = 133.641661725 * G / (x0 + 0.000374193137473672);
        double x2 = pp * pp;
        double x3 = pp * pp * pp;
        double x4 = pp * pp * pp * pp;
        double x5 = pp * pp * pp * pp * pp;
        double x6 = x1 * (-6.27567493976828e-5 * pp + 4.32959552123787e-10 * x2 -
                          1.36184599900831e-15 * x3 + 1.96424682992079e-21 * x4 -
                          1.10088063828155e-27 * x5 + 14.462054466034);
        double x7 = G * (beta * l2 + l1);
        int x8 = x7 >= 0.0;
        double x9 = x1 * (0.000105066798974752 * pp - 4.59203872247267e-10 * x2 +
                          1.14855459001783e-15 * x3 - 1.38377833290208e-21 * x4 +
                          6.42465433109007e-28 * x5 - 6.59951125132568);
        int x10 = x7 < 0.0;
        double x11 = 133.641661725 * G / (-x0 + 0.000374193137473672);
        double x12 = pn * pn;
        double x13 = pn * pn * pn;
        double x14 = pn * pn * pn * pn;
        double x15 = pn * pn * pn * pn * pn;
        double x16 = x11 * (0.000105066798974752 * pn - 4.59203872247267e-10 * x12 +
                            1.14855459001783e-15 * x13 - 1.38377833290208e-21 * x14 +
                            6.42465433109007e-28 * x15 - 6.59951125132568);
        double x17 = x11 * (-6.27567493976828e-5 * pn + 4.32959552123787e-10 * x12 -
                            1.36184599900831e-15 * x13 + 1.96424682992079e-21 * x14 -
                            1.10088063828155e-27 * x15 + 14.462054466034);
        result.setZero();
        if (x8) {
          result(0, 0) = x6;
          result(0, 1) = beta * x6;

          result(1, 0) = -x16;
          result(1, 1) = -beta * x16;
        };
        if (x10) {
          result(0, 0) = x9;
          result(0, 1) = beta * x9;
          result(1, 0) = -x17;
          result(1, 1) = -beta * x17;
        };
      });

  simLsodar->addSensor(sensor, hControl);
  simLsodar->addActuator(twisting, hControl);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ...\n\n";
  // initialise the process and the ControlManager
  simLsodar->initialize();

  //  (std::static_pointer_cast<LsodarOSI>(simLsodar->integrator()))->setJT(1);
  //  (std::static_pointer_cast<LsodarOSI>(simLsodar->integrator()))->setMaxOrder(0, 5);
  cout << "====> Simulation run ...\n\n";
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("SMCElectroPneumaticItw.dat", data,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-8;
  if ((error = siconos::algebra::io::compareRefFile(data, "SMCElectroPneumaticItw.ref", eps)) >
      eps)
    return 1;
  else
    return 0;
}
