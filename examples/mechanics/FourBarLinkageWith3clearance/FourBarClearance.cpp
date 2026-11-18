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

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>  // for pi

#include "tools.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

// parameters according to Table 1
// geometrical characteristics

using namespace user;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 7;     // degrees of freedom for robot arm
    double t0 = 0;    // initial computation time
    double T = 1.0;   // final computation time
    double h = 1e-4;  // time step
    // T = 2*h;
    // double criterion = 1e-6;
    // unsigned int maxIter = 2000;
    char filename[50] = "simu_";
    double eN = 0.0;
    // eN1 = 0.1;
    double eT = 0.0;
    double mu = 0.1;

    double r1, r3, r5, Kp, lmd;

    std::cout << "argc :" << argc << "\n";
    if (argc < 2) {
      std::cout << "Using default arguments\n";

      r1 = 0.055;
      r3 = 0.055;
      r5 = 0.055;
      Kp = 200;
      lmd = 10;
      int sizeofargv1 = strlen("0.055");
      int sizeofargv2 = strlen("0.055");
      int sizeofargv3 = strlen("0.055");
      int sizeofargv4 = strlen("200");
      int sizeofargv5 = strlen("10");
      strncpy(&filename[5], "0.055", sizeofargv1);
      strncpy(&filename[11], "0.055", sizeofargv2);
      strncpy(&filename[17], "0.055", sizeofargv3);
      strncpy(&filename[23], "200", sizeofargv4);
      strncpy(&filename[26], "10", sizeofargv5);

    } else if (argc == 6) {
      r1 = atof(argv[1]);
      r3 = atof(argv[2]);
      r5 = atof(argv[3]);
      Kp = atof(argv[4]);
      lmd = atof(argv[5]);
      int sizeofargv1 = strlen(argv[1]);
      int sizeofargv2 = strlen(argv[2]);
      int sizeofargv3 = strlen(argv[3]);
      int sizeofargv4 = strlen(argv[4]);
      int sizeofargv5 = strlen(argv[5]);
      strncpy(&filename[5], argv[1], sizeofargv1);
      strncpy(&filename[11], argv[2], sizeofargv2);
      strncpy(&filename[17], argv[3], sizeofargv3);
      strncpy(&filename[23], argv[4], sizeofargv4);
      strncpy(&filename[26], argv[5], sizeofargv5);
    } else {
      std::cout << "Wrong number of arguments\n";
      return 1;
    }

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: slidercrank ---

    // Initial position (angles in radian)
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    q0(0) = 1.570823772407980;   // 1.5708;
    q0(1) = 0.3532842020624460;  // 0.3533;
    q0(2) = 1.264872058968431;   // 1.2649;
    q0(3) = 1.876454585097650;   // 1.87647;
    q0(4) = 1.691962091335582;   // 1.69199;
    q0(5) = 0.3764686197082958;  // 0.3764+3.5e-5;
    q0(6) = 1.191962183453451;   // 1.19197;
    v0(0) = 0.0;
    v0(1) = 0.0;
    v0(2) = 0.0;

    auto fourbar =
        std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    fourbar->setComputeMassFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
           Eigen::Ref<siconos::algebra::MapType> mass) {
          mass.setZero();

          mass(0, 0) = J1 + (0.25 * m1) * l1 * l1;
          mass(1, 1) = J2;
          mass(2, 2) = J3;
          mass(3, 3) = m2;
          mass(4, 4) = m2;
          mass(5, 5) = m3;
          mass(6, 6) = m3;
        });

    // In the plugin functions, we use 'params' variable to save user-defined parameters, with:
    // params = [params0, params1, params2, r1, r3, r5, Kp, lmd]
    // params values are supposed to be initialized in main driver file.

    std::vector<double> params = {0.0, 0.0, 0.0, r1, r3, r5, Kp, lmd};
    fourbar->setComputeFintFunction(
        [&params](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                  const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
                  Eigen::Ref<siconos::algebra::MapVectorType> fint) {
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
          double mass11 = MASS1(q);
          double nonnl11 = NonNL1(q, velocity);
          double Kp = params[6];
          double lmd = params[7];
          fint.setZero();
          fint(0) =
              (0.5 * m1) * gravity * l1 * cos(q(0)) -
              (2 * (Jx1 + (Jx2 * s11 * s11) + (Jx3 * s21 * s21) + (P1 * c11 * s11)) *
               (-6.0 * 0.75 * 0.75 * std::numbers::pi * std::numbers::pi *
                    sin(0.75 * std::numbers::pi * time) -
                lmd * (velocity(0) -
                       0.75 * std::numbers::pi * 6.0 * cos(0.75 * std::numbers::pi * time)))) -
              ((2 * Jx2 * s11 * T11 + 2 * J3 * s21 * T13 +
                P1 * (c11 * T11 + s11 * (T7 + s11 * T8))) *
               velocity(0)) *
                  (velocity(0) - lmd * (q(0) - 6.0 * sin(0.75 * std::numbers::pi * time))) +
              Kp * (velocity(0) -
                    0.75 * std::numbers::pi * 6.0 * cos(0.75 * std::numbers::pi * time)) +
              Kp * lmd * (q(0) - 6.0 * sin(0.75 * std::numbers::pi * time)) -
              (-gt1 - s11 * ga1 - s21 * gp1);

          fint(4) = m2 * gravity;
          fint(5) = 0.0;
          fint(6) = m3 * gravity;
          params[0] = mass11;
          params[1] = nonnl11;
          params[2] =
              (2 * (Jx1 + (Jx2 * s11 * s11) + (Jx3 * s21 * s21) + (P1 * c11 * s11)) *
               (-6.0 * 0.75 * 0.75 * std::numbers::pi * std::numbers::pi *
                    sin(0.75 * std::numbers::pi * time) -
                lmd * (velocity(0) -
                       0.75 * std::numbers::pi * 6.0 * cos(0.75 * std::numbers::pi * time)))) +
              ((2 * Jx2 * s11 * T11 + 2 * J3 * s21 * T13 +
                P1 * (c11 * T11 + s11 * (T7 + s11 * T8))) *
               velocity(0)) *
                  (velocity(0) - lmd * (q(0) - 6.0 * sin(0.75 * std::numbers::pi * time))) -
              Kp * (velocity(0) -
                    0.75 * std::numbers::pi * 6.0 * cos(0.75 * std::numbers::pi * time)) -
              Kp * lmd * (q(0) - 6.0 * sin(0.75 * std::numbers::pi * time)) +
              (-gt1 - s11 * ga1 - s21 * gp1);
        });

    // fourbar->setComputeFgyrFunction(
    //     [](const Eigen::Ref<const siconos::algebra::SiconosVector>& vel,
    //        const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
    //        Eigen::Ref<siconos::algebra::MapVectorType> result) { result.setZero(); });

    // fourbar->->setComputeJacobianFgyrOver_qFunction(
    //   [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
    //      const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
    //      Eigen::Ref<siconos::algebra::MapType> result) { result.setZero(); });

    fourbar->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& v,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) = -(0.5 * m1) * gravity * l1 * sin(q(0));
        });

    // -------------------
    // --- Interactions---
    // -------------------

    // InteractionsSet allInteractions;

    // -- relations --

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        eN, eT, mu, 2);  // EqualityConditionNSL NewtonImpactNSL
    auto relation = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation->setComputehFunction([params](const siconos::algebra::BlockVector& q,
                                           Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y.setZero();
      double v1 = fcnExpression1(*q.vector(0));
      double r1 = params[3];
      y(0) = l1 * (r2 - r1) - v1;
    });

    relation->setComputeHFunction(
        [params](const siconos::algebra::BlockVector& q,
                 Eigen::Ref<siconos::algebra::MapType> result) {
          double v1 = fcnExpression1(*q.vector(0));
          double r1 = params[3];
          result.setZero();
          result(0, 0) = (-q(3) * l1 * sin(q(0)) + 0.5 * l1 * l2 * sin(q(0) - q(1)) +
                          q(4) * l1 * cos(q(0))) /
                         v1;
          result(1, 0) = ((q(3) * l1 * cos(q(0)) - 0.5 * l1 * l2 * cos(q(0) - q(1)) +
                           q(4) * l1 * sin(q(0)) - l1 * l1) /
                          v1) -
                         r1;

          result(0, 1) = (-0.5 * q(3) * l2 * sin(q(1)) - 0.5 * l1 * l2 * sin(q(0) - q(1)) +
                          0.5 * q(4) * l2 * cos(q(1))) /
                         v1;
          result(1, 1) = ((0.5 * q(3) * l2 * cos(q(1)) - 0.5 * l1 * l2 * cos(q(0) - q(1)) +
                           0.5 * q(4) * l2 * sin(q(1)) - 0.25 * l2 * l2) /
                          v1) +
                         r2;

          result(0, 3) = (-q(3) + l1 * cos(q(0)) + 0.5 * l2 * cos(q(1))) / v1;
          result(1, 3) = (q(4) - 0.5 * l2 * sin(q(1)) - l1 * sin(q(0))) / v1;

          result(0, 4) = (-q(4) + l1 * sin(q(0)) + 0.5 * l2 * sin(q(1))) / v1;
          result(1, 4) = -(q(3) - 0.5 * l2 * cos(q(1)) - l1 * cos(q(0))) / v1;
        });

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        eN, eT, mu, 2);  // EqualityConditionNSL NewtonImpactNSL
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

    relation1->setComputehFunction([params](const siconos::algebra::BlockVector& q,
                                            Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y.setZero();
      double v2 = fcnExpression2(*q.vector(0));
      double r3 = params[4];
      y(0) = (r4 - r3) - v2;
    });

    relation1->setComputeHFunction(
        [params](const siconos::algebra::BlockVector& q,
                 Eigen::Ref<siconos::algebra::MapType> result) {
          double v2 = fcnExpression2(*q.vector(0));
          double r3 = params[4];

          result.setZero();
          result(1, 0) = 0.0;

          result(0, 1) = (-0.5 * l0 * l2 * sin(q(1)) - 0.5 * q(5) * l2 * sin(q(1)) +
                          0.5 * q(3) * l2 * sin(q(1)) - 0.25 * l2 * l3 * sin(q(1) - q(2)) -
                          0.5 * q(4) * l2 * cos(q(1)) + 0.5 * q(6) * l2 * cos(q(1))) /
                         v2;
          result(1, 1) =
              ((0.5 * l0 * l2 * cos(q(1)) + 0.5 * q(5) * l2 * cos(q(1)) -
                0.5 * q(3) * l2 * cos(q(1)) + 0.25 * l2 * l3 * cos(q(1) - q(2)) -
                0.5 * q(4) * l2 * sin(q(1)) + 0.5 * q(6) * l2 * sin(q(1)) - 0.25 * l2 * l2) /
               v2) -
              r3;

          result(0, 2) = (0.5 * l0 * l3 * sin(q(2)) + 0.5 * q(5) * l3 * sin(q(2)) -
                          0.5 * q(3) * l3 * sin(q(2)) + 0.25 * l2 * l3 * sin(q(1) - q(2)) +
                          0.5 * q(4) * l3 * cos(q(2)) - 0.5 * q(6) * l3 * cos(q(2))) /
                         v2;
          result(1, 2) =
              ((-0.5 * l0 * l3 * cos(q(2)) - 0.5 * q(5) * l3 * cos(q(2)) +
                0.5 * q(3) * l3 * cos(q(2)) + 0.25 * l2 * l3 * cos(q(1) - q(2)) +
                0.5 * q(4) * l3 * sin(q(2)) - 0.5 * q(6) * l3 * sin(q(2)) - 0.25 * l3 * l3) /
               v2) +
              r4;

          result(0, 3) =
              (-q(3) + l0 + q(5) + 0.5 * l3 * cos(q(2)) - 0.5 * l2 * cos(q(1))) / v2;
          result(1, 3) = (q(4) - q(6) - 0.5 * l3 * sin(q(2)) + 0.5 * l2 * sin(q(1))) / v2;

          result(0, 4) = (-q(4) + q(6) + 0.5 * l3 * sin(q(2)) - 0.5 * l2 * sin(q(1))) / v2;
          result(1, 4) = (q(5) + l0 - q(3) + 0.5 * l3 * cos(q(2)) - 0.5 * l2 * cos(q(1))) / v2;

          result(0, 5) =
              (-q(5) - l0 - 0.5 * l3 * cos(q(2)) + 0.5 * l2 * cos(q(1)) + q(3)) / v2;
          result(1, 5) = (q(6) - q(4) + 0.5 * l3 * sin(q(2)) - 0.5 * l2 * sin(q(1))) / v2;

          result(0, 6) = (-q(6) + q(4) - 0.5 * l3 * sin(q(2)) + 0.5 * l2 * sin(q(1))) / v2;
          result(1, 6) =
              (-q(5) - l0 - 0.5 * l3 * cos(q(2)) + 0.5 * l2 * cos(q(1)) + q(3)) / v2;
        });
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        eN, eT, mu, 2);  // EqualityConditionNSL NewtonImpactNSL
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

    relation2->setComputehFunction([params](const siconos::algebra::BlockVector& q,
                                            Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y.setZero();
      double v3 = fcnExpression3(*q.vector(0));
      double r5 = params[5];
      y(0) = (r6 - r5) - v3;
    });

    relation2->setComputeHFunction(
        [params](const siconos::algebra::BlockVector& q,
                 Eigen::Ref<siconos::algebra::MapType> result) {
          double v3 = fcnExpression3(*q.vector(0));
          double r5 = params[5];
          result.setZero();
          result(0, 2) = (-0.5 * q(5) * l3 * sin(q(2)) + 0.5 * q(6) * l3 * cos(q(2))) / v3;
          result(1, 2) =
              ((0.5 * q(5) * l3 * cos(q(2)) + 0.5 * q(6) * l3 * sin(q(2)) - 0.25 * l3 * l3) /
               v3) -
              r5;

          result(0, 5) = (-q(5) + 0.5 * l3 * cos(q(2))) / v3;
          result(1, 5) = (q(6) - 0.5 * l3 * sin(q(2))) / v3;

          result(0, 6) = (-q(6) + 0.5 * l3 * sin(q(2))) / v3;
          result(1, 6) = (-q(5) + 0.5 * l3 * cos(q(2))) / v3;
        });

    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // auto nslaw2= std::make_shared<siconos::modeling::EqualityConditionNSL>(e1);
    // //EqualityConditionNSL NewtonImpactNSL auto relation2=
    // std::make_shared<siconos::modeling::LagrangianScleronomousR>("FourBarClearancePlugin:g3",
    // "FourBarClearancePlugin:W3"); auto inter2=
    // std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // auto nslaw3= std::make_shared<siconos::modeling::EqualityConditionNSL>(e1);
    // //EqualityConditionNSL NewtonImpactNSL auto relation3=
    // std::make_shared<siconos::modeling::LagrangianScleronomousR>("FourBarClearancePlugin:g4",
    // "FourBarClearancePlugin:W4"); auto inter3=
    // std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // allInteractions.insert(inter);
    // allInteractions.insert(inter1);
    // allInteractions.insert(inter2);
    // allInteractions.insert(inter3);

    // -------------
    // --- Model ---
    // -------------

    auto FourBarIdeal = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    FourBarIdeal->insertDynamicalSystem(fourbar);
    FourBarIdeal->link(inter, fourbar);
    FourBarIdeal->link(inter1, fourbar);
    FourBarIdeal->link(inter2, fourbar);
    // FourBarIdeal->link(inter3, fourbar);
    //  ----------------
    //  --- Simulation ---
    //  ----------------
    fourbar->computeTotalForces(fourbar->velocity_read(), fourbar->q_read(), t0);

    inter->computeOutput(t0, 0);
    inter1->computeOutput(t0, 0);
    inter2->computeOutput(t0, 0);

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto OSI = std::make_shared<siconos::integrators::MoreauJeanCombinedProjectionOSI>(0.5);
    // -- set the integrator for the four bar linkage --

    auto impact = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(
        2, SICONOS_FRICTION_2D_NSGS); /*,SICONOS_FRICTION_2D_ENUM
        SICONOS_FRICTION_2D_LEMKE notworking //
        SICONOS_FRICTION_2D_PGS=>working 0.75PI */
    impact->numericsSolverOptions()->dparam[0] = 1e-6;
    impact->numericsSolverOptions()->iparam[0] = 2000;
    impact->numericsSolverOptions()->iparam[2] = 1;  // random
    auto position =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM);

    auto s = std::make_shared<siconos::simulation::TimeSteppingCombinedProjection>(
        FourBarIdeal, t, OSI, impact, position, 2);
    s->setProjectionMaxIteration(2000);
    s->setConstraintTolUnilateral(1e-6);
    s->setConstraintTol(1e-6);

    std::cout << "=== End of model loading === \n";
    /////////////////////////////////////////////////////////////////////////
    // -- Time discretisation --
    // auto t= std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // auto s= std::make_shared<siconos::simulation::TimeStepping>(t);

    // -- OneStepIntegrators --//

    // double theta = 0.500001;

    // auto OSI(new Moreau(fourbar, theta);
    // s->insertIntegrator(OSI);

    // -- OneStepNsProblem --//

    // auto osnspb=
    // std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);

    // s->insertNonSmoothProblem(osnspb);

    std::cout << "=== End of model loading === \n";

    ////////////////////////////////////////////////////////////////////////

    // =========================== End of model definition ===========================
    // dataPlot(k,7) = (*inter->y(0))(0);

    // ================================= Computation =================================

    int k = 1;
    int N = ceil((T - t0) / h) + 1;
    std::cout << "Number of time step   " << N << "\n";
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 42;
    Matrix dataPlot((N / 500) + 1, outputSize);
    Matrix beam1Plot(2, 3 * ((N / 500) + 1));
    Matrix beam2Plot(2, 3 * ((N / 500) + 1));
    Matrix beam3Plot(2, 3 * ((N / 500) + 1));
    Matrix beam4Plot(2, 3 * ((N / 500) + 1));
    Matrix beam5Plot(1, 2 * ((N / 500) + 1));
    Matrix beam6Plot(1, 1 * ((N / 500) + 1));
    Matrix beam7Plot(1, 2 * ((N / 500) + 1));
    Matrix beam8Plot(1, 2 * ((N / 500) + 1));
    Matrix beam9Plot(1, 2 * ((N / 500) + 1));
    Matrix beam10Plot(1, 2 * ((N / 500) + 1));
    Matrix beam11Plot(1, 2 * ((N / 500) + 1));
    Matrix beam12Plot(1, 2 * ((N / 500) + 1));
    Matrix beam13Plot(1, 2 * ((N / 500) + 1));
    // std::cout << "size here " << (N/20) + 1 << "\n";
    //  For the initial time step:
    //  time
    auto q = fourbar->q();
    auto v = fourbar->velocity();
    dataPlot(0, 0) = FourBarIdeal->t0();
    dataPlot(0, 1) = (*q)(0);  // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*q)(3);
    dataPlot(0, 5) = (*q)(4);
    dataPlot(0, 6) = (*q)(5);
    dataPlot(0, 7) = (*q)(6);
    dataPlot(0, 8) = (*v)(0);
    dataPlot(0, 9) = (*v)(1);
    dataPlot(0, 10) = (*v)(2);
    dataPlot(0, 11) = (*inter1->y(1))(0);
    dataPlot(0, 12) = (*inter1->y(1))(1);
    dataPlot(0, 13) = (*q)(3) - l1 * cos((*q)(0)) - 0.5 * l2 * cos((*q)(1));
    dataPlot(0, 14) = (*q)(4) - l1 * sin((*q)(0)) - 0.5 * l2 * sin((*q)(1));
    dataPlot(0, 15) =
        l0 + (*q)(5) + 0.5 * l3 * cos((*q)(2)) - 0.5 * l2 * cos((*q)(1)) - (*q)(3);
    dataPlot(0, 16) = (*q)(6) + 0.5 * l3 * sin((*q)(2)) - 0.5 * l2 * sin((*q)(1)) - (*q)(4);
    dataPlot(0, 17) = (*q)(5) - 0.5 * l3 * cos((*q)(2));
    dataPlot(0, 18) = (*q)(6) - 0.5 * l3 * sin((*q)(2));
    dataPlot(0, 19) = (*inter->y(1))(0);
    dataPlot(0, 20) = (*inter->y(1))(1);
    dataPlot(0, 21) = (*inter->lambda(1))(0);   // lambda1
    dataPlot(0, 22) = (*inter1->lambda(1))(0);  // lambda2
    dataPlot(0, 23) = (*inter2->y(1))(0);
    dataPlot(0, 24) = (*inter2->y(1))(1);
    dataPlot(0, 25) = (*inter2->lambda(1))(0);  // lambda2
    dataPlot(0, 26) = (*inter2->lambda(1))(1);  // lambda2
    dataPlot(0, 27) = (*inter1->lambda(1))(1);  // lambda2
    dataPlot(0, 28) = (*inter->lambda(1))(1);   // lambda2
    dataPlot(0, 29) = sqrt(pow(((*q)(3) - 0.5 * l2 * cos((*q)(1)) - l1 * cos((*q)(0))), 2) +
                           pow(((*q)(4) - 0.5 * l2 * sin((*q)(1)) - l1 * sin((*q)(0))), 2));
    dataPlot(0, 30) = sqrt(
        pow((l0 + (*q)(5) + 0.5 * l3 * cos((*q)(2)) - 0.5 * l2 * cos((*q)(1)) - (*q)(3)), 2) +
        pow((-(*q)(4) + (*q)(6) + 0.5 * l3 * sin((*q)(2)) - 0.5 * l2 * sin((*q)(1))), 2));
    dataPlot(0, 31) = sqrt((pow(((*q)(5) - 0.5 * l3 * cos((*q)(2))), 2) +
                            pow(((*q)(6) - 0.5 * l3 * sin((*q)(2))), 2)));
    dataPlot(0, 32) = (*inter->y(0))(0);
    dataPlot(0, 33) = (*inter1->y(0))(0);
    dataPlot(0, 34) = (*inter2->y(0))(0);

    dataPlot(0, 35) = fourbar->fint()(0) - (0.5 * m1) * gravity * l1 * cos((*q)(0));
    dataPlot(0, 36) =
        0.5 * ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h)) *
            ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h)) +
        0.5 * 3000.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * h)) *
            ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * h)) +
        0.5 * 10.0 *
            ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h)) *
            ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * h));
    dataPlot(0, 37) = 6.0 * sin(0.75 * std::numbers::pi * h);
    dataPlot(0, 38) = 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h);
    dataPlot(0, 39) = params[2];
    dataPlot(0, 40) =
        0.5 *
        (((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h)) +
         10.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * h))) *
        (((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * h)) +
         10.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * h)));
    dataPlot(0, 41) = 0.5 * ((*inter->y(1))(0) * (*inter->y(1))(0));
    auto start = std::chrono::system_clock::now();
    // --- Time loop ---
    std::cout << "Start computation ... \n";

    double tt = 0;
    int kk = 1;
    while (s->hasNextEvent()) {
      // k++;
      //
      //  if (!(div(k,1000).rem))  std::cout <<"Step number "<< k << "\n";
      s->advanceToEvent();
      // Solve problem
      // s->newtonSolve(criterion, maxIter); //2000000
      // Data Output

      tt = s->nextTime();
      if (((k % 500) == 0) && (k != 0)) {
        // std::cout << "size here " << kk << endl;
        dataPlot(kk, 0) = tt;
        dataPlot(kk, 1) = (*q)(0);  // crank revolution
        dataPlot(kk, 2) = (*q)(1);
        dataPlot(kk, 3) = (*q)(2);
        dataPlot(kk, 4) = (*q)(3);
        dataPlot(kk, 5) = (*q)(4);
        dataPlot(kk, 6) = (*q)(5);
        dataPlot(kk, 7) = (*q)(6);
        dataPlot(kk, 8) = (*v)(0);
        dataPlot(kk, 9) = (*v)(1);
        dataPlot(kk, 10) = (*v)(2);
        dataPlot(kk, 11) = (*inter1->y(1))(0);
        dataPlot(kk, 12) = (*inter1->y(1))(1);
        dataPlot(kk, 13) = (*q)(3) - l1 * cos((*q)(0)) - 0.5 * l2 * cos((*q)(1));
        dataPlot(kk, 14) = (*q)(4) - l1 * sin((*q)(0)) - 0.5 * l2 * sin((*q)(1));
        dataPlot(kk, 15) =
            l0 + (*q)(5) + 0.5 * l3 * cos((*q)(2)) - 0.5 * l2 * cos((*q)(1)) - (*q)(3);
        dataPlot(kk, 16) =
            (*q)(6) + 0.5 * l3 * sin((*q)(2)) - 0.5 * l2 * sin((*q)(1)) - (*q)(4);
        dataPlot(kk, 17) = (*q)(5) - 0.5 * l3 * cos((*q)(2));
        dataPlot(kk, 18) = (*q)(6) - 0.5 * l3 * sin((*q)(2));
        dataPlot(kk, 19) = (*inter->y(1))(0);
        dataPlot(kk, 20) = (*inter->y(1))(1);
        dataPlot(kk, 21) = (*inter->lambda(1))(0);   // lambda1
        dataPlot(kk, 22) = (*inter1->lambda(1))(0);  // lambda2
        dataPlot(kk, 23) = (*inter2->y(1))(0);
        dataPlot(kk, 24) = (*inter2->y(1))(1);
        dataPlot(kk, 25) = (*inter2->lambda(1))(0);  // lambda2
        dataPlot(kk, 26) = (*inter2->lambda(1))(1);  // lambda2
        dataPlot(kk, 27) = (*inter1->lambda(1))(1);  // lambda2
        dataPlot(kk, 28) = (*inter->lambda(1))(1);   // lambda2
        dataPlot(kk, 29) =
            sqrt(pow(((*q)(3) - 0.5 * l2 * cos((*q)(1)) - l1 * cos((*q)(0))), 2) +
                 pow(((*q)(4) - 0.5 * l2 * sin((*q)(1)) - l1 * sin((*q)(0))), 2));
        dataPlot(kk, 30) = sqrt(
            pow((l0 + (*q)(5) + 0.5 * l3 * cos((*q)(2)) - 0.5 * l2 * cos((*q)(1)) - (*q)(3)),
                2) +
            pow((-(*q)(4) + (*q)(6) + 0.5 * l3 * sin((*q)(2)) - 0.5 * l2 * sin((*q)(1))), 2));
        dataPlot(kk, 31) = sqrt((pow(((*q)(5) - 0.5 * l3 * cos((*q)(2))), 2) +
                                 pow(((*q)(6) - 0.5 * l3 * sin((*q)(2))), 2)));
        dataPlot(kk, 32) = (*inter->y(0))(0);
        dataPlot(kk, 33) = (*inter1->y(0))(0);
        dataPlot(kk, 34) = (*inter2->y(0))(0);
        dataPlot(kk, 35) = fourbar->fint()(0) - (0.5 * m1) * gravity * l1 * cos((*q)(0));
        dataPlot(kk, 36) =
            0.5 *
                ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt)) *
                ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt)) +
            0.5 * 3000.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * tt)) *
                ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * tt)) +
            0.5 * 10.0 *
                ((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt)) *
                ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * tt));
        dataPlot(kk, 37) = 6.0 * sin(0.75 * std::numbers::pi * tt);
        dataPlot(kk, 38) = 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt);
        dataPlot(kk, 39) = params[2];
        dataPlot(kk, 40) =
            0.5 *
            (((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt)) +
             10.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * tt))) *
            (((*v)(0) - 6.0 * 0.75 * std::numbers::pi * cos(0.75 * std::numbers::pi * tt)) +
             10.0 * ((*q)(0) - 6.0 * sin(0.75 * std::numbers::pi * tt)));
        dataPlot(kk, 41) = 0.5 * ((*inter->y(1))(0) * (*inter->y(1))(0));

        beam1Plot(0, 3 * kk) = 0.0;
        beam1Plot(0, 3 * kk + 1) = 0.0;
        beam1Plot(0, 3 * kk + 2) = 0.0;
        beam1Plot(1, 3 * kk) = l1 * cos((*q)(0));
        beam1Plot(1, 3 * kk + 1) = l1 * sin((*q)(0));
        beam1Plot(1, 3 * kk + 2) = 0.0;

        beam2Plot(0, 3 * kk) = ((*q)(3) - 0.5 * l2 * cos((*q)(1))) +
                               (-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
        beam2Plot(0, 3 * kk + 1) = ((*q)(4) - 0.5 * l2 * sin((*q)(1))) +
                                   (-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));
        beam2Plot(0, 3 * kk + 2) = 0.0;
        beam2Plot(1, 3 * kk) = (*q)(3) + 0.5 * l2 * cos((*q)(1));
        beam2Plot(1, 3 * kk + 1) = (*q)(4) + 0.5 * l2 * sin((*q)(1));
        beam2Plot(1, 3 * kk + 2) = 0.0;

        beam3Plot(0, 3 * kk) = l0 + l3 * cos((*q)(2));
        beam3Plot(0, 3 * kk + 1) = l3 * sin((*q)(2));
        beam3Plot(0, 3 * kk + 2) = 0.0;
        beam3Plot(1, 3 * kk) = l0;
        beam3Plot(1, 3 * kk + 1) = 0.0;
        beam3Plot(1, 3 * kk + 2) = 0.0;

        beam4Plot(0, 3 * kk) = 0.0;
        beam4Plot(0, 3 * kk + 1) = 0.0;
        beam4Plot(0, 3 * kk + 2) = 0.0;
        beam4Plot(1, 3 * kk) = l0;
        beam4Plot(1, 3 * kk + 1) = 0.0;
        beam4Plot(1, 3 * kk + 2) = 0.0;

        beam5Plot(0, 2 * kk) = l1 * cos((*q)(0));
        beam5Plot(0, 2 * kk + 1) = l1 * sin((*q)(0));
        // beam5Plot(0,3*kk+2) = 0.0;
        // beam5Plot(1,3*kk) = 0.0;
        // beam5Plot(1,3*kk+1) =0.0;
        // beam5Plot(1,3*kk+2) = 0.0;

        beam6Plot(0, 1 * kk) = tt;

        beam7Plot(0, 2 * kk) = ((*q)(3) - 0.5 * l2 * cos((*q)(1))) -
                               (-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
        beam7Plot(0, 2 * kk + 1) = ((*q)(4) - 0.5 * l2 * sin((*q)(1))) -
                                   (-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));

        beam8Plot(0, 2 * kk) =
            -2.0 - 10 * (-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
        beam8Plot(0, 2 * kk + 1) =
            2.5 - 10 * (-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));

        beam9Plot(0, 2 * kk) = -2.0;
        beam9Plot(0, 2 * kk + 1) = 2.5;

        beam10Plot(0, 2 * kk) = (l0 + l3 * cos((*q)(2)));
        beam10Plot(0, 2 * kk + 1) = (l3 * sin((*q)(2)));

        beam11Plot(0, 2 * kk) = ((*q)(3) + 0.5 * l2 * cos((*q)(1))) -
                                (-l0 - l3 * cos((*q)(2)) + 0.5 * l2 * cos((*q)(1)) + (*q)(3));
        beam11Plot(0, 2 * kk + 1) = ((*q)(4) + 0.5 * l2 * sin((*q)(1))) -
                                    ((*q)(4) - l3 * sin((*q)(2)) + 0.5 * l2 * sin((*q)(1)));

        beam12Plot(0, 2 * kk) =
            -2.0 + 10 * (-l0 - l3 * cos((*q)(2)) + 0.5 * l2 * cos((*q)(1)) + (*q)(3));
        beam12Plot(0, 2 * kk + 1) =
            -1.5 + 10 * ((*q)(4) - l3 * sin((*q)(2)) + 0.5 * l2 * sin((*q)(1)));

        beam13Plot(0, 2 * kk) = -2.0;
        beam13Plot(0, 2 * kk + 1) = -1.5;

        kk++;
      }
      // dataPlot(k, 17) = s->getNewtonNbSteps();
      // dataPlot(k, 18) = s->nbProjectionIteration();
      // dataPlot(k, 19) = s->maxViolationUnilateral();
      // dataPlot(k, 20) = s->nbIndexSetsIteration();
      // dataPlot(k, 21) = s->cumulatedNewtonNbSteps();
      // dataPlot(k, 22) = s->nbCumulatedProjectionIteration();
      s->processEvents();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write(filename, dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("Link1.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("Link2.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("Link3.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("Link4.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex_ey.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("time.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1_ey1.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1x_ey1x.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy_ey1xy.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy1_ey1xy1.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy2_ey1xy2.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy3_ey1xy3..dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy_ey1xy.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("ex1xy4_ey1xy4.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 2.e-05;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "FourBarClearance.ref", eps)) >
        eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
