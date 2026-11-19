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

// =============================== Robot arm sample (HuMAnsPa10)
// ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include <SiconosKernel.hpp>
#include <chrono>

#include "TwoLinksMultiPlugins.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace two_links_multi_plugins;  // Where plugin functions (lambdas) are defined

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    double t0 = 0;    // initial computation time
    double T = 30;    // final computation time
    double h = 1e-3;  // time step
    double criterion = 1e-8;
    unsigned int maxIter = 200000;
    double e = 0.7;  // nslaw
    double e2 = 0.0;
    // double L = 0.0;
    int test = 0;
    int nimpact = 0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: manipulator arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See
    // corresponding .pdf for more details)

    // Initial position (angles in radian)
    Vector q0(nDof), v0(nDof);
    q0.setZero();
    v0.setZero();
    q0(0) = 1.5;
    q0(1) = -0.9;

    param[0] = q0(0);
    param[1] = q0(1);
    param[2] = v0(0);
    param[3] = v0(1);
    param[11] = std::numbers::pi;

    auto arm =
        std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);

    arm->setComputeMassFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                                   Eigen::Ref<siconos::algebra::MapType> mass) {
      mass(0, 0) =
          m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q(1)));
      mass(1, 0) = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q(1)) / 2;
      mass(0, 1) = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q(1)) / 2;
      mass(1, 1) = I2 + m2 * l2 * l2 / 4;
    });

    arm->setComputeFgyrFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapVectorType> result) {
          result(0) = -m2 * l1 * l2 * sin(q(1)) *
                      (velocity(0) * velocity(1) + velocity(1) * velocity(1) / 2);
          result(1) = m2 * l1 * l2 * sin(q(1)) * velocity(0) * velocity(0) / 2;
        });

    // set 'random' value for jacobians, whatever fgyr is, just for tests
    arm->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = 0;
          result(1, 0) = 0;
          result(0, 1) = -m2 * l1 * l2 * cos(q(1)) *
                         (velocity(0) * velocity(1) + velocity(1) * velocity(1) / 2);
          result(1, 1) = m2 * l1 * l2 * cos(q(1)) * velocity(0) * velocity(0) / 2;
        });

    arm->setComputeJacobianFgyrOver_velocityFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = -m2 * l1 * l2 * sin(q(1)) * velocity(1);
          result(1, 0) = m2 * l1 * l2 * sin(q(1)) * velocity(0);
          result(0, 1) = -m2 * l1 * l2 * sin(q(1)) * (velocity(0) + velocity(1));
          result(1, 1) = 0;
        });

    arm->setComputeFintFunction(u_func);

    arm->setComputeJacobianFintOver_qFunction(jacobian_fint_over_q_func);

    arm->setComputeJacobianFintOver_velocityFunction(jacobian_fint_over_v_func);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.
    // -- relations --

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation01 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    auto inter01 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation01);
    relation01->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                       Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = 2 + l1 * cos(q(0)) + l2 * cos(q(0) + q(1));
    });

    auto relation02 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

    relation02->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                       Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = l1 * sin(q(0)) + l2 * sin(q(0) + q(1));
    });

    relation02->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = l1 * cos(q(0)) + l2 * cos(q(0) + q(1));
          result(0, 1) = l2 * cos(q(0) + q(1));
        });
    auto inter02 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation02);

    auto relation31 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation31->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                       Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = 0.7 - l1 * cos(q(0)) - l2 * cos(q(0) + q(1));
    });

    relation31->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = l1 * sin(q(0)) + l2 * sin(q(0) + q(1));
          result(0, 1) = l2 * sin(q(0) + q(1));
        });

    auto inter31 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation31);
    auto relation32 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation32->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                       Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = 2 + l1 * sin(q(0)) + l2 * sin(q(0) + q(1));
    });

    auto inter32 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation32);

    auto H10 = std::make_shared<Matrix>(1, 2);
    auto b10 = std::make_shared<Vector>(1);
    H10->setZero();
    (*H10)(0, 0) = -1;
    (*b10)(0) = std::numbers::pi;

    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e2);
    auto relation10 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H10, *b10);
    auto inter10 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation10);

    auto H11 = std::make_shared<Matrix>(1, 2);
    H11->setZero();
    (*H11)(0, 0) = 1;
    auto relation11 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H11);
    auto inter11 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation11);

    auto H20 = std::make_shared<Matrix>(1, 2);
    auto b20 = std::make_shared<Vector>(1);
    H20->setZero();
    (*H20)(0, 1) = -1;
    (*b20)(0) = 0.0001;

    auto relation20 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H20, *b20);
    auto inter20 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation20);

    auto H21 = std::make_shared<Matrix>(1, 2);
    auto b21 = std::make_shared<Vector>(1);
    H21->setZero();
    (*H21)(0, 1) = 1;
    (*b21)(0) = std::numbers::pi - 0.0001;

    auto relation21 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H21, *b21);
    auto inter21 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation21);

    // -------------
    // --- Model ---
    // -------------

    auto Manipulator = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    Manipulator->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->link(inter01, arm);
    Manipulator->link(inter02, arm);
    Manipulator->link(inter31, arm);
    Manipulator->link(inter32, arm);
    Manipulator->link(inter10, arm);
    Manipulator->link(inter20, arm);
    Manipulator->link(inter11, arm);
    Manipulator->link(inter21, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto s = std::make_shared<siconos::simulation::TimeStepping>(Manipulator, t);
    s->setNewtonTolerance(criterion);
    s->setNewtonMaxIteration(maxIter);
    // s->setResetAllLambda(false);
    //  -- OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(0.500001);
    s->insertIntegrator(OSI);
    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    s->insertNonSmoothProblem(osnspb);
    std::cout << "=== End of model loading === \n";

    // =========================== End of model definition ===========================

    // ================================= Computation
    unsigned int k = 0;
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15;
    Matrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    auto q = arm->q_read();
    auto v = arm->velocity_read();
    auto p = arm->p(1);

    // Initialization of the dicrete parameter z needs the followinf first computations
    arm->computeTotalForces(v, q, t0);
    arm->computeJacobianTotalForcesOver_velocity(v, q, t0);
    arm->computeJacobianTotalForcesOver_q(v, q, t0);
    inter01->computeOutput(t0, 0);
    inter02->computeOutput(t0, 0);
    inter31->computeOutput(t0, 0);
    inter32->computeOutput(t0, 0);
    inter10->computeOutput(t0, 0);
    inter20->computeOutput(t0, 0);
    inter11->computeOutput(t0, 0);
    inter21->computeOutput(t0, 0);

    // EventsManager * eventsManager = s->eventsManager();

    dataPlot(k, 0) = Manipulator->t0();
    dataPlot(k, 1) = q(0);
    dataPlot(k, 2) = q(1);
    dataPlot(k, 3) = (*inter02->y(0))(0);
    dataPlot(k, 4) = v(0);
    dataPlot(k, 5) = v(1);
    dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
    dataPlot(k, 7) = nimpact;  //(*inter->y(1))(1);
    dataPlot(k, 8) = param[6];
    dataPlot(k, 9) = param[4];
    dataPlot(k, 10) = param[24];
    dataPlot(k, 11) = test;
    dataPlot(k, 12) = (*p)(1);
    dataPlot(k, 13) = param[22];
    dataPlot(k, 14) = param[23];

    std::cout << "====> Start computation ... \n\n";

    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      param[0] = q(0);
      param[1] = q(1);
      param[2] = v(0);
      param[3] = v(1);
      param[16] = param[14];
      param[17] = param[15];
      param[20] = param[18];
      param[21] = param[19];

      // get current time step
      k++;

      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = q(1);
      dataPlot(k, 3) = (*inter02->y(0))(0);
      dataPlot(k, 4) = v(0);
      dataPlot(k, 5) = v(1);
      dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
      dataPlot(k, 7) = nimpact;  //(*inter->y(1))(1);
      dataPlot(k, 8) = param[6];
      if (test == 3)
        dataPlot(k, 9) = param[4] / h;
      else
        dataPlot(k, 9) = param[4];
      if (test == 5)
        dataPlot(k, 10) = param[24] / h;
      else
        dataPlot(k, 10) = param[24];
      dataPlot(k, 11) = test;
      dataPlot(k, 13) = param[22];
      dataPlot(k, 14) = param[23];

      s->advanceToEvent();
      // siconos::algebra::print(*osnspb);
      dataPlot(k, 12) = (*p)(1);
      param[4] = (inter02->getLambda(1))(0);
      param[24] = (inter31->getLambda(1))(0);
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k, 3) <= 0.01) && (test == 0) && (dataPlot(k, 6) < 0.3)) {
        param[8] = dataPlot(k, 0);
        param[5] = 0.7 + 0.5 * cos(std::numbers::pi / 2 +
                                   2 * std::numbers::pi * (param[8]) / param[11]);
        param[7] = param[9];
        arm->setComputeFintFunction(u10_func);
        test = 1;
      }

      //  controller during impacts accumulation phase after the first impact
      if ((dataPlot(k, 12) > 0) && (test == 1)) {
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u11_func);
        test = 2;
      }
      if ((dataPlot(k, 12) > 0) && (test == 2)) nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((dataPlot(k, 12) > 0) && (test == 2) && (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1) &&
          (fabs((*inter02->y(0))(0)) < 1e-6)) {
        // L= dataPlot(k,0)-(*z)(8);
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u20_func);
        test = 3;
        nimpact = 0;
      }
      //  controller during impacts accumulation phase after the first impact
      if ((param[24] > 0) && (test == 3)) {
        arm->setComputeFintFunction(u21_func);
        test = 4;
      }
      if ((param[24] > 0) && (test == 4)) nimpact = nimpact + 1;
      // controller during constraint-motion phase.
      if ((param[24] > 0) && (test == 4) &&
          (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1))  // && (fabs((*inter0->y(1))(0))<1e-6))
      {
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u22_func);
        test = 5;
        nimpact = 0;
      }
      // change of control law with a particular design of the desired trajectory that
      // guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / param[11]) > trunc((dataPlot(k, 0)) / param[11])) &&
          (test == 5)) {
        param[8] = dataPlot(k, 0) + h;
        param[10] = param[12];
        arm->setComputeFintFunction(u3_func);
        test = 6;
        // L = 0;
      }

      //  controller during free-motion phase
      if ((param[13] - 0.1 >= 0) && (test == 6)) {
        arm->setComputeFintFunction(u_func);
        test = 0;
        param[13] = 0;
      }
    }
    std::cout << "\nEnd of computation - Number of iterations done: " << k << "\n";
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    siconos::algebra::io::write("TwoLinkManipulator_MultiConstraints.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "TwoLinkManipulator_MultiConstraints.ref", eps)) > eps)
      return 1;
    else
      return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
