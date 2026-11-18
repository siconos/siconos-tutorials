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

#include "TwoLinksFlexibleMultiPlugins.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace two_links_flexible_multi_plugins;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    double t0 = 0;    // initial computation time
    double T = 30;    // final computation time
    double h = 1e-3;  // time step
    double criterion = 1e-8;
    unsigned int maxIter = 20000;
    double e = 0.0;
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
    q0(2) = 1.5;
    q0(3) = -0.9;

    param[0] = q0(0);
    param[1] = q0(1);
    param[2] = v0(0);
    param[3] = v0(1);
    param[11] = std::numbers::pi;

    auto arm = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    arm->setComputeMassFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                                   Eigen::Ref<siconos::algebra::MapType> mass) {
      mass(0, 0) =
          m1 * (l1 * l1 / 4) + I1 + I2 + m2 * (l1 * l1 + (l2 * l2 / 4) + l1 * l2 * cos(q(1)));
      mass(1, 0) = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q(1)) / 2;
      mass(0, 1) = I2 + m2 * l2 * l2 / 4 + m2 * l1 * l2 * cos(q(1)) / 2;
      mass(1, 1) = I2 + m2 * l2 * l2 / 4;
      mass(2, 2) = J1;
      mass(3, 3) = J2;
    });

    arm->setComputeFgyrFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapVectorType> result) {
          result(0) = -m2 * l1 * l2 * sin(q(1)) *
                          (velocity(0) * velocity(1) + velocity(1) * velocity(1) / 2) +
                      g * (l1 * cos(q(0)) * (m2 + m1 / 2) + m2 * l2 * cos(q(0) + q(1)) / 2) +
                      K1 * (q(0) - q(2));
          result(1) = m2 * l1 * l2 * sin(q(1)) * velocity(0) * velocity(0) / 2 +
                      g * m2 * l2 * cos(q(0) + q(1)) / 2 + K2 * (q(1) - q(3));
          result(2) = K1 * (q(2) - q(0));
          result(3) = K2 * (q(3) - q(1));
        });

    arm->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) =
              -g * (l1 * sin(q(0)) * (m2 + m1 / 2) + m2 * l2 * sin(q(0) + q(1)) / 2) + K1;
          result(1, 0) = -g * m2 * l2 * sin(q(0) + q(1)) / 2;

          result(2, 0) = -K1;
          result(0, 1) = -m2 * l1 * l2 * cos(q(1)) *
                             (velocity(0) * velocity(1) + velocity(1) * velocity(1) / 2) -
                         g * m2 * l2 * sin(q(0) + q(1)) / 2;

          result(1, 1) = m2 * l1 * l2 * cos(q(1)) * velocity(0) * velocity(0) / 2 -
                         g * m2 * l2 * sin(q(0) + q(1)) / 2 + K2;
          result(3, 1) = -K2;
          result(0, 2) = -K1;
          result(2, 2) = K1;
          result(1, 3) = -K2;
          result(3, 3) = K2;
        });

    arm->setComputeJacobianFgyrOver_velocityFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) = -m2 * l1 * l2 * sin(q(1)) * velocity(1);
          result(1, 0) = m2 * l1 * l2 * sin(q(1)) * velocity(0);
          result(0, 1) = -m2 * l1 * l2 * sin(q(1)) * (velocity(0) + velocity(1));
        });

    arm->setComputeFintFunction(u_func);
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

    relation02->setComputeHFunction(
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

    relation31->setComputeHFunction(
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

    // -------------
    // --- Model ---
    // -------------

    auto Manipulator = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    Manipulator->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->link(inter01, arm);
    Manipulator->link(inter02, arm);
    // link the interaction and the dynamical system
    Manipulator->link(inter31, arm);
    Manipulator->link(inter32, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto s = std::make_shared<siconos::simulation::TimeStepping>(Manipulator, t);
    s->setNewtonTolerance(criterion);
    s->setNewtonMaxIteration(maxIter);
    // -- OneStepIntegrators --
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
    unsigned int outputSize = 16;
    Matrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    auto q = arm->q_read();
    auto v = arm->velocity_read();
    auto p = arm->p(1);

    // Initialization of the dicrete parameter z needs the followinf first computations
    arm->computeTotalForces(v, q, t0);
    arm->computeJacobianTotalForcesOver_q(v, q, t0);
    arm->computeJacobianTotalForcesOver_velocity(v, q, t0);
    inter01->computeOutput(t0, 0);
    inter02->computeOutput(t0, 0);
    inter31->computeOutput(t0, 0);
    inter32->computeOutput(t0, 0);

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
    dataPlot(k, 10) = param[18];
    dataPlot(k, 11) = test;
    dataPlot(k, 12) = (*p)(1);
    dataPlot(k, 13) = param[14];
    dataPlot(k, 14) = param[15];
    dataPlot(k, 15) = param[16];

    bool stop = 0;

    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      param[0] = q(0);
      param[1] = q(1);
      param[2] = v(0);
      param[3] = v(1);
      // param[14] = (*q)(2);
      //         param[15] = (*q)(3);
      //  param[16] = (*v)(2);
      //         param[17] = (*v)(3);

      // get current time step
      k++;

      //      relation->computeOutput(s->nextTime());
      //  if(k==1106) stop = 1;

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
        dataPlot(k, 10) = param[18] / h;
      else
        dataPlot(k, 10) = param[18];
      dataPlot(k, 11) = test;
      dataPlot(k, 13) = param[14];
      dataPlot(k, 14) = param[15];
      dataPlot(k, 15) = param[16];  //(19)*param[19];

      s->advanceToEvent();
      dataPlot(k, 12) = (*p)(1);
      param[4] = (inter02->getLambda(1))(0);
      param[18] = (inter31->getLambda(1))(0);
      //  if(k==41000)
      //    {
      //      (*v)(0) = 0.0;(*v)(1) = 0.0;(*v)(2) = 0.0;(*v)(3) = 0.0;
      //    }
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k - 1, 14) <= 0.1) && (test == 0) && (dataPlot(k, 13) < 0.6)) {
        param[8] = dataPlot(k, 0);
        param[5] = param[14];
        param[10] = dataPlot(k, 3);
        param[7] = param[9];
        arm->setComputeFintFunction(u1_func);
        test = 1;
      }

      // controller during impacts accumulation phase after the first impact
      if ((param[4] > 0) && (test == 1)) {
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u2_func);
        test = 2;
      }
      if ((param[4] > 0) && (test == 2)) nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((param[4] > 0) && (test == 2) &&
          (dataPlot(k, 7) - dataPlot(k - 3, 7) == 3))  //  && (fabs((*inter0->y(1))(1))<1e-6))
      {
        // L= dataPlot(k,0)-param[8];
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u3_func);
        test = 3;
        nimpact = 0;
      }
      //  controller during impacts accumulation phase after the first impact
      if ((dataPlot(k, 10) > 0) && (test == 3)) {
        arm->setComputeFintFunction(u4_func);
        test = 4;
      }

      if ((param[18] > 0) && (test == 4)) nimpact = nimpact + 1;
      // controller during constraint-motion phase.
      if ((param[18] > 0) && (test == 4) &&
          (dataPlot(k, 7) - dataPlot(k - 3, 7) == 3))  // && (fabs((*inter0->y(1))(0))<1e-6))
      {
        param[8] = dataPlot(k, 0);
        arm->setComputeFintFunction(u5_func);
        test = 5;
        nimpact = 0;
      }
      // change of control law with a particular design of the desired trajectory that
      // guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / param[11]) > trunc((dataPlot(k, 0)) / param[11])) &&
          (test == 5)) {
        param[8] = dataPlot(k, 0) + h;
        param[10] = param[12];
        arm->setComputeFintFunction(u6_func);
        test = 6;
        // L = 0;
      }

      //  controller during free-motion phase
      if ((param[13] >= 0) && (test == 6)) {
        arm->setComputeFintFunction(u_func);
        test = 0;
        param[13] = 0;
      }

      if (stop) break;
    }
    std::cout << "\nEnd of computation - Number of iterations done: " << k << "\n";
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    siconos::algebra::io::write("TwoLinkManipulator_MultiConstraintsFlexible.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "TwoLinkManipulator_MultiConstraintsFlexible.ref", eps)) > eps)
      return 1;
    else
      return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
