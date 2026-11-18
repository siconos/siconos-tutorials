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

// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

double gravity = 10.0;
double m1 = 1.0;
double m2 = 1.0;
double l1 = 1.0;
double l2 = 1.0;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 2;  // degrees of freedom for robot arm
    double t0 = 0;          // initial computation time
    double T = 5.0;         // final computation time
    double h = 0.0005;      // time step
    double criterion = 0.05;
    unsigned int maxIter = 20000;
    double e = 1.0;  // nslaw
    double e1 = 0.0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: Double Pendulum ---

    // Initial position (angles in radian)
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    // q0(0) = 1.5;
    // q0(1) = 1.5;

    // for sympy plugins uncomment below (relative parametrization)
    // Note, we have the relation :
    // absolute[q0(0)(0)] + relative[q0(0)(1)] = absolute[q0(0)(1)]
    // q0(0) = 0.1;
    // q0(1) = 0.1;

    q0(0) = 0.1;
    q0(1) = 0.2;

    /*REGULAR PLUGINS - uncomment to use*/
    auto doublependulum = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    doublependulum->setComputeMassFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> mass) {
          mass(0, 0) = (m1 + m2) * l1;
          mass(1, 0) = m2 * l1 * cos(q(0) - q(1));
          mass(0, 1) = m2 * l2 * cos(q(0) - q(1));
          mass(1, 1) = m2 * l2;
        });

    doublependulum->setComputeFgyrFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapVectorType> fgyr) {
          fgyr(0) = m2 * l2 * velocity(1) * velocity(1) * sin(q(0) - q(1));
          fgyr(1) = -m2 * l1 * velocity(0) * velocity(0) * sin(q(0) - q(1));
        });

    // set 'random' value for jacobians, whatever fgyr is, just for tests
    doublependulum->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = m2 * l2 * velocity(1) * velocity(0) * cos(q(0) - q(1));
          jacob(1, 0) = -m2 * l1 * velocity(0) * velocity(0) * cos(q(0) - q(1));
          jacob(0, 1) = -m2 * l2 * velocity(1) * velocity(1) * cos(q(0) - q(1));
          jacob(1, 1) = m2 * l1 * velocity(0) * velocity(0) * cos(q(0) - q(1));
        });

    doublependulum->setComputeJacobianFgyrOver_velocityFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = 0.0;
          jacob(1, 0) = -2.0 * m2 * l1 * velocity(0) * sin(q(0) - q(1));
          jacob(0, 1) = 2.0 * m2 * l2 * velocity(1) * sin(q(0) - q(1));
          jacob(1, 1) = 0.0;
        });

    doublependulum->setComputeFintFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapVectorType> fint) {
          fint(0) = sin(q(0)) * gravity * (m1 + m2);
          fint(1) = sin(q(1)) * gravity * m2;
        });

    doublependulum->setComputeJacobianFintOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = cos(q(0)) * gravity * (m1 + m2);
          jacob(1, 0) = 0.0;
          jacob(0, 1) = 0.0;
          jacob(1, 1) = cos(q(1)) * gravity * (m2);
          ;
        });

    /*SYMPY PLUGINS - uncomment to use*/
    // auto doublependulum= std::make_shared<siconos::modeling::LagrangianDS>(q0, v0,
    // "DoublePendulumSymPyPlugin:mass");
    // doublependulum->setComputeFGyrFunction("DoublePendulumSymPyPlugin", "FGyr");
    // doublependulum->setComputeJacobianFGyrqDotFunction("DoublePendulumSymPyPlugin",
    // "jacobianVFGyr");
    // doublependulum->setComputeJacobianFGyrqFunction("DoublePendulumSymPyPlugin",
    // "jacobianFGyrq");

    // -------------------
    // --- Interactions---
    // -------------------

    // -- relations --

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

    relation->setComputehFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::SiconosVector> y) { y(0) = l1 * sin(q(0)); });

    relation->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = l1 * cos(q(0));
          result(0, 1) = 0.0;
        });

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e1);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();

    relation1->setComputehFunction([](const siconos::algebra::BlockVector &q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = l1 * sin(q(0)) + l2 * sin(q(1));
    });

    relation1->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = l1 * cos(q(0));
          result(0, 1) = l2 * cos(q(1));
        });

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -------------
    // --- Model ---
    // -------------

    auto Pendulum = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Pendulum->insertDynamicalSystem(doublependulum);
    Pendulum->link(inter, doublependulum);
    Pendulum->link(inter1, doublependulum);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto s = std::make_shared<siconos::simulation::TimeStepping>(Pendulum, t);
    s->setNewtonTolerance(criterion);
    s->setNewtonMaxIteration(maxIter);
    // -- OneStepIntegrators --

    // double theta=0.500001;
    double theta = 0.500001;

    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    s->insertIntegrator(OSI);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    s->insertNonSmoothProblem(osnspb);
    std::cout << "=== End of model loading === \n";

    // =========================== End of model definition ===========================
    // dataPlot(k,7) = (*inter->y(0))(0);

    // ================================= Computation =================================

    int k = 0;
    int N = ceil((T - t0) / h);
    std::cout << "Number of time step   " << N << "\n";
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    Matrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    auto q = doublependulum->q();
    auto v = doublependulum->velocity();

    dataPlot(k, 0) = t0;
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*v)(1);
    dataPlot(k, 5) = l1 * sin((*q)(0));
    dataPlot(k, 6) = -l1 * cos((*q)(0));
    dataPlot(k, 7) = l1 * sin((*q)(0)) + l2 * sin((*q)(1));
    dataPlot(k, 8) = -l1 * cos((*q)(0)) - l2 * cos((*q)(1));
    dataPlot(k, 9) = l1 * cos((*q)(0)) * ((*v)(0));
    dataPlot(k, 10) = l1 * cos((*q)(0)) * ((*v)(0)) + l2 * cos((*q)(1)) * ((*v)(1));

    auto start = std::chrono::system_clock::now();

    // --- Time loop ---
    std::cout << "Start computation ... \n";

    while (s->hasNextEvent()) {
      k++;

      //  if (!(div(k,1000).rem))  cout <<"Step number "<< k << "\n";

      // Solve problem
      s->advanceToEvent();
      // Data Output
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*q)(1);
      // sympy plugin with relative parametrization:
      // dataPlot(k, 3) = (*q)(0) + (*q)(1);
      dataPlot(k, 4) = (*v)(1);
      // sympy plugin with relative parametrization:
      // dataPlot(k, 4) = (*v)(0) + (*v)(1);
      dataPlot(k, 5) = l1 * sin((*q)(0));
      dataPlot(k, 6) = -l1 * cos((*q)(0));
      dataPlot(k, 7) = l1 * sin((*q)(0)) + l2 * sin((*q)(1));
      dataPlot(k, 8) = -l1 * cos((*q)(0)) - l2 * cos((*q)(1));
      dataPlot(k, 9) = l1 * cos((*q)(0)) * ((*v)(0));
      dataPlot(k, 10) = l1 * cos((*q)(0)) * ((*v)(0)) + l2 * cos((*q)(1)) * ((*v)(1));
      s->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("DoublePendulumResult.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "DoublePendulumResult.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
