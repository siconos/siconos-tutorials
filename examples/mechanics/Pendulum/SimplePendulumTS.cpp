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

// =============================== Simple Pendulum Example ===============================
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

constexpr double gravity = 10.;
double m1 = 1.0;
double l1 = 1.0;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 1;  // degrees of freedom for robot arm
    double t0 = 0;          // initial computation time
    double T = 50.0;        // final computation time
    double h = 0.0005;      // time step
    double criterion = 0.00005;
    unsigned int maxIter = 2000;
    double e = 0.9;  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: Simple Pendulum ---

    // Initial position (angles in radian)
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();
    q0(0) = 1;

    auto simplependulum = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    Matrix mass{nDof, nDof};
    mass(0, 0) = m1 * l1;
    simplependulum->setConstantMass(mass, siconos::algebra::alias_t);

    simplependulum->setComputeFintFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapVectorType> fint) {
          fint(0) = m1 * sin(q(0)) * gravity;
        });

    simplependulum->setComputeJacobianFintOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &v,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = cos(q(0)) * gravity * (m1);
        });

    // -------------------
    // --- Interactions---
    // -------------------
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation->setComputehFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::SiconosVector> y) { y(0) = l1 * sin(q(0)); });

    relation->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> result) { result(0, 0) = l1 * cos(q(0)); });

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------

    auto Pendulum = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Pendulum->insertDynamicalSystem(simplependulum);
    Pendulum->link(inter, simplependulum);

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

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    s->insertNonSmoothProblem(osnspb);
    std::cout << "=== End of model loading === \n";

    // --- Simulation initialization ---
    int k = 0;
    int N = ceil((T - t0) / h);
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    Matrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    auto q = simplependulum->q_read();
    auto vel = simplependulum->velocity_read();
    dataPlot(k, 0) = Pendulum->t0();
    dataPlot(k, 1) = q(0);
    dataPlot(k, 2) = vel(0);
    dataPlot(k, 3) = l1 * sin(q(0));
    dataPlot(k, 4) = -l1 * cos(q(0));
    dataPlot(k, 5) = l1 * cos(q(0)) * vel(0);
    // --- Compute elapsed time ---
    auto start = std::chrono::system_clock::now();
    //    EventsManager eventsManager = s->eventsManager();
    // --- Time loop ---
    std::cout << "Start computation ... \n";
    std::cout << "Number of time steps " << N << "\n";
    while (s->hasNextEvent()) {
      k++;
      if (!(div(k, 10000).rem)) std::cout << "Step number " << k << "\n";

      // Solve problem
      s->advanceToEvent();
      // Data Output
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = vel(0);
      dataPlot(k, 3) = l1 * sin(q(0));
      dataPlot(k, 4) = -l1 * cos(q(0));
      dataPlot(k, 5) = l1 * cos(q(0)) * vel(0);

      s->nextStep();
      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    siconos::algebra::io::write("SimplePendulumResult.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SimplePendulumResult.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
