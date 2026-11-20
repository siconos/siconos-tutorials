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

/*
 *C++ input file, MoreauJeanOSI-Time-Stepping version
 */

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

#define WITH_FRICTION
// #define DISPLAY_INTER

// Inertial parameters
double m1 = 1.0;
double m2 = 1.0;
double l = 1.0;
double a = 0.025;
double d = 0.25;

// force elements
double gravity = 9.81;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    int nDof = 3;     // degrees of freedom for robot arm
    double t0 = 0;    // initial computation time
    double T = 3.0;   // final computation time
    double h = 1e-3;  // time step : do not decrease, because of strong penetrations

    // contact parameters
    double eN = 0.5;
    double eT = 0.;
    double mu = 0.8;

    // initial conditions
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    q0(0) = 0.1;
    q0(2) = 0.1;
    v0(0) = 2.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    std::cout << "====> Model loading ...\n\n";

    auto pendulum =
        std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    pendulum->setComputeMassFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> mass) {
          mass.setZero();
          mass(0, 0) = m1 + m2;
          mass(2, 0) = cos(q(2));
          mass(1, 1) = 1.;
          mass(0, 2) = m2 * l * cos(q(2));
          mass(2, 2) = l;
        });

    pendulum->setComputeFgyrFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapVectorType> fgyr) {
          fgyr.setZero();
          fgyr(0) = -m2 * l * velocity(2) * velocity(2) * sin(q(2));
        });

    // set 'random' value for jacobians, whatever fgyr is, just for tests
    pendulum->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(2, 0) = -m2 * l * velocity(2) * velocity(2) * cos(q(2));
        });

    pendulum->setComputeJacobianFgyrOver_velocityFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(2, 0) = -2.0 * m2 * velocity(2) * sin(q(2));
        });

    pendulum->setComputeFintFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
           Eigen::Ref<siconos::algebra::MapVectorType> fint) {
          fint(0) = 0;
          fint(1) = gravity * q(1);
          fint(2) = gravity * sin(q(2));
        });

    pendulum->setComputeJacobianFintOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob.setZero();
          jacob(1, 0) = gravity;
          jacob(2, 2) = gravity * cos(q(2));
        });

    // -------------------
    // --- Interactions---
    // -------------------
    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN, eT, mu, 2);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation1->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(1);
      y(1) = 0.;
    });

    Matrix Jhq1{2, nDof};
    Jhq1.setZero();
    Jhq1(1, 0) = 1.;
    Jhq1(0, 1) = 1.;
    //    relation1->setConstantContactMatrix(Jhq1);
    relation1->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(1, 0) = 1.;
          result(0, 1) = 1.;
        });

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    auto nslaw2 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN, 0.0, 0.0, 2);
    // auto nslaw2= std::make_shared<siconos::modeling::NewtonImpactNSL>(eN);
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation2->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = d - q(0) - 0.5 * a;  // normal
      y(1) = 0.;
    });
    Matrix Jhq2{2, nDof};
    Jhq2.setZero();
    Jhq2(0, 0) = -1.;
    // relation2->setConstantContactMatrix(Jhq2);
    relation2->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) = -1.;
        });
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    auto nslaw3 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN, 0.0, 0.0, 2);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation3->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(0) - 0.5 * a;  // normal
      y(1) = 0.;
    });
    Matrix Jhq3{2, nDof};
    Jhq3.setZero();
    Jhq3(0, 0) = 1.;
    //    relation3->setConstantContactMatrix(Jhq3);
    relation3->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) = 1.;
        });
    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // -------------
    // --- Model ---
    // -------------
    auto pendulumWithSlider =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    pendulumWithSlider->insertDynamicalSystem(pendulum);
    pendulumWithSlider->link(inter1, pendulum);
    pendulumWithSlider->link(inter2, pendulum);
    pendulumWithSlider->link(inter3, pendulum);

    // ----------------
    // --- Simulation ---
    // ----------------
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(0.5, 1.);

    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    auto impact = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(
        2, SICONOS_FRICTION_2D_ENUM);
    // auto impact= std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_ENUM);

    impact->numericsSolverOptions()->dparam[0] = 1e-08;
    impact->numericsSolverOptions()->iparam[0] = 100;
    impact->numericsSolverOptions()->iparam[2] = 1;  // random
    auto s = std::make_shared<siconos::simulation::TimeStepping>(pendulumWithSlider, t);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    auto topo = pendulumWithSlider->topology();

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 14;
    Matrix dataPlot(N + 1, outputSize);

    auto q = pendulum->q();
    auto v = pendulum->velocity();

    dataPlot(0, 0) = pendulumWithSlider->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) = (*inter1->y(0))(0);  // g1
    dataPlot(0, 8) = (*inter2->y(0))(0);  // g2
    dataPlot(0, 9) = (*inter3->y(0))(0);  // g3

    // --- Time loop ---
    std::cout << "====> Start computation ... \n\n";

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();

    while ((s->hasNextEvent()) && (k <= 3000))
    //    while ((s->hasNextEvent()))
    {
      // std::cout <<"t = " <<s->nextTime()-h  <<std::endl;
      // //std::cout <<"=====================================================" <<std::endl;
      // cout << "q[0] = " << (*q)(0)  << endl;
      // cout << "q[1] = " << (*q)(1)  << endl;
      // cout << "q[2] = " << (*q)(2)  << endl;
      // cout << "v[0] = " << (*v)(0)  << endl;
      // cout << "v[1] = " << (*v)(1)  << endl;
      // cout << "v[2] = " << (*v)(2)  << endl;

      // std::cout << "=============== Step k ="<< k<< std::endl;
      s->advanceToEvent();
      impact->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = (*inter1->y(0))(0);  // g1
      dataPlot(k, 8) = (*inter2->y(0))(0);  // g2
      dataPlot(k, 9) = (*inter3->y(0))(0);  // g3
      dataPlot(k, 10) = s->getNewtonNbIterations();
      auto indexSet1 = topo->indexSet(1);
      dataPlot(k, 11) = indexSet1->size();

      // if (indexSet1->size() > 5)
      // {
      //   siconos::algebra::print(*impact);
      // }
      //      if (s->nextTime() > 0.035 and (*inter1->lambda(1))(0) >0.0)
#ifdef DISPLAY_INTER
      std::cout << "=============== Step k =" << k << std::endl;
      std::cout << "Time " << s->nextTime() << std::endl;

      siconos::algebra::print(*impact);
      std::cout << " (*inter1->lambda(1))(0) " << (*inter1->lambda(1))(0) << std::endl;
      std::cout << " (*inter2->lambda(1))(0) " << (*inter2->lambda(1))(0) << std::endl;
#endif

      s->processEvents();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("PendulumSlider.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    // all the columns but 11 (size of indexSet1): a contact of the pendulum stays at a
    // gap of +/-1e-17 (round-off), so the number of active contacts depends on the platform
    std::vector<int> compared_columns;
    for (int j = 0; j < static_cast<int>(outputSize); ++j)
      if (j != 11) compared_columns.push_back(j);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "PendulumSlider.ref", eps,
                                                      compared_columns)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
