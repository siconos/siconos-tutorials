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

/*
  C++ input file, D1MinusLinearOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a D1MinusLinearOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

#include "SCConst.h"  // Simulation parameters

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace parameters;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    int nDof = 3;  // degrees of freedom for the slider crank
    double t0 = 0.;         // initial computation time
    double T = 0.2;         // final computation time
    double h = 1e-5;        // time step : do not decrease, because of strong penetrations

    // contact parameters
    double eN1 = 0.4;
    double eN2 = 0.4;
    double eN3 = 0.4;
    double eN4 = 0.4;
    // eN1 = 0.1;
    // eN2 = 0.1;
    // eN3 = 0.1;
    // eN4 = 0.1;
    // double mu1 = 0.01;
    // double mu2 = 0.01;
    // double mu3 = 0.01;
    // double mu4 = 0.01;

    // initial conditions
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    v0(0) = 150.;
    v0(1) = -75.;

    // t0 = 7e-5;
    // q0(0)=  1.129178e-02;
    // q0(1)= -5.777764e-03;
    // q0(2)=  0.000000e+00;

    // v0(0) = 1.971606e+02 ;
    // v0(1) = -1.064301e+02;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    std::cout << "====> Model loading ...\n\n";

    auto slider = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
    slider->setComputeMassFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> mass) {
          mass.setZero();
          mass(0, 0) = J1 + (0.25 * m1 + m2 + m3) * l1 * l1;
          mass(1, 0) = (0.5 * m2 + m3) * l1 * l2 * cos(q(1) - q(0));

          mass(0, 1) = (0.5 * m2 + m3) * l1 * l2 * cos(q(1) - q(0));
          mass(1, 1) = J2 + (0.25 * m2 + m3) * l2 * l2;
          mass(2, 2) = J3;
        });

    slider->setComputeFgyrFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapVectorType> fgyr) {
          fgyr(0) = (0.5 * m2 + m3) * l1 * l2 * sin(q(0) - q(1)) * velocity(1) * velocity(1);
          fgyr(1) = -(0.5 * m2 + m3) * l1 * l2 * sin(q(0) - q(1)) * velocity(0) * velocity(0);
          fgyr(2) = 0.;
        });

    // set 'random' value for jacobians, whatever fgyr is, just for tests
    slider->setComputeJacobianFgyrOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob.setZero();
          jacob(0, 0) =
              (0.5 * m2 + m3) * l1 * l2 * cos(q(0) - q(1)) * velocity(1) * velocity(1);
          jacob(1, 0) =
              -(0.5 * m2 + m3) * l1 * l2 * cos(q(0) - q(1)) * velocity(0) * velocity(0);

          jacob(0, 1) =
              -(0.5 * m2 + m3) * l1 * l2 * cos(q(0) - q(1)) * velocity(1) * velocity(1);
          jacob(1, 1) =
              (0.5 * m2 + m3) * l1 * l2 * cos(q(0) - q(1)) * velocity(0) * velocity(0);
        });

    slider->setComputeJacobianFgyrOver_velocityFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob.setZero();
          jacob(1, 0) = -2. * (0.5 * m2 + m3) * l1 * l2 * sin(q(0) - q(1)) * velocity(0);
          jacob(0, 1) = 2. * (0.5 * m2 + m3) * l1 * l2 * sin(q(0) - q(1)) * velocity(1);
        });

    slider->setComputeFintFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapVectorType> fint) {
          fint(0) = (0.5 * m1 + m2 + m3) * gravity * l1 * cos(q(0));
          fint(1) = (0.5 * m2 + m3) * gravity * l2 * cos(q(1));
        });

    slider->setComputeJacobianFintOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
           const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob.setZero();
          jacob(0, 0) = -(0.5 * m1 + m2 + m3) * gravity * l1 * sin(q(0));
          jacob(1, 1) = -(0.5 * m2 + m3) * gravity * l2 * sin(q(1));
        });

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN1);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation1->setComputehFunction([](const siconos::algebra::BlockVector &q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = 0.5 * d -
             (l1 * sin(q(0)) + l2 * sin(q(1)) - a * sin(q(2)) + b * cos(q(2)));  // normal
    });

    relation1->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = -l1 * cos(q(0));
          jacob(0, 1) = -l2 * cos(q(1));
          jacob(0, 2) = a * cos(q(2)) + b * sin(q(2));
        });

    relation1->setComputeHdotFunction(
        [](const siconos::algebra::BlockVector &q, const siconos::algebra::BlockVector &qdot,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = parameters::l1 * sin(q(0)) * qdot(0);
          jacob(0, 1) = parameters::l2 * sin(q(1)) * qdot(1);
          jacob(0, 2) = (-parameters::a * sin(q(2)) + parameters::b * cos(q(2))) * qdot(2);
        });

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -- corner 2 --
    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN2);
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation2->setComputehFunction([](const siconos::algebra::BlockVector &q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = 0.5 * d -
             (l1 * sin(q(0)) + l2 * sin(q(1)) + a * sin(q(2)) + b * cos(q(2)));  // normal
    });

    relation2->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = -l1 * cos(q(0));
          jacob(0, 1) = -l2 * cos(q(1));
          jacob(0, 2) = -a * cos(q(2)) + b * sin(q(2));
        });

    relation2->setComputeHdotFunction(
        [](const siconos::algebra::BlockVector &q, const siconos::algebra::BlockVector &qdot,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = parameters::l1 * sin(q(0)) * qdot(0);
          jacob(0, 1) = parameters::l2 * sin(q(1)) * qdot(1);
          jacob(0, 2) = (-parameters::a * sin(q(2)) + parameters::b * cos(q(2))) * qdot(2);
        });

    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // -- corner 3 --
    auto nslaw3 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN3);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation3->setComputehFunction([](const siconos::algebra::BlockVector &q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) =
          0.5 * d + l1 * sin(q(0)) + l2 * sin(q(1)) - a * sin(q(2)) - b * cos(q(2));  // normal
    });

    relation3->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = l1 * cos(q(0));
          jacob(0, 1) = l2 * cos(q(1));
          jacob(0, 2) = -a * cos(q(2)) + b * sin(q(2));
        });

    relation3->setComputeHdotFunction(
        [](const siconos::algebra::BlockVector &q, const siconos::algebra::BlockVector &qdot,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = -parameters::l1 * sin(q(0)) * qdot(0);
          jacob(0, 1) = -parameters::l2 * sin(q(1)) * qdot(1);
          jacob(0, 2) = (parameters::a * sin(q(2)) + parameters::b * cos(q(2))) * qdot(2);
        });

    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // -- corner 4 --
    auto nslaw4 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN4);
    auto relation4 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation4->setComputehFunction([](const siconos::algebra::BlockVector &q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) =
          0.5 * d + l1 * sin(q(0)) + l2 * sin(q(1)) + a * sin(q(2)) - b * cos(q(2));  // normal
    });

    relation4->setComputeHFunction(
        [](const siconos::algebra::BlockVector &q,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = l1 * cos(q(0));
          jacob(0, 1) = l2 * cos(q(1));
          jacob(0, 2) = a * cos(q(2)) + b * sin(q(2));
        });
    relation4->setComputeHdotFunction(
        [](const siconos::algebra::BlockVector &q, const siconos::algebra::BlockVector &qdot,
           Eigen::Ref<siconos::algebra::MapType> jacob) {
          jacob(0, 0) = -parameters::l1 * sin(q(0)) * qdot(0);
          jacob(0, 1) = -parameters::l2 * sin(q(1)) * qdot(1);
          jacob(0, 2) = (parameters::a * sin(q(2)) + parameters::b * cos(q(2))) * qdot(2);
        });

    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);

    // -------------
    // --- Model ---
    // -------------
    auto sliderWithClearance =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    sliderWithClearance->insertDynamicalSystem(slider);
    sliderWithClearance->link(inter1, slider);
    sliderWithClearance->link(inter2, slider);
    sliderWithClearance->link(inter3, slider);
    sliderWithClearance->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    auto ositype = siconos::integrators::D1MinusLinearOSI::Type::halfexplicit_velocity_level;
    auto OSI = std::make_shared<siconos::integrators::D1MinusLinearOSI>(ositype);
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto force = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    auto s =
        std::make_shared<siconos::simulation::TimeSteppingD1Minus>(sliderWithClearance, t, 2);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 35;
    Matrix dataPlot(N, outputSize);

    auto q = slider->q();
    auto v = slider->velocity();

    // computation for a first consistent output
    inter1->computeOutput(t0, 0);
    inter2->computeOutput(t0, 0);
    inter3->computeOutput(t0, 0);
    inter4->computeOutput(t0, 0);

    dataPlot(0, 0) = sliderWithClearance->t0();
    dataPlot(0, 1) = (*q)(0) / (2. * std::numbers::pi);  // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
        c;  // y corner 1 (normalized)
    dataPlot(0, 8) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
        c;  // y corner 2 (normalized)
    dataPlot(0, 9) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
        (c);  // y corner 3 (normalized)
    dataPlot(0, 10) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
        (c);  // y corner 4 (normalized)
    dataPlot(0, 11) =
        (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1;          // x slider (normalized)
    dataPlot(0, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c;  // y slider (normalized
    dataPlot(0, 13) = (*inter1->y(0))(0);                           // g1
    dataPlot(0, 14) = (*inter2->y(0))(0);                           // g2
    dataPlot(0, 15) = (*inter3->y(0))(0);                           // g3
    dataPlot(0, 16) = (*inter4->y(0))(0);                           // g4
    dataPlot(0, 17) = (*inter1->y(1))(0);                           // dot g1
    dataPlot(0, 18) = (*inter2->y(1))(0);                           // dot g2
    dataPlot(0, 19) = (*inter3->y(1))(0);                           // dot g3
    dataPlot(0, 20) = (*inter4->y(1))(0);                           // dot g4
    dataPlot(0, 21) = (*inter1->lambda(1))(0);                      // lambda1
    dataPlot(0, 22) = (*inter2->lambda(1))(0);                      // lambda1
    dataPlot(0, 23) = (*inter3->lambda(1))(0);                      // lambda3
    dataPlot(0, 24) = (*inter4->lambda(1))(0);                      // lambda4
    dataPlot(0, 25) = 0;
    dataPlot(0, 26) = 0;
    // dataPlot(k, 27) = (*inter1->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 28) = (*inter2->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 29) = (*inter3->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 30) = (*inter4->lambda(2))(0) ; // lambda1_{k+1}^-

    // not yet allocated
    // dataPlot(k, 31) = ( inter1->lambdaMemory(2).getSiconosVector(0) )(0) ; // lambda1old
    // dataPlot(k, 32) = ( inter2->lambdaMemory(2).getSiconosVector(0) )(0) ; // lambda1old
    // dataPlot(k, 33) = ( inter3->lambdaMemory(2).getSiconosVector(0) )(0) ; // lambda1old
    // dataPlot(k, 34) = ( inter4->lambdaMemory(2).getSiconosVector(0) )(0) ; // lambda1old

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    while ((s->hasNextEvent())) {
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"Iteration k = " << k <<std::endl;
      // std::cout <<"s->nextTime() = " <<s->nextTime()  <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;

      s->advanceToEvent();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2. * std::numbers::pi);  // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
          c;  // y corner 1 (normalized)
      dataPlot(k, 8) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
          c;  // y corner 2 (normalized)
      dataPlot(k, 9) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
          (c);  // y corner 3 (normalized)
      dataPlot(k, 10) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
          (c);  // y corner 4 (normalized)
      dataPlot(k, 11) =
          (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1;          // x slider (normalized)
      dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c;  // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0);                           // g1
      dataPlot(k, 14) = (*inter2->y(0))(0);                           // g2
      dataPlot(k, 15) = (*inter3->y(0))(0);                           // g3
      dataPlot(k, 16) = (*inter4->y(0))(0);                           // g4
      dataPlot(k, 17) = (*inter1->y(1))(0);                           // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0);                           // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0);                           // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0);                           // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0);                      // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0);                      // lambda1
      dataPlot(k, 23) = (*inter3->lambda(1))(0);                      // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0);                      // lambda4
      dataPlot(k, 25) = 0;
      dataPlot(k, 26) = 0;
      dataPlot(k, 27) = (*inter1->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 28) = (*inter2->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 29) = (*inter3->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 30) = (*inter4->lambda(2))(0);  // lambda1_{k+1}^-

      dataPlot(k, 31) = (inter1->lambdaMemory(2).getSiconosVector(0))(0);  // lambda1_k^+
      dataPlot(k, 32) = (inter2->lambdaMemory(2).getSiconosVector(0))(0);  // lambda2_k^+
      dataPlot(k, 33) = (inter3->lambdaMemory(2).getSiconosVector(0))(0);  // lambda3_k^+
      dataPlot(k, 34) = (inter4->lambdaMemory(2).getSiconosVector(0))(0);  // lambda4_k^+

      // std::cout << "dataPlot(k, 27)" << dataPlot(k, 27)  << std::endl;
      // std::cout << "dataPlot(k, 31)" << dataPlot(k, 31)  << std::endl;

      // siconos::algebra::print(*q);
      // siconos::algebra::print(*v);

      s->processEvents();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("SliderCrankD1MinusLinearOSIVelocityLevel.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-9;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "SliderCrankD1MinusLinearOSIVelocityLevel.ref", eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
